#!/usr/bin/env python3
"""Load one pinned local checkpoint and run all three monomer contact predictions."""
import argparse
import importlib.metadata as metadata
import io
import json
import os
import platform
import resource
import subprocess
import sys
import time
import traceback
from pathlib import Path
import numpy as np
from common import atomic_bytes, fingerprint, read_fasta, require, sha256, write_json


def validate_matrix(matrix, length):
    require(matrix.shape == (length, length), f"Unexpected contact shape: {matrix.shape}")
    require(np.isfinite(matrix).all(), "Nonfinite contact scores")
    require((matrix >= 0).all() and (matrix <= 1).all(), "Scores outside [0,1]")
    require(np.allclose(matrix, matrix.T, rtol=0, atol=1e-6), "Asymmetric contact matrix")


def validated_checkpoint(folder, expected_fingerprint, length):
    marker = folder / "COMPLETE.json"
    if not marker.exists():
        return False
    saved = json.loads(marker.read_text())
    require(saved["run_fingerprint"] == expected_fingerprint, f"Changed settings at {folder}")
    require(saved["contacts_sha256"] == sha256(folder / "contacts.npy"), f"Corrupt contacts at {folder}")
    validate_matrix(np.load(folder / "contacts.npy", allow_pickle=False), length)
    return True


def run(root):
    require(bool(os.environ.get("SLURM_JOB_ID")), "Inference requires a Slurm compute allocation")
    require(not platform.node().split('.')[0].startswith("biowulf"), "Do not infer on a login node")
    import torch
    require(torch.cuda.is_available(), "CUDA GPU unavailable; no CPU fallback")
    require(torch.cuda.is_bf16_supported(), "Pilot requires a GPU with bfloat16 support")
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    torch.set_num_threads(threads)
    torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cudnn.benchmark = False
    # Import only after the CUDA check; upstream initializes GPU kernels at import.
    import MSA_Pairformer
    from MSA_Pairformer.model import MSAPairformer
    from MSA_Pairformer.dataset import aa2tok_d, prepare_msa_masks

    manifest = json.loads((root / "inputs/manifest.json").read_text())
    config = manifest["config"]
    config_now = json.loads((root / "config.json").read_text())
    require(config == config_now, "Configuration differs from frozen manifest")
    ready = json.loads((root / "env/READY.json").read_text())
    require(ready["config_sha256"] == sha256(root / "config.json"), "Environment/config mismatch")
    require(ready["constraints_sha256"] == sha256(root / "constraints.txt"), "Environment/constraints mismatch")
    direct_url = json.loads(metadata.distribution("msa-pairformer").read_text("direct_url.json") or "{}")
    require(direct_url.get("vcs_info", {}).get("commit_id") == config["source_commit"],
            "Installed Pairformer is not the pinned Git commit; run setup_env.sh")
    snapshot = root / "weights" / config["weights_repo_directory"] / "snapshots" / config["weights_revision"]
    filenames = ("model_cuex.bin", "contact.bin", "confind_contact.bin")
    for name in filenames:
        require((snapshot / name).is_file(), f"Missing checkpoint file: {snapshot / name}")
    hashes = {name: sha256(snapshot / name) for name in filenames}
    require(hashes == config["weights_sha256"], "Checkpoint hashes differ from inspected shared snapshot")
    package_root = Path(MSA_Pairformer.__file__).resolve().parent
    package_hashes = {str(p.relative_to(package_root)): sha256(p) for p in package_root.rglob("*.py")}
    installed = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
    scripts = {p.name: sha256(p) for p in (root / "common.py", root / "run_inference.py", root / "prepare_inputs.py")}
    provenance = dict(manifest_fingerprint=manifest["fingerprint"], weights=hashes,
                      weights_revision=config["weights_revision"], source_commit=config["source_commit"],
                      package_hashes=package_hashes, scripts=scripts, pip_freeze=installed,
                      python=sys.version, torch=torch.__version__, cuda=torch.version.cuda,
                      gpu=torch.cuda.get_device_name(0), precision="bfloat16", seed=config["seed"],
                      threads=threads, head="cb", query_biasing=True, chain_breaks=None)
    run_key = fingerprint(provenance)
    output = root / "results/all_rows_cb"
    output.mkdir(parents=True, exist_ok=True)
    if (output / "RUN.json").exists():
        old = json.loads((output / "RUN.json").read_text())
        require(old["run_fingerprint"] == run_key,
                "Inputs/software/hardware/settings changed; preserve results and use a new pilot directory")
    else:
        write_json(output / "RUN.json", dict(provenance, run_fingerprint=run_key))
    pending = []
    for spec in manifest["inputs"]:
        require(sha256(root / spec["frozen_path"]) == spec["frozen_sha256"], "Frozen FASTA changed")
        if validated_checkpoint(output / spec["id"], run_key, spec["expected_length"]):
            print("Validated; skipping:", spec["id"], flush=True)
        else:
            pending.append(spec)
    if not pending:
        return
    np.random.seed(config["seed"])
    torch.manual_seed(config["seed"])
    torch.cuda.manual_seed_all(config["seed"])
    started = time.monotonic()
    # Mirror upstream's checkpoint merge, but bind to this exact local snapshot.
    # No network calls and no glob across potentially different HF revisions.
    state = torch.load(snapshot / "model_cuex.bin", map_location="cpu", weights_only=True)
    state.update(torch.load(snapshot / "confind_contact.bin", map_location="cpu", weights_only=True))
    state.update(torch.load(snapshot / "contact.bin", map_location="cpu", weights_only=True))
    model = MSAPairformer()
    model.load_state_dict(state, strict=True)
    del state
    device = torch.device("cuda:0")
    model = model.to(device).eval()
    model.turn_on_query_biasing()
    torch.cuda.synchronize()
    load_seconds = time.monotonic() - started
    print(f"Loaded checkpoint once in {load_seconds:.2f}s", flush=True)
    failures = []
    first_id = manifest["inputs"][0]["id"]
    for spec in pending:
        folder = output / spec["id"]
        folder.mkdir(parents=True, exist_ok=True)
        try:
            rows = read_fasta(root / spec["frozen_path"])
            require(len(rows) == spec["expected_depth"] and rows[0][1] == spec["query"], "Input mismatch")
            tokens = torch.tensor([[aa2tok_d[a] for a in s] for _, s in rows], dtype=torch.long).unsqueeze(0)
            masks = prepare_msa_masks(tokens, device=device)
            msa = torch.nn.functional.one_hot(tokens, num_classes=len(aa2tok_d)).to(device, dtype=torch.bfloat16)

            def predict():
                with torch.inference_mode(), torch.amp.autocast(device_type="cuda", dtype=torch.bfloat16):
                    result = model.predict_cb_contacts(
                        msa=msa, mask=masks[0], msa_mask=masks[1], full_mask=masks[2],
                        pairwise_mask=masks[3], complex_chain_break_indices=None,
                        return_seq_weights=False)
                return result["predicted_cb_contacts"][0].float().cpu().numpy()

            torch.cuda.reset_peak_memory_stats()
            torch.cuda.synchronize()
            start = time.monotonic()
            contacts = predict()
            torch.cuda.synchronize()
            elapsed = time.monotonic() - start
            validate_matrix(contacts, spec["expected_length"])
            repeat_difference = None
            if spec["id"] == first_id:
                repeated = predict()
                validate_matrix(repeated, spec["expected_length"])
                repeat_difference = float(np.abs(contacts - repeated).max())
                require(repeat_difference <= 1e-3, "Repeated inference differs by >1e-3")
            buffer = io.BytesIO()
            np.save(buffer, contacts, allow_pickle=False)
            atomic_bytes(folder / "contacts.npy", buffer.getvalue())
            csv = io.StringIO()
            np.savetxt(csv, contacts, delimiter=",")
            atomic_bytes(folder / "full_contacts.csv", csv.getvalue().encode())
            details = dict(run_fingerprint=run_key, contacts_sha256=sha256(folder / "contacts.npy"),
                           input_sha256=spec["frozen_sha256"], depth=len(rows), length=contacts.shape[0],
                           inference_seconds=elapsed, model_load_seconds=load_seconds,
                           peak_gpu_allocated_bytes=torch.cuda.max_memory_allocated(),
                           peak_gpu_reserved_bytes=torch.cuda.max_memory_reserved(),
                           process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                           repeat_max_abs_difference=repeat_difference, slurm_job_id=os.environ["SLURM_JOB_ID"],
                           host=platform.node())
            write_json(folder / "COMPLETE.json", details)
            if (folder / "FAILED.json").exists():
                (folder / "FAILED.json").unlink()
            print(f"{spec['id']}: {contacts.shape}, {elapsed:.2f}s; checkpoint saved", flush=True)
            del tokens, masks, msa, contacts
        except Exception:
            error = traceback.format_exc()
            write_json(folder / "FAILED.json", dict(run_fingerprint=run_key, error=error))
            print(error, file=sys.stderr, flush=True)
            failures.append(spec["id"])
            # CUDA failures can poison the context. Stop safely; completed items persist.
            break
    require(not failures, f"Pilot failed: {failures}. Inspect FAILED.json and logs before retrying.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    run(parser.parse_args().root.resolve())
