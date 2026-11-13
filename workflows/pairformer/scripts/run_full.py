#!/usr/bin/env python3
"""One-model, sequential inference with a larger-input smoke gate and per-MSA checkpoints."""
import argparse
import importlib.metadata as metadata
import io
import json
import os
import platform
import subprocess
import sys
import time
import traceback
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent/'pilot_lib'))
from common import require, sha256, fingerprint, read_fasta, atomic_bytes, write_json
from run_inference import validate_matrix, validated_checkpoint


def save_matrix(folder, matrix, record):
    folder.mkdir(parents=True, exist_ok=True)
    buf=io.BytesIO(); np.save(buf, matrix, allow_pickle=False)
    atomic_bytes(folder/'contacts.npy', buf.getvalue())
    csv=io.StringIO(); np.savetxt(csv, matrix, delimiter=',')
    atomic_bytes(folder/'full_contacts.csv', csv.getvalue().encode())
    write_json(folder/'COMPLETE.json', dict(record, contacts_sha256=sha256(folder/'contacts.npy')))
    if (folder/'FAILED.json').exists():
        (folder/'FAILED.json').unlink()


def run(root, pilot):
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Need a Slurm allocation')
    require(not platform.node().split('.')[0].startswith('biowulf'), 'Never run on login node')
    import torch
    require(torch.cuda.is_available() and torch.cuda.is_bf16_supported(), 'Need a bfloat16-capable GPU')
    torch.set_num_threads(int(os.environ.get('SLURM_CPUS_PER_TASK','1')))
    torch.set_num_interop_threads(1)
    torch.backends.cuda.matmul.allow_tf32=False
    torch.backends.cudnn.allow_tf32=False
    torch.backends.cudnn.benchmark=False
    import MSA_Pairformer
    from MSA_Pairformer.model import MSAPairformer
    from MSA_Pairformer.dataset import aa2tok_d, prepare_msa_masks

    config=json.loads((root/'model_config.json').read_text())
    manifest=json.loads((root/'inputs/manifest.json').read_text())
    require(len(manifest['inputs'])==396, 'Incomplete manifest')
    require(config['max_depth']==512 and config['query_biasing'] is False, 'Wrong target configuration')
    require(manifest['selection']['cap']==512, 'Wrong prepared depth cap')
    require(all(s['selected_depth']==min(s['depth'],512) for s in manifest['inputs']), 'Wrong selected depth')
    pilot_run=json.loads((pilot/'results/all_rows_cb/RUN.json').read_text())
    direct=json.loads(metadata.distribution('msa-pairformer').read_text('direct_url.json') or '{}')
    require(direct.get('vcs_info',{}).get('commit_id')==config['source_commit'], 'Wrong source commit')
    freeze=subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True)
    require(freeze==pilot_run['pip_freeze'], 'Environment differs from successful pilot')
    package=Path(MSA_Pairformer.__file__).resolve().parent
    package_hashes={str(p.relative_to(package)):sha256(p) for p in package.rglob('*.py')}
    require(package_hashes==pilot_run['package_hashes'], 'Installed source changed since pilot')
    snapshot=pilot/'weights'/config['weights_repo_directory']/'snapshots'/config['weights_revision']
    weights={name:sha256(snapshot/name) for name in config['weights_sha256']}
    require(weights==config['weights_sha256']==pilot_run['weights'], 'Checkpoint mismatch')
    require(pilot_run['precision']=='bfloat16' and pilot_run['head']=='cb' and
            pilot_run['query_biasing'] is True and pilot_run['chain_breaks'] is None, 'Pilot inference settings differ')
    settings=dict(source_commit=config['source_commit'],weights=weights,precision='bfloat16',head='cb',
                  seed=config['seed'],max_depth=config['max_depth'],query_biasing=False,sequence_weighting='uniform',chain_breaks=None,package_hashes=package_hashes,
                  pip_freeze=freeze,manifest_fingerprint=manifest['fingerprint'],
                  runner_sha256=sha256(Path(__file__)),common_sha256=sha256(Path(__file__).resolve().parent/'pilot_lib/common.py'),
                  validator_sha256=sha256(Path(__file__).resolve().parent/'pilot_lib/run_inference.py'))
    key=fingerprint(settings)
    out=root/'results'; out.mkdir(exist_ok=True)
    if (out/'RUN.json').exists():
        require(json.loads((out/'RUN.json').read_text())['run_fingerprint']==key,
                'Experiment changed; use a separate results directory')
    else:
        write_json(out/'RUN.json',dict(settings,run_fingerprint=key))
    pending=[]; reused=0; existing=0
    for item in manifest['inputs']:
        require(sha256(root/item['frozen_path'])==item['selected_sha256'], 'Frozen input changed')
        folder=out/item['id']
        if validated_checkpoint(folder,key,item['length']):
            existing+=1
        else:
            pending.append(item)
    print(f'Inventory: validated existing={existing}, reused pilot={reused}, pending={len(pending)}',flush=True)
    if not pending:
        return
    np.random.seed(config['seed']); torch.manual_seed(config['seed']); torch.cuda.manual_seed_all(config['seed'])
    start=time.monotonic()
    state=torch.load(snapshot/'model_cuex.bin',map_location='cpu',weights_only=True)
    state.update(torch.load(snapshot/'confind_contact.bin',map_location='cpu',weights_only=True))
    state.update(torch.load(snapshot/'contact.bin',map_location='cpu',weights_only=True))
    model=MSAPairformer(); model.load_state_dict(state,strict=True); del state
    device=torch.device('cuda:0'); model=model.to(device).eval(); model.turn_on_query_biasing()
    torch.cuda.synchronize(); load_seconds=time.monotonic()-start

    weight_audit=[]
    def predict(item, audit_weights=False):
        rows=read_fasta(root/item['frozen_path'])
        require(len(rows)==item['selected_depth'] and rows[0][1]==item['query_sequence'], 'Bad selected MSA')
        tokens=torch.tensor([[aa2tok_d[a] for a in s] for _,s in rows],dtype=torch.long).unsqueeze(0)
        masks=prepare_msa_masks(tokens,device=device)
        msa=torch.nn.functional.one_hot(tokens,num_classes=len(aa2tok_d)).to(device,dtype=torch.bfloat16)
        torch.cuda.reset_peak_memory_stats(); torch.cuda.synchronize(); began=time.monotonic()
        with torch.inference_mode(),torch.amp.autocast(device_type='cuda',dtype=torch.bfloat16):
            result=model.predict_cb_contacts(msa=msa,mask=masks[0],msa_mask=masks[1],full_mask=masks[2],
                    pairwise_mask=masks[3],complex_chain_break_indices=None,return_seq_weights=audit_weights)
        if audit_weights:
            audits=[]
            for layer, weights in result['seq_weights_list_d'].items():
                error=float((weights.float()-1.0/len(rows)).abs().max())
                require(error<1e-7, f'Nonuniform sequence weighting in {layer}: {error}')
                audits.append(dict(layer=layer,shape=list(weights.shape),max_abs_error_from_uniform=error))
            require(len(audits)==model.contact_layer+1, 'Missing weighting audit layers')
            weight_audit.append(dict(id=item['id'],depth=len(rows),layers=audits))
        matrix=result['predicted_cb_contacts'][0].float().cpu().numpy()
        torch.cuda.synchronize(); seconds=time.monotonic()-began
        validate_matrix(matrix,item['length'])
        stats=dict(inference_seconds=seconds,peak_gpu_allocated_bytes=torch.cuda.max_memory_allocated(),
                   peak_gpu_reserved_bytes=torch.cuda.max_memory_reserved())
        require(stats['peak_gpu_reserved_bytes'] < .85*torch.cuda.get_device_properties(0).total_memory,
                'Insufficient GPU memory headroom; stop before scaling')
        return matrix,stats

    # One control replay checks the new driver against the preserved pilot, not new biological work.
    control=next(s for s in manifest['inputs'] if s['id']=='proper_1gb1')
    control_matrix,control_stats=predict(control)
    original=np.load(pilot/'results/all_rows_cb/proper_1gb1/contacts.npy',allow_pickle=False)
    control_delta=float(abs(original-control_matrix).max())
    require(control_delta<=1e-3, 'Control replay differs from pilot by >1e-3')
    # The old biased replay is a technical control only; none of its maps enter this run.
    model.turn_off_query_biasing()
    from MSA_Pairformer.outer_product import OuterProduct
    flags=[m.opm.seq_attn for m in model.modules() if isinstance(m,OuterProduct)]
    require(len(flags)==len(model.core_stack.layers) and not any(flags), 'Query bias was not disabled in every layer')
    smoke_outputs={}
    for name in ('proper_1gb1','proper_2mm4','swap_1gb1_2mm4'):
        probe=next(s for s in manifest['inputs'] if s['id']==name)
        smoke_outputs[name]=predict(probe,audit_weights=True)
    write_json(out/'WEIGHTING_AUDIT.json',dict(query_biasing=False,all_layer_flags=flags,
               sequence_weighting='uniform',control_biased_max_abs_difference=control_delta,
               audited_inputs=weight_audit))
    # Largest selected workloads first; their validated outputs count toward production.
    ordered=sorted(pending,key=lambda s:(s['length']**3*s['selected_depth'],s['id']),reverse=True)
    benchmark=[]; start=time.monotonic()
    for number,item in enumerate(ordered,1):
        folder=out/item['id']
        try:
            matrix,stats=smoke_outputs.pop(item['id']) if item['id'] in smoke_outputs else predict(item, audit_weights=(number==1))
            info=dict(stats,run_fingerprint=key,selected_sha256=item['selected_sha256'],
                      source_sha256=item['source_sha256'],selected_depth=item['selected_depth'],length=item['length'],
                      kind=item['kind'],query_biasing=False,max_depth=512,sequence_weighting='uniform',job_id=os.environ['SLURM_JOB_ID'],host=platform.node(),
                      gpu=torch.cuda.get_device_name(0),model_load_seconds=load_seconds)
            save_matrix(folder,matrix,info)
            if number<=6:
                benchmark.append(dict(id=item['id'],**stats))
            if number==min(6,len(ordered)):
                write_json(out/'BENCHMARK.json',dict(control_max_abs_difference=control_delta,control=control_stats,
                           largest_inputs=benchmark,model_load_seconds=load_seconds,
                           remaining=len(ordered)-number,
                           estimated_remaining_seconds=(len(ordered)-number)*max(b['inference_seconds'] for b in benchmark),
                           scheduling='All remaining inputs grouped in this one process; no per-MSA jobs'))
                print('Larger-input smoke gate passed; proceeding in the same allocation',flush=True)
            print(f"{number}/{len(ordered)} {item['id']} L={item['length']} N={item['selected_depth']} {stats['inference_seconds']:.3f}s",flush=True)
        except Exception:
            error=traceback.format_exc(); write_json(folder/'FAILED.json',dict(run_fingerprint=key,error=error))
            raise
    write_json(out/'WEIGHTING_AUDIT.json',dict(query_biasing=False,all_layer_flags=flags,sequence_weighting='uniform',control_biased_max_abs_difference=control_delta,audited_inputs=weight_audit))
    for item in manifest['inputs']:
        require(validated_checkpoint(out/item['id'],key,item['length']), 'Final output validation failed')
    write_json(out/'INFERENCE_COMPLETE.json',dict(run_fingerprint=key,validated=396,new=len(ordered),
               reused_pilot=reused,previously_complete=existing,seconds=time.monotonic()-start,
               job_id=os.environ['SLURM_JOB_ID'],control_max_abs_difference=control_delta))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--root',type=Path,required=True); p.add_argument('--pilot',type=Path,required=True)
    a=p.parse_args(); run(a.root.resolve(),a.pilot.resolve())
