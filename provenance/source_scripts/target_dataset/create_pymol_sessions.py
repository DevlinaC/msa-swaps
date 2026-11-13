"""Create editable PyMOL galleries for the 60-, 100-, and 150-aa target groups.

Run with the local PyMOL/Biopython/NumPy environment:
    python create_pymol_sessions.py
    python create_pymol_sessions.py --groups 60 --width 2400

Outputs: one self-contained folds_<size>.pse session and ray-traced PNG per group,
filtered first-model PDBs, a manifest, and small .pml launchers. Each session has
an overview scene plus an individual scene per target. Structures are arranged
as separate objects; there is NO structural alignment or morphing. A proper
rigid PCA rotation provides a reproducible starting orientation, which can be
adjusted interactively. Original coordinate files are never changed.

1BOR uses the downloaded CIF and its deposited secondary-structure annotations.
Its C-alpha coordinates are checked against the scoring PDB before replacement.
Other structures use PyMOL DSS.

Residues match the TM-score coordinate loader (same first model, mapped chain,
MODRES policy, and observed C-alpha residues). Full atoms of those residues are
retained, using Biopython's selected alternate conformer. Incomplete residues
without CA and separate ligands are excluded. No missing residues are rebuilt.

Colors represent secondary structure (CIF for 1BOR, PyMOL DSS otherwise): helix orange, strand blue,
loop gray. Color does not denote CATH fold class or prove fold uniqueness. Labels
show PDB ID and target FASTA length, not resolved-coordinate length. The manifest
records both. PCA can obscure some features; the stored individual scenes allow
manual rotation before final publication. Rendering is bounded to one CPU.
On Biowulf, execute inside a Slurm compute allocation, not the login node.
"""

import argparse
import csv
import json
import math
from pathlib import Path

from create_tmscore_matrix import BASE_DIR, checksum, load_coordinates
from create_table_S1 import EXPECTED_COUNTS, read_fasta, read_csv, unique_lookup
import numpy as np
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, Select
import pymol2


class TargetAtoms(Select):
    """Retain the selected model/chain/residues and one alternate per atom."""
    def __init__(self, chain, residue_ids):
        self.chain, self.residue_ids = chain, set(residue_ids)

    def accept_model(self, model):
        return model.id == 0

    def accept_chain(self, chain):
        return chain.id == self.chain

    def accept_residue(self, residue):
        return residue.id in self.residue_ids

    def accept_atom(self, atom):
        selected = atom.get_parent()[atom.name]
        return atom.altloc == selected.get_altloc()


def rotation(coords):
    """Right-handed PCA basis, longest axis vertical; never reflect a fold."""
    center = coords.mean(axis=0)
    _, _, axes = np.linalg.svd(coords-center, full_matrices=False)
    for vector in axes:
        if vector[np.argmax(np.abs(vector))] < 0:
            vector *= -1
    basis = axes[[1, 0, 2]].T
    if np.linalg.det(basis) < 0:
        basis[:, 2] *= -1
    return center, basis


def prepare_targets(args):
    mappings = unique_lookup(read_csv(args.chain_map, ['pdb_id', 'chain_id']), args.chain_map)
    prepared = []
    pdb_dir = args.output_dir/'structures'
    pdb_dir.mkdir(parents=True, exist_ok=True)
    for size in EXPECTED_COUNTS:
        records = list(read_fasta(BASE_DIR/f'targets_{size}.fa'))
        if len(records) != EXPECTED_COUNTS[size]:
            raise ValueError(f'Unexpected target count in {size}-aa group')
        for pid, header_chain, sequence in sorted(records, key=lambda r: (len(r[2]), r[0])):
            chain = mappings[pid]['chain_id']
            if header_chain and header_chain != chain:
                raise ValueError(f'{pid}: conflicting chain IDs')
            source = args.pdb_dir/f'{pid.lower()}.pdb'
            ca, observed, ids, _ = load_coordinates(source, chain, sequence)
            cif_source = args.bor_cif if pid == '1BOR' else None
            if cif_source:
                structure = MMCIFParser(QUIET=True).get_structure(pid, cif_source)
                cif_residues = [r for r in structure[0][chain] if r.id[0] == ' ' and 'CA' in r]
                cif_ca = np.asarray([r['CA'].coord for r in cif_residues])
                # This replacement is for display metadata, not a new analysis.
                if cif_ca.shape != ca.shape or not np.allclose(cif_ca, ca, atol=.001, rtol=0):
                    raise ValueError('1BOR CIF differs from scoring coordinates; review the replacement before plotting')
            else:
                structure = PDBParser(PERMISSIVE=False, QUIET=True).get_structure(pid, source)
            # Use complete residue IDs so a separate HET ligand sharing a number
            # cannot accidentally join the protein selection.
            ca_ids = set(ids)
            modified = set()
            for line in source.read_text().splitlines():
                if line.startswith('MODRES') and line[16] == chain:
                    modified.add((int(line[18:22]), line[22], line[12:15].strip()))
            kept = [r.id for r in structure[0][chain]
                    if f'{r.id[1]}{r.id[2].strip()}' in ca_ids and 'CA' in r
                    and (r.id[0] == ' ' or (r.id[1], r.id[2], r.resname) in modified)]
            output = pdb_dir/f'{pid}_{chain}.pdb'
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(output), TargetAtoms(chain, kept))
            center, basis = rotation(ca)
            oriented = (ca-center)@basis
            prepared.append(dict(pdb_id=pid, chain_id=chain, length_class=size,
                                 fasta_length=len(sequence), n_ca=len(ca),
                                 source=str((cif_source or source).resolve()), source_sha256=checksum(cif_source or source),
                                 secondary_structure='CIF annotations' if cif_source else 'PyMOL DSS',
                                 display_cif=str(cif_source.resolve()) if cif_source else None,
                                 prepared_pdb=str(output.resolve()), center=center.tolist(),
                                 rotation=basis.tolist(), extent=np.ptp(oriented,axis=0).tolist()))
    return prepared


def build_session(records, args, cell_width):
    size = records[0]['length_class']
    if size == 60:
        cell_width = max(r['extent'][0] for r in records)+9
    columns = 4
    nrows = math.ceil(len(records)/columns)
    # Compact each row around its actual tallest member, not the longest protein
    # anywhere in the dataset. No coordinate scaling or geometry changes.
    row_heights = [max(r['extent'][1] for r in records[j:j+columns])+17
                   for j in range(0, len(records), columns)]
    row_tops = -np.r_[0, np.cumsum(row_heights[:-1])]
    body_height = sum(row_heights)
    width_world = columns*cell_width
    height_world = body_height + 26
    height = round(args.width*height_world/width_world)
    stem = args.output_dir/f'folds_{size}'
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        cmd.set('max_threads', 1)
        cmd.bg_color('white')
        cmd.set('orthoscopic', 1)
        cmd.set('depth_cue', 0)
        cmd.set('ray_shadows', 0)
        cmd.set('ray_opaque_background', 1)
        cmd.set('antialias', 2)
        cmd.set('specular', .15)
        cmd.set('cartoon_fancy_helices', 1)
        cmd.set('cartoon_flat_sheets', 1)
        cmd.set('cartoon_loop_radius', .28)
        cmd.set('label_font_id', 7)  # Sans Bold
        cmd.set('label_size', -6.2 if size == 60 else -4.8)  # World units: stable size across exports.
        cmd.set('label_color', 'black')
        cmd.set('label_outline_color', 'white')
        cmd.set('label_position', [0, 0, 0])
        cmd.set('label_shadow_mode', 0)
        cmd.set_color('fold_helix', [0.90, .48, .05])
        cmd.set_color('fold_strand', [0, .45, .70])
        cmd.set_color('fold_loop', [.48, .48, .48])
        object_names = []
        for i, record in enumerate(records):
            obj = f"p_{record['pdb_id']}"
            object_names.append(obj)
            if record['display_cif']:
                cmd.load(record['display_cif'], 'source_cif')
                cmd.create(obj, f"source_cif and chain {record['chain_id']} and polymer.protein", 1, 1)
                cmd.delete('source_cif')
            else:
                cmd.load(record['prepared_pdb'], obj)
            if cmd.count_states(obj) != 1 or cmd.count_atoms(f'{obj} and name CA') != record['n_ca']:
                raise ValueError(f'{obj}: exported coordinates differ from scoring residues')
            xyz = cmd.get_coords(obj, state=1)
            row = i // columns
            # Center an incomplete final row instead of leaving unused edge panels.
            row_count = min(columns, len(records)-row*columns)
            offset = np.array([(i % columns + .5 + (columns-row_count)/2)*cell_width,
                               row_tops[row]-row_heights[row]/2-3, 0.])
            coords = (xyz-np.array(record['center']))@np.array(record['rotation'])
            # Bounding-box centering improves whitespace without distorting shape.
            correction = (coords.min(axis=0)+coords.max(axis=0))/2
            coords = coords-correction+offset
            record['bounding_box_center_correction'] = correction.tolist()
            cmd.load_coords(coords, obj, state=1)
            record['gallery_translation'] = offset.tolist()
            cmd.hide('everything', obj)
            if not record['display_cif']:
                cmd.dss(obj)
            cmd.show('cartoon', obj)
            cmd.color('fold_loop', obj)
            cmd.color('fold_helix', f'{obj} and ss H')
            cmd.color('fold_strand', f'{obj} and ss S')
            label = record['pdb_id'] if size == 60 else f"{record['pdb_id']}  ({record['fasta_length']} aa)"
            cmd.pseudoatom(f"label_{record['pdb_id']}", pos=[offset[0], row_tops[row]-3, 0], label=label)
            cmd.hide('nonbonded', f"label_{record['pdb_id']}")
        title = f'~{size} aa targets | {len(records)} proteins'
        cmd.pseudoatom('gallery_title', pos=[width_world/2, 12, 0], label=title)
        cmd.set('label_size', -5.8, 'gallery_title')
        cmd.pseudoatom('gallery_legend', pos=[width_world/2, -body_height-6, 0],
                       label='Helix: orange   Strand: blue   Loop: gray')
        cmd.set('label_size', -3.7, 'gallery_legend')
        cmd.hide('nonbonded', 'gallery_title or gallery_legend')
        cmd.viewport(args.width, height)
        # Set orthographic framing explicitly: spherical zoom adds excessive
        # margins for these tall galleries. Keep identical physical pixel scale.
        distance = (height_world + 6) / (2*math.tan(math.radians(10)))
        center_y = (18 - body_height - 10)/2
        cmd.set_view((1,0,0, 0,1,0, 0,0,1, 0,0,-distance,
                      width_world/2,center_y,0, distance-150,distance+150,20))
        cmd.scene('overview', 'store')
        cmd.png(str(stem.with_suffix('.png')), width=args.width, height=height, dpi=300, ray=1, quiet=1)
        for record, obj in zip(records, object_names):
            cmd.disable('all')
            cmd.enable(obj)
            # An independent camera fills each publication panel without
            # scaling molecular coordinates or changing the cartoon styling.
            panel_dir = args.output_dir/f'panels_{size}'
            panel_dir.mkdir(exist_ok=True)
            cmd.viewport(750, 750)
            cmd.zoom(obj, buffer=1, complete=0)
            cmd.png(str(panel_dir/f"{record['pdb_id']}.png"), width=750, height=750,
                    dpi=300, ray=1, quiet=1)
            cmd.enable(f"label_{record['pdb_id']}")
            cmd.zoom(f"{obj} or label_{record['pdb_id']}", buffer=6, complete=1)
            cmd.scene(record['pdb_id'], 'store')
        cmd.viewport(args.width, height)
        cmd.scene('overview', 'recall')
        cmd.set('scene_buttons', 1)
        cmd.save(str(stem.with_suffix('.pse')))
        version = cmd.get_version()[0]
    # Reopen the actual session to check embedded objects, chain/state counts,
    # and scene persistence rather than trusting save success alone.
    with pymol2.PyMOL() as check:
        cmd = check.cmd
        cmd.load(str(stem.with_suffix('.pse')))
        if len(cmd.get_scene_list()) != len(records)+1:
            raise ValueError(f'{stem}: missing saved scenes')
        for record, obj in zip(records, object_names):
            if cmd.count_atoms(f'{obj} and name CA') != record['n_ca'] or cmd.count_states(obj) != 1:
                raise ValueError(f'{stem}: saved structure validation failed')
            if cmd.get_chains(obj) != [record['chain_id']]:
                raise ValueError(f'{stem}: unexpected chain')
    session_path = str(stem.with_suffix('.pse').resolve())
    stem.with_suffix('.pml').write_text(
        '# Open this launcher from any working directory.\npython\n'
        f'cmd.load({session_path!r})\n'
        "cmd.scene('overview', 'recall')\npython end\n")
    compose_panels(records, args.output_dir)
    print(f'Created and verified {stem.name}.pse ({len(records)} proteins, {len(records)+1} scenes)', flush=True)
    return version


def compose_panels(records, output):
    """Publication montage with independent camera zoom; coordinates unchanged.

    Saved separately from the session overview. Same helix/strand styling and
    PCA orientation as the session; panels differ only in camera magnification.
    """
    import os
    import tempfile
    os.environ.setdefault('MPLCONFIGDIR', tempfile.mkdtemp(prefix='pymol_panels_', dir=os.environ.get('TMPDIR')))
    os.environ.setdefault('XDG_CACHE_HOME', os.environ['MPLCONFIGDIR'])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.weight':'bold',
                         'pdf.fonttype':42, 'axes.titleweight':'bold'})
    size = records[0]['length_class']
    nrows = math.ceil(len(records)/4)
    figure_height = nrows*3.5+1.5
    fig = plt.figure(figsize=(14, figure_height), facecolor='white')
    grid = fig.add_gridspec(nrows, 8, left=.015, right=.985, bottom=.95/figure_height,
                           top=1-1.25/figure_height,
                           wspace=.02, hspace=.12)
    for i, record in enumerate(records):
        row = i//4
        row_count = min(4, len(records)-row*4)
        col = (i%4)*2 + (4-row_count)
        ax = fig.add_subplot(grid[row,col:col+2])
        pixels = plt.imread(output/f'panels_{size}'/f"{record['pdb_id']}.png")
        foreground = np.any(pixels[..., :3] < .97, axis=2)
        ys, xs = np.where(foreground)
        if not len(xs):
            raise ValueError(f"Blank panel for {record['pdb_id']}")
        cropped = pixels[max(0,ys.min()-8):ys.max()+9, max(0,xs.min()-8):xs.max()+9]
        # Square white padding preserves projection aspect ratio and consistent
        # label positions while removing ray-tracer framing margins.
        side = max(cropped.shape[:2])+24
        canvas = np.ones((side,side,pixels.shape[2]), dtype=pixels.dtype)
        y, x = (side-cropped.shape[0])//2, (side-cropped.shape[1])//2
        canvas[y:y+cropped.shape[0],x:x+cropped.shape[1]] = cropped
        ax.imshow(canvas)
        ax.set_axis_off()
        ax.set_title(record['pdb_id'], fontsize=22, fontweight='bold', pad=0)
    fig.suptitle(f'~{size} aa targets: {len(records)} protein structures', fontsize=27,
                 fontweight='bold', y=1-.18/figure_height)
    fig.legend(handles=[Patch(color='#E67A0D',label='Helix'),Patch(color='#0073B3',label='Strand'),
                        Patch(color='#7A7A7A',label='Loop')],
               loc='lower center',bbox_to_anchor=(.5,.38/figure_height),ncol=3,frameon=False,
               prop={'size':18,'weight':'bold'})
    fig.text(.5,.16/figure_height,'Each protein is independently zoomed; panel scales differ.',
             ha='center',fontsize=14,fontweight='bold')
    for extension in ('png','pdf'):
        fig.savefig(output/f'folds_{size}_panels.{extension}',dpi=300,facecolor='white')
    plt.close(fig)


def panels_from_session(records, output, reuse=False):
    """Render independently zoomed panels from an existing verified session.

    Reuse is only appropriate for text/layout-only changes. Rebuild sessions
    normally after structural, orientation, or cartoon-style changes.
    """
    size = records[0]['length_class']
    directory = output/f'panels_{size}'
    directory.mkdir(exist_ok=True)
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        cmd.load(str(output/f'folds_{size}.pse'))
        cmd.set('max_threads',1)
        cmd.scene('overview','recall')
        cmd.viewport(750,750)
        for record in records:
            path = directory/f"{record['pdb_id']}.png"
            if reuse and path.exists():
                continue
            obj = f"p_{record['pdb_id']}"
            if cmd.count_atoms(f'{obj} and name CA') != record['n_ca']:
                raise ValueError(f'{obj}: session/input mismatch')
            cmd.disable('all')
            cmd.enable(obj)
            cmd.zoom(obj,buffer=1,complete=0)
            cmd.png(str(path),width=750,height=750,dpi=300,ray=1,quiet=1)
    compose_panels(records,output)
    print(f'Updated folds_{size}_panels PNG and PDF',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--pdb-dir', type=Path, default=BASE_DIR/'original_pdbs')
    parser.add_argument('--chain-map', type=Path, default=BASE_DIR/'target_chains.csv')
    parser.add_argument('--output-dir', type=Path, default=BASE_DIR/'outputs/pymol')
    parser.add_argument('--bor-cif', type=Path, default=BASE_DIR/'outputs/pymol/structures/1BOR-assembly1.cif',
                        help='Downloaded 1BOR CIF; preserve deposited secondary-structure annotations')
    parser.add_argument('--groups', type=int, nargs='+', choices=(60,100,150), default=[60,100,150])
    parser.add_argument('--width', type=int, default=3000, help='Gallery PNG width in pixels')
    parser.add_argument('--panels-only', action='store_true',
                        help='Use existing sessions to update publication montages without rebuilding galleries')
    parser.add_argument('--reuse-panel-renders', action='store_true',
                        help='With --panels-only, reuse existing per-protein PNGs for text/layout edits')
    args = parser.parse_args()
    if args.width < 800:
        parser.error('--width must be at least 800')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.panels_only:
        manifest = json.loads((args.output_dir/'manifest.json').read_text())
        for size in args.groups:
            records = [r for r in manifest['targets'] if r['length_class']==size]
            if len(records) != EXPECTED_COUNTS[size]:
                parser.error(f'Incomplete manifest for {size}-aa group')
            panels_from_session(records,args.output_dir,args.reuse_panel_renders)
        panel_manifest = dict(script_sha256=checksum(Path(__file__)), groups=args.groups,
                              pdb_label_points=22, title_points=27, independent_zoom=True,
                              session_hashes={str(size):checksum(args.output_dir/f'folds_{size}.pse') for size in args.groups})
        (args.output_dir/'panels_manifest.json').write_text(json.dumps(panel_manifest,indent=2))
        return
    records = prepare_targets(args)
    # One scale across all groups: rotations/translations only, no resizing.
    cell_width = max(r['extent'][0] for r in records)+12
    manifest_path = args.output_dir/'manifest.json'
    previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    versions = previous.get('pymol_versions', {}).copy()
    for size in args.groups:
        versions[str(size)] = build_session([r for r in records if r['length_class']==size], args, cell_width)
    manifest = dict(pymol_versions=versions, source_script_sha256=checksum(Path(__file__)),
                    chain_map_sha256=checksum(args.chain_map),
                    cartoon_colors={'helix':'orange','strand':'blue','loop':'gray'},
                    model='first', orientation='right-handed PCA; longest axis vertical',
                    layout='4 columns; compact per-row heights; no molecular scaling',
                    notes='Independent structures, not aligned; labels give FASTA lengths; missing residues not rebuilt.',
                    targets=([r for r in previous.get('targets', []) if r['length_class'] not in args.groups]
                             + [r for r in records if r['length_class'] in args.groups]),
                    panel_layout='all groups independently zoomed; 4 columns; 22-pt bold PDB labels; separated title',
                    last_regenerated_groups=args.groups)
    (args.output_dir/'manifest.json').write_text(json.dumps(manifest,indent=2))


if __name__ == '__main__':
    main()
