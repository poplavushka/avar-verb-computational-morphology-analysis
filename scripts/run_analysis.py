# SPDX-License-Identifier: GPL-3.0-only
"""Extract fresh Qumin edits patterns and fit bottom-up MDL models."""
import argparse
from collections import defaultdict
import importlib.metadata
from itertools import combinations
import logging
import math
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
from types import SimpleNamespace

from common import (ROOT, groups_from_assignments, partition_hash, read_csv,
                    read_json, reconstruct_mdl, summarize, write_csv, write_json)


class OutputFiles:
    """Minimal file sink for Qumin's clustering export API."""
    def __init__(self, root):
        self.root = root

    def get_path(self, path):
        target = self.root / path
        target.parent.mkdir(parents=True, exist_ok=True)
        return str(target)

    def register_file(self, *args, **kwargs):
        pass


def label_groups(groups, prefix):
    ordered = sorted((sorted(group) for group in groups), key=lambda group: (-len(group), group))
    return {lid:f'{prefix}{i:02}' for i, group in enumerate(ordered, 1) for lid in group}


def run_model(config_path, target):
    if os.environ.get('PYTHONHASHSEED') != '0':
        raise RuntimeError('Start this computation through the main entry point')
    if importlib.metadata.version('qumin') != '4.0':
        raise RuntimeError('This analysis requires qumin==4.0')
    import numpy as np
    import pandas as pd
    from qumin.representations.patternstore import PatternStore, _find_cellpair_patterns
    from qumin.representations.segments import Inventory, Form
    from qumin.clustering import algorithms, descriptionlength, find_min_attribute

    logging.basicConfig(level=logging.WARNING)
    config = read_json(config_path)
    random.seed(1)
    np.random.seed(1)
    dataset = read_csv(ROOT / 'data/dataset.csv')
    ids = {row['lexeme_id'] for row in dataset}
    if config['sample'] == 'sample-593':
        ids = set((ROOT / 'data/sample-593.txt').read_text(encoding='utf-8').splitlines())
    if len(ids) != config['n_lexemes']:
        raise ValueError('Sample size does not match the configuration')
    inventory = Inventory.from_file(ROOT / 'data/segments' / (config['sample']+'.csv'))
    inputs = defaultdict(list)
    pairs = list(combinations(config['cells'], 2))
    for row in read_csv(ROOT / 'data/model-inputs.csv'):
        pair = row['cell_x'], row['cell_y']
        if row['lexeme_id'] in ids and pair in pairs:
            inputs[pair].append(row)
    store = PatternStore()
    store.cells = list(config['cells'])
    store.pairs = pairs
    store.algorithm = 'edits'
    pattern_rows = []
    for a, b in pairs:
        pair_rows = sorted(inputs[(a, b)], key=lambda row: row['lexeme_id'])
        if len(pair_rows) != len(ids) or {r['lexeme_id'] for r in pair_rows} != ids:
            raise ValueError(f'Incomplete input pair {a}/{b}')
        forms = []
        for row in pair_rows:
            lid = row['lexeme_id']
            for value in [row['form_x'], row['form_y']]:
                if not value.strip() or not set(value.split()) <= set(inventory._classes):
                    raise ValueError(f'Unknown or missing model tokens: {lid}')
            forms.append(dict(lexeme=lid, form_x=Form(row['form_x'], form_id=f'{lid}_{a}_1'), form_y=Form(row['form_y'], form_id=f'{lid}_{b}_1')))
        _, frame = _find_cellpair_patterns(inventory, pd.DataFrame(forms), (a, b), 'edits', False)
        if frame.pattern.isna().any():
            raise ValueError(f'Pattern inference failed for {a}/{b}')
        for row in frame.itertuples(index=False):
            pattern_rows.append(dict(lexeme_id=row.lexeme, cell_x=a, cell_y=b, pattern=re.sub(r' <[^>]+>$', '', repr(row.pattern))))
        # Qumin's macroclasses command also converts patterns to str before fitting.
        frame.pattern = frame.pattern.map(str)
        store[(a, b)] = frame

    result = target / 'results'
    write_csv(result / 'patterns.csv', pattern_rows)
    micro = store.to_microclasses(store.to_incidence_table())
    microlabel = label_groups(micro.values(), 'm')
    seed_rows = []
    for seed in config['seeds']:
        # Extraction and MDL use separately initialized random generators.
        random.seed(seed)
        np.random.seed(seed)
        folder = result if seed == 1 else result / 'seeds' / f'seed-{seed}'
        folder.mkdir(parents=True, exist_ok=True)
        node = algorithms.hierarchical_clustering(store, SimpleNamespace(frequencies=None), descriptionlength.BUDLClustersBuilder, md=OutputFiles(folder))
        macro = node.macroclasses()
        if not macro:
            raise ValueError('Qumin returned no macroclasses')
        macro_groups = [sorted(lid for rep in reps for lid in micro[rep]) for reps in macro.values()]
        macrolabel = label_groups(macro_groups, 'M')
        assignments = [dict(lexeme_id=lid, microclass=microlabel[lid], macroclass=macrolabel[lid]) for lid in sorted(ids)]
        write_csv(folder / 'assignments.csv', assignments)
        components = reconstruct_mdl(pattern_rows, assignments)
        minimum = float(find_min_attribute(node, 'DL'))
        if not math.isclose(components['mdl_bits'], minimum, rel_tol=0, abs_tol=1e-8):
            raise ValueError('Independent MDL reconstruction differs from Qumin')
        (folder / 'tree.txt').write_text(node.tree_string()+'\nBottom-up DL clustering on Min :'+str(minimum)+'\n', encoding='utf-8')
        write_csv(folder / 'mdl.csv', [dict(seed=seed, n_lexemes=len(ids), n_cells=len(config['cells']), n_pairs=len(pairs), microclasses=len(micro), macroclasses=len(macro_groups), **components)])
        seed_rows.append(dict(seed=seed, microclasses=len(micro), macroclasses=len(macro_groups), mdl_bits=minimum, micro_partition_sha256=partition_hash(micro.values()), macro_partition_sha256=partition_hash(macro_groups)))
    write_csv(result / 'seed-checks.csv', seed_rows)
    print(f"Completed {config['experiment']}/{config['sample']}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'reproduced', help='New output directory; it must not exist')
    parser.add_argument('--experiment', choices=['four-forms', 'five-forms-masdar', 'five-forms-converb', 'six-forms'], help='Run one form set (both samples)')
    parser.add_argument('--worker', nargs=2, metavar=('CONFIG', 'OUTPUT'), help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        run_model(Path(args.worker[0]), Path(args.worker[1]))
        return
    target = args.output.resolve()
    if target.exists():
        parser.error('The output directory already exists; choose a fresh directory')
    configs = sorted((ROOT / 'experiments').glob('*/*/config.json'))
    if args.experiment:
        configs = [p for p in configs if p.parent.parent.name == args.experiment]
    # Validate supplied data and reference outputs before starting new fits.
    from verify_results import verify_dataset
    verify_dataset()
    target.mkdir(parents=True)
    env = {**os.environ, 'PYTHONHASHSEED':'0', 'PYTHONDONTWRITEBYTECODE':'1', 'MPLBACKEND':'Agg'}
    for config in configs:
        relative = config.parent.relative_to(ROOT)
        destination = target / relative
        destination.mkdir(parents=True)
        shutil.copy2(config, destination / 'config.json')
        print('Running '+str(relative.relative_to('experiments')), flush=True)
        with (destination / 'calculation.log').open('w', encoding='utf-8') as log:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker', str(config), str(destination)], env=env, cwd=target, stdout=log, stderr=subprocess.STDOUT, check=True)
    summarize(target)
    from verify_results import verify_results
    report = verify_results(target, allow_partial=bool(args.experiment))
    report['patterns_freshly_extracted'] = True
    report['dependency_versions'] = {name:importlib.metadata.version(name) for name in ['qumin', 'numpy', 'pandas', 'paralex']}
    report['python'] = '.'.join(map(str, sys.version_info[:3]))
    write_json(target / 'verification.json', report)
    print(f"Verified {report['models_checked']} models and {report['mdl_runs_checked']} MDL runs.", flush=True)


if __name__ == '__main__':
    main()
