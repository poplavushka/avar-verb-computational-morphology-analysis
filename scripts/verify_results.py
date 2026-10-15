# SPDX-License-Identifier: GPL-3.0-only
"""Verify data, complete partitions, fresh pattern equivalence and MDL scores."""
import argparse
from collections import defaultdict
from itertools import combinations
import math
from pathlib import Path
import re

from common import (ROOT, groups_from_assignments, partition_hash, pattern_hash,
                    read_csv, read_json, reconstruct_mdl, write_json)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_dataset():
    dataset = read_csv(ROOT / 'data/dataset.csv')
    cells = ['inf', 'prs', 'fut', 'aor', 'cvb.pfv', 'vn']
    ids = {row['lexeme_id'] for row in dataset}
    require(len(ids) == len(dataset) == 700, 'Invalid dataset identifiers')
    require(list(dataset[0]) == ['lexeme_id']+cells, 'Unexpected dataset columns')
    subset_list = (ROOT / 'data/sample-593.txt').read_text(encoding='utf-8').splitlines()
    subset = set(subset_list)
    require(len(subset_list) == len(subset) == 593 and subset <= ids, 'Invalid reduced sample')
    inputs = read_csv(ROOT / 'data/model-inputs.csv')
    require(list(inputs[0]) == ['lexeme_id', 'cell_x', 'cell_y', 'form_x', 'form_y'], 'Unexpected model input columns')
    expected = {(lid, a, b) for lid in ids for a,b in combinations(cells,2)}
    keys = [(r['lexeme_id'],r['cell_x'],r['cell_y']) for r in inputs]
    require(len(keys) == len(expected) and set(keys) == expected, 'Incomplete or duplicate model inputs')
    inventory = {r['sound_id'] for r in read_csv(ROOT / 'data/segments/sample-700.csv')}
    for row in inputs:
        for field in ['form_x', 'form_y']:
            require(bool(row[field].strip()) and set(row[field].split()) <= inventory, 'Unknown or empty model token sequence')
    smaller_inventory = {r['sound_id'] for r in read_csv(ROOT / 'data/segments/sample-593.csv')}
    require(smaller_inventory <= inventory, 'Inconsistent segment inventories')
    require(all(set(r[field].split()) <= smaller_inventory for r in inputs if r['lexeme_id'] in subset for field in ['form_x','form_y']), 'Unknown token in smaller sample')
    return ids, subset


def verify_results(output, allow_partial=False):
    full, reduced = verify_dataset()
    output = Path(output)
    configs = sorted((output / 'experiments').glob('*/*/config.json'))
    require(bool(configs), 'No experiment configurations found')
    require(allow_partial or len(configs) == 8, 'Expected eight experiments')
    errors = []
    run_count = 0
    pattern_records = 0
    for path in configs:
        config = read_json(path)
        name = config['experiment']+'/'+config['sample']
        baseline = ROOT / 'experiments' / config['experiment'] / config['sample']
        require(config == read_json(baseline / 'config.json'), f'Configuration changed: {name}')
        expected = read_json(baseline / 'expected.json')
        ids = full if config['sample'] == 'sample-700' else reduced
        result = path.parent / 'results'
        patterns = read_csv(result / 'patterns.csv')
        pairs = set(combinations(config['cells'], 2))
        keys = [(r['lexeme_id'],r['cell_x'],r['cell_y']) for r in patterns]
        require(len(keys) == len(ids)*len(pairs) and set(keys) == {(lid,a,b) for lid in ids for a,b in pairs}, f'Incomplete patterns: {name}')
        require(pattern_hash(patterns) == expected['patterns_sha256'], f'Patterns differ from the reference calculation: {name}')
        pattern_records += len(patterns)
        seed_checks = read_csv(result / 'seed-checks.csv')
        require([int(r['seed']) for r in seed_checks] == config['seeds'], f'Incomplete seed checks: {name}')
        for reference in expected['seeds']:
            seed = reference['seed']
            folder = result if seed == 1 else result / 'seeds' / f'seed-{seed}'
            assignments = read_csv(folder / 'assignments.csv')
            require(len(assignments) == len(ids) and {r['lexeme_id'] for r in assignments} == ids, f'Incomplete assignments: {name}, seed={seed}')
            require(list(assignments[0]) == ['lexeme_id','microclass','macroclass'], f'Unexpected assignment columns: {name}')
            components = reconstruct_mdl(patterns, assignments)
            saved_rows = read_csv(folder / 'mdl.csv')
            require(len(saved_rows) == 1, f'Expected one MDL record: {name}, seed={seed}')
            saved = saved_rows[0]
            require(int(saved['seed']) == seed and int(saved['n_lexemes']) == len(ids), f'Incorrect MDL metadata: {name}, seed={seed}')
            require(int(saved['n_cells']) == len(config['cells']) and int(saved['n_pairs']) == len(pairs), f'Incorrect model dimensions: {name}, seed={seed}')
            for level in ['microclass','macroclass']:
                groups = groups_from_assignments(assignments, level)
                require(len(groups) == expected[level+'es'], f'Incorrect class count: {name}, seed={seed}')
                require(int(saved[level+'es']) == len(groups), f'Incorrect saved class count: {name}, seed={seed}')
                require(partition_hash(groups) == reference[level.removesuffix('class')+'_partition_sha256'], f'Changed {level} partition: {name}, seed={seed}')
            minimum = re.search(r'Min\s*:\s*([-+0-9.eE]+)',(folder / 'tree.txt').read_text(encoding='utf-8'))
            require(minimum is not None, f'Missing Qumin tree score: {name}, seed={seed}')
            scores = [float(saved['mdl_bits']), float(reference['mdl_bits']), float(minimum.group(1))]
            for reported in scores:
                error = abs(components['mdl_bits']-reported)
                errors.append(error)
                require(error < 1e-8, f'MDL reconstruction failed: {name}, seed={seed}')
            for component in ['M','C','P','R']:
                require(math.isclose(components[component],float(saved[component]),rel_tol=0,abs_tol=1e-8), f'Incorrect component {component}: {name}, seed={seed}')
            check = next(r for r in seed_checks if int(r['seed']) == seed)
            require(check['micro_partition_sha256'] == reference['micro_partition_sha256'] and check['macro_partition_sha256'] == reference['macro_partition_sha256'], f'Incorrect seed table: {name}, seed={seed}')
            require(math.isclose(float(check['mdl_bits']),float(saved['mdl_bits']),rel_tol=0,abs_tol=1e-8), f'Incorrect seed score: {name}, seed={seed}')
            run_count += 1
    # Check that aggregate tables are fully derived from the experiment CSVs.
    from common import summarize
    import tempfile
    import shutil
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        shutil.copytree(output / 'experiments', temp / 'experiments')
        summarize(temp)
        for name in ['model-summary.csv','comparisons.csv']:
            require(read_csv(temp / 'results' / name) == read_csv(output / 'results' / name), f'Aggregate table differs from complete memberships: {name}')
    return dict(passed=True, models_checked=len(configs), mdl_runs_checked=run_count, pattern_records_checked=pattern_records, pattern_assignments_match_reference=True, complete_partitions_match_reference=True, maximum_mdl_absolute_error=max(errors), tolerance=1e-8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT, help='Repository or reproduced output root')
    parser.add_argument('--report', type=Path, help='Write a JSON verification report')
    args = parser.parse_args()
    report = verify_results(args.results.resolve())
    if args.report:
        write_json(args.report, report)
    import json
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
