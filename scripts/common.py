# SPDX-License-Identifier: GPL-3.0-only
"""File formats, partition comparisons and independent MDL calculations."""
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


def partition_hash(groups):
    data = sorted(sorted(group) for group in groups)
    return hashlib.sha256(json.dumps(data, separators=(',', ':')).encode()).hexdigest()


def pattern_hash(rows):
    data = sorted((r['lexeme_id'], r['cell_x'], r['cell_y'], r['pattern']) for r in rows)
    return hashlib.sha256(json.dumps(data, ensure_ascii=False, separators=(',', ':')).encode()).hexdigest()


def groups_from_assignments(rows, level):
    groups = defaultdict(list)
    for row in rows:
        groups[row[level]].append(row['lexeme_id'])
    return [sorted(group) for group in groups.values()]


def entropy_cost(count, total):
    return -count * math.log2(count / total) if count else 0.0


def reconstruct_mdl(patterns, assignments):
    """Reconstruct M+C+P+R using only CSVs and the Python standard library."""
    ids = [row['lexeme_id'] for row in assignments]
    if len(ids) != len(set(ids)) or not ids:
        raise ValueError('Assignments must contain distinct lexical identifiers')
    micro = defaultdict(list)
    macro_of = {}
    for row in assignments:
        micro[row['microclass']].append(row['lexeme_id'])
        macro_of[row['lexeme_id']] = row['macroclass']
    macro = defaultdict(list)
    for members in micro.values():
        labels = {macro_of[lid] for lid in members}
        if len(labels) != 1:
            raise ValueError('A microclass crosses macroclass boundaries')
        macro[next(iter(labels))].append(sorted(members)[0])
    pairs = defaultdict(dict)
    for row in patterns:
        key = row['cell_x'], row['cell_y']
        if row['lexeme_id'] in pairs[key]:
            raise ValueError('Duplicate lexical identifier in a cell pair')
        pairs[key][row['lexeme_id']] = row['pattern']
    if not pairs or any(set(table) != set(ids) for table in pairs.values()):
        raise ValueError('Pattern memberships are incomplete')
    # Check both directions: signatures define exactly the reported microclasses.
    signatures = defaultdict(list)
    for lid in ids:
        signatures[tuple(pairs[pair][lid] for pair in sorted(pairs))].append(lid)
    if partition_hash(signatures.values()) != partition_hash(micro.values()):
        raise ValueError('Microclasses differ from complete pattern signatures')
    M = sum(entropy_cost(len(members), len(ids)) for members in micro.values())
    C = sum(entropy_cost(len(reps), len(micro)) for reps in macro.values())
    P = R = 0.0
    for table in pairs.values():
        counters = [Counter(table[lid] for lid in reps) for reps in macro.values()]
        R += sum(entropy_cost(n, sum(counter.values())) for counter in counters for n in counter.values())
        distribution = Counter(pattern for counter in counters for pattern in counter)
        total = sum(distribution.values())
        P += sum(entropy_cost(n, total) for n in distribution.values())
        P += sum(entropy_cost(len(counter), total) for counter in counters)
    return dict(M=M, C=C, P=P, R=R, mdl_bits=M+C+P+R)


def adjusted_rand(left, right):
    ids = sorted(set(left) & set(right))
    choose2 = lambda n: n * (n - 1) / 2
    if len(ids) < 2:
        return 1.0
    intersections = Counter((left[lid], right[lid]) for lid in ids)
    rows = Counter(left[lid] for lid in ids)
    cols = Counter(right[lid] for lid in ids)
    observed = sum(choose2(n) for n in intersections.values())
    a = sum(choose2(n) for n in rows.values())
    b = sum(choose2(n) for n in cols.values())
    expected = a * b / choose2(len(ids))
    maximum = (a + b) / 2
    return (observed - expected) / (maximum - expected) if maximum != expected else 1.0


def summarize(output):
    output = Path(output)
    summaries = []
    assignments = {}
    for path in sorted((output / 'experiments').glob('*/*/config.json')):
        config = read_json(path)
        key = config['experiment'], config['sample']
        result = path.parent / 'results'
        summaries.append(dict(experiment=key[0], sample=key[1], **read_csv(result / 'mdl.csv')[0]))
        assignments[key] = read_csv(result / 'assignments.csv')
    write_csv(output / 'results/model-summary.csv', summaries)
    comparisons = []
    for key_a, key_b in combinations(sorted(assignments), 2):
        if key_a[0] != key_b[0] and key_a[1] != key_b[1]:
            continue
        left = {row['lexeme_id']: row for row in assignments[key_a]}
        right = {row['lexeme_id']: row for row in assignments[key_b]}
        common = sorted(set(left) & set(right))
        row = dict(experiment_a=key_a[0], sample_a=key_a[1], experiment_b=key_b[0], sample_b=key_b[1], n_common=len(common))
        for level in ['microclass', 'macroclass']:
            a = {lid:left[lid][level] for lid in common}
            b = {lid:right[lid][level] for lid in common}
            row[level+'_ari'] = adjusted_rand(a, b)
            row[level+'_same_partition'] = partition_hash(groups_from_assignments([left[lid] for lid in common], level)) == partition_hash(groups_from_assignments([right[lid] for lid in common], level))
        comparisons.append(row)
    write_csv(output / 'results/comparisons.csv', comparisons)
