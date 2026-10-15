# Avar Verb Computational Morphology Analysis

Dataset, code and calculation outputs for pairwise pattern extraction, microclass identification and bottom-up minimum-description-length (MDL) clustering of Avar verbs using Qumín 4.0.

## Experiments

Each experiment contains `sample-700/` and `sample-593/`, with a configuration and results.

| Directory | Forms |
| --- | --- |
| `experiments/four-forms/` | INF, PRS, FUT, AOR |
| `experiments/five-forms-masdar/` | INF, PRS, FUT, AOR, VN |
| `experiments/five-forms-converb/` | INF, PRS, FUT, AOR, CVB.PFV |
| `experiments/six-forms/` | INF, PRS, FUT, AOR, CVB.PFV, VN |

All eight models use freshly extracted `edits` patterns. The six-form models additionally include MDL runs with seeds 2–5. Calculations include class assignments, comparisons between experiments and samples, and an independent reconstruction of MDL components.

## Files

- `data/dataset.csv`: lexical identifiers and six recorded forms.
- `data/model-inputs.csv`: prepared, segmented inputs for each cell pair.
- `data/sample-593.txt`: identifiers for the smaller sample.
- `data/segments/`: identity-based segment inventories for the two samples.
- `experiments/*/*/results/`: patterns, complete assignments, MDL components, trees and seed checks.
- `results/`: aggregate model tables and partition comparisons.
- `scripts/`: analysis and verification code.
- [METHODS.md](METHODS.md): computational procedure.

## Reproduction

Python 3.12:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock.txt
python scripts/run_analysis.py --output reproduced
python scripts/verify_results.py --results reproduced
```

The output directory must be new. The command recomputes patterns and MDL from the supplied model inputs and checks them against the reference calculations.

To verify the supplied results without installing additional packages:

```sh
python3 scripts/verify_results.py
```

To run one form set for both samples:

```sh
python scripts/run_analysis.py --experiment six-forms --output reproduced-six-forms
```

## License and dependencies

Analysis code: GNU GPL v3; see [LICENSE](LICENSE). Third-party authorship and licensing are listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). This anonymous submission package does not identify the authors of the study.
