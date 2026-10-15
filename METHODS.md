# Computational methods

The supplied dataset contains 700 lexical identifiers and six forms per identifier. Four cell sets are evaluated on each of the two supplied samples (700 and 593 identifiers), giving eight models. The smaller sample is specified by its identifier list. Lexical identifiers are retained separately throughout the calculation.

The computational starting point is `data/model-inputs.csv`: one row per lexical identifier and unordered pair of cells, with a segmented input form for each side. These prepared pairwise inputs are used as supplied; preprocessing is outside the scope of this package. The recorded forms in `data/dataset.csv` identify the lexical material. Model inputs contain no class labels or annotation features.

For each experiment, the code selects the requested cells and identifiers, then extracts new pairwise alternation patterns with Qumín 4.0's `edits` algorithm and the supplied identity-based segment inventory. Inference is sequential. Patterns are not shared between fitted experiments and no saved pattern assignments are used for extraction.

Microclasses are equivalence classes of lexical identifiers with identical complete vectors of pairwise patterns. Qumín's `BUDLClustersBuilder` performs bottom-up hierarchical clustering to obtain macroclasses. The default Qumín description-length objective is retained: lexical membership sizes contribute to M, and microclasses are weighted uniformly in C, P and R. Reported class sizes are expanded to complete lexical memberships.

Python and NumPy random generators are initialized explicitly, with `PYTHONHASHSEED=0`. Each primary MDL calculation uses seed 1. For the six-form experiments, seeds 2–5 repeat MDL clustering on the freshly extracted patterns from the same experiment. Full assignments and a tree are saved for every seed.

Verification compares every lexical pattern assignment and both complete partitions with the reference calculations, independently of arbitrary class labels. A separate standard-library calculation reconstructs M, C, P, R and their sum from the exported patterns and memberships. Its result is checked against the Qumín tree minimum and the saved score with an absolute tolerance of 1e-8 bits.

Partition comparisons use the adjusted Rand index on shared lexical identifiers. Aggregate tables are recomputed from complete experiment outputs during verification.
