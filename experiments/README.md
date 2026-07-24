# Row-sweep split-test experiments

This folder evaluates `SpeciesTree::row_sweep_test_idx`, which predicts whether
a candidate bipartition `A|B` is a split in the **tree of blobs** (ToB). The
experiments use the CAMUS simulated datasets under `data/camus-dataset`.

## What the implementation is testing

For each four-taxon set `{x,y,rho,r}`, the theoretical oracle returns an exact
labelled quarnet. The algorithm counts a contradiction whenever
`Query(x,y,rho,r) != xy|rho r`. Therefore a 4-blob answer is counted, but a
different quartet-tree topology is also counted.

The C++ implementation has gene-tree counts rather than an exact oracle. It
uses MSCquartets' T1 test as a statistical surrogate:

* `true`: the specific quartet tree `xy|rho r` was not rejected;
* `false`: T1 rejected it (or no resolved quartet was observed), so row sweep
  counts one contradiction.

A rejection is not a unique 4-blob diagnosis, and finite data can fail to
reject a real 4-blob. The experiment measures how well this surrogate supports
the full row-sweep decision.

The experiment flow is:

1. Compute trusted ToB splits from the true network with PhyloNetworks.
2. Generate true splits plus hard and random non-splits.
3. Run `tree-qmc --rowsweep` on every candidate.
4. Compare ACCEPT/REJECT predictions with the trusted labels.

## Ground truth

The true non-trivial splits are computed from `true_net.nwk` by
PhyloNetworks' `treeofblobs` through `compute_tree_of_blob.jl`. The resulting
tree's internal edges define the reference splits. `tree_of_blobs.py` provides
a dependency-free bridge-based fallback and structural diagnostics.

* **Positives** — the true non-trivial ToB splits (should be ACCEPTed).
* **Negatives** — a *mixed* set of non-splits (should be REJECTed): half
  **hard** (a true split perturbed by moving one taxon) and half **random**.

`--rowsweep_out` writes `1=ACCEPT` and `0=REJECT`:

| ground truth | correct | error | name |
|--------------|---------|-------|------|
| positive (true split) | 1 | 0 | **false negative** |
| negative (non-split)  | 0 | 1 | **false positive** |

## Parameters

| parameter | meaning |
|---|---|
| `delta` | Oracle-noise bound used by row sweep to set `theta=(1+delta)/4`; the proof assumes `delta < 1/3`. |
| `--query-alpha` | Per-quartet T1 rejection level used by the empirical query surrogate. |

The default `query-alpha` is `0.0005`.
| `--eps` | Desired error bound used only in the report's theoretical sample-size calculation. |

`delta` and `query-alpha` are not interchangeable. Lower `query-alpha` yields
fewer T1 rejections: this protects true splits but reduces power against
4-blobs and other contradictions. Calibrate it on held-out simulations. For a
conservative family-wise level of 0.05 over `m` planned quartet tests, a simple
starting point is `0.05/m`.

## One dataset

```bash
./run_split_test.sh n25/10 0.15                 # one delta
./run_split_test.sh n25/10 0.05,0.15,0.3        # sweep deltas
./run_split_test.sh n15/00 0.15 g_true.nwk      # choose gene trees
./run_split_test.sh n25/10 0.15 iqtree_500.nwk --neg-mult 4 --seed 7
./run_split_test.sh n25/10 0.15 iqtree_500.nwk --query-alpha 0.001
```

Arguments are
`<dataset-path> <delta[,delta2,...]> [gene_trees=iqtree_500.nwk] [extra flags]`.
The report includes FP/FN rates,
each misclassified split, blob adjacency, and the theorem's requirement on the
larger-side size `s`.

To compile, run the classifier metrics, construct the row-sweep ToB, and
compare it with the trusted ToB in one command:

```bash
./compile_and_test_rowsweep.sh n15/00 0.15 0.0005
```

For a paired comparison, pass the row-sweep refinement to 3f1a so both methods
start from exactly the same candidate edges:

```bash
REF=results/compile_test/n15_00_d0.15_a0.0005/refinement.nwk
./compile_and_test_rowsweep.sh n15/00 0.15 0.0005 iqtree_500.nwk "$REF"
./compile_and_test_3f1a.sh n15/00 1e-7 0.95 iqtree_500.nwk "$REF"
```

The tree comparison reports false negatives, false positives, and normalized
unrooted Robinson–Foulds distance. The normalization denominator is the maximum
RF distance for unrooted trees on `n` taxa, `2(n-3)`.

To tune row-sweep on a seed-controlled random sample of 5 n15 and 3 n25
datasets, run:

```bash
./sweep_rowsweep_hyperparameters.sh
```

The default grid has 7 delta values from 0.01 through 0.30 and 10 log-spaced
query-alpha values from 0.000001 through 0.05. The script builds once and uses
one refinement per dataset. Its only persistent output is `comparisons.tsv`,
with dataset, delta, query-alpha, FN, FP, and normalized RF columns. All trees
and logs are temporary. Set `SEED`, `DELTAS`, `QUERY_ALPHAS`, or `GENE_TREES`
to override the defaults.

## Many datasets × many deltas

```bash
export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"; export R_HOME="$(R RHOME)"
python3 run_experiment.py --group n15 --datasets 10          # pooled FNR/FPR per delta
python3 run_experiment.py --group n25 --datasets 10 --deltas 0.05,0.15,0.3
```

Outputs under `results/…` include:

* `detail_a<alpha>.csv`: one row per dataset and delta;
* `summary_by_delta_a<alpha>.csv`: pooled rates by delta;
* `bips/`, `preds/`, `logs/`, and `tob/`: intermediate and diagnostic files.

Alpha is part of output names so runs at different T1 thresholds do not
silently overwrite one another.

## Files

| file | role |
|------|------|
| `compute_tree_of_blob.jl` | PhyloNetworks `treeofblobs` → ToB tree (trusted ground truth). |
| `tree_of_blobs.py` | Newick parsing, bridge-based ToB splits (fallback), blob-adjacency diagnostics. |
| `gen_bipartitions.py` | positives + mixed negatives → labelled TSV. |
| `common.py` | path/environment resolution, Julia ToB call, TSV parsing, execution, and named score metrics. |
| `analyze_split.py` | single dataset+delta analysis with weakness report. |
| `run_experiment.py` | multi-dataset × multi-delta sweep. |
| `run_split_test.sh` | shell entry point (the one to use day-to-day). |
| `compile_and_test_rowsweep.sh` | compile, classify candidate splits, construct a ToB, and report metrics. |
| `compile_and_test_3f1a.sh` | compile, construct a 3f1a ToB, and report comparable split metrics. |
| `compare_inferred_tob.jl` | use PhyloNetworks to report FN, FP, and normalized unrooted Robinson–Foulds distance. |
| `sweep_rowsweep_hyperparameters.sh` | tune delta and query-alpha on sampled n15/n25 datasets. |

## Prerequisites

* `tree-qmc` built at `../build/tree-qmc` (see build notes below).
* R with `Rcpp`, `RInside`, `MSCquartets` (linked by the binary).
* Julia with `PhyloNetworks` (`import Pkg; Pkg.add("PhyloNetworks")`).

### Building tree-qmc on macOS / Apple clang

```bash
export PATH="/opt/homebrew/bin:$PATH"
cmake -S TREE-QMC -B TREE-QMC/build -DCMAKE_BUILD_TYPE=Release && cmake --build TREE-QMC/build -j4
```

If R package compilation fails with `invalid value 'gnu23'`, add `~/.R/Makevars`
with `CC=clang -std=gnu17` and reinstall the packages.

## Scope of the C++ changes

The topology-specific `pvalue_t1` helper and row-sweep-only qCF/T1 caches are
called through `query_pairs_together -> row_sweep_test_idx`, from either the
split experiment or the row-sweep `SpeciesTree` constructor.
Existing T3 p-values used by 3f1a, 2f2a, and the other tree-of-blobs searches
remain unchanged.
