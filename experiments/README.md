# row_sweep_test_idx assessment harness

Evaluates `SpeciesTree::row_sweep_test_idx` — which decides whether a taxon
bipartition `A|B` is a split of the **tree of blobs** of a network — against
ground truth from the camus simulated datasets (`data/camus-dataset`).

## Ground truth

The set of true non-trivial splits is computed from `true_net.nwk` by
**PhyloNetworks' `treeofblobs`** (`compute_tree_of_blob.jl`) — the trusted
reference. We parse the resulting tree's internal edges into splits. A pure
Python bridge-finder (`tree_of_blobs.py`) is kept as a fast fallback and was
cross-validated to give identical splits.

* **Positives** — the true non-trivial ToB splits (should be ACCEPTed).
* **Negatives** — a *mixed* set of non-splits (should be REJECTed): half
  **hard** (a true split perturbed by moving 1–2 taxa) and half **random**.

Scoring (`--rowsweep_out` gives `1=ACCEPT` / `0=REJECT`):

| ground truth | correct | error | name |
|--------------|---------|-------|------|
| positive (true split) | 1 | 0 | **false negative** |
| negative (non-split)  | 0 | 1 | **false positive** |

## Quick start — one dataset, with a weakness report

```bash
./run_split_test.sh n25/10 0.15                 # one delta
./run_split_test.sh n25/10 0.05,0.15,0.3        # sweep deltas
./run_split_test.sh n15/00 0.15 g_true.nwk      # choose gene trees
./run_split_test.sh n25/10 0.15 iqtree_500.nwk --neg-mult 4 --seed 7
```

Args: `<dataset-path> <delta[,delta2,...]> [gene_trees=iqtree_500.nwk] [extra flags]`.
The script sets `PATH`/`R_HOME` and calls `analyze_split.py`, which prints
FP/FN rates plus a structural breakdown of every error (blob adjacency, split
size, hard-vs-random negatives, nearest true split) and a `theory reference`
line comparing the larger-side size `s` to the Theorem's requirement.

## Many datasets × many deltas

```bash
export PATH="$HOME/.juliaup/bin:/opt/homebrew/bin:$PATH"; export R_HOME="$(R RHOME)"
python3 run_experiment.py --group n15 --datasets 10          # pooled FNR/FPR per delta
python3 run_experiment.py --group n25 --datasets 10 --deltas 0.05,0.15,0.3
```

Outputs (`results/…`): `detail.csv` (per dataset×delta), `summary_by_delta.csv`
(pooled), plus `bips/`, `preds/`, `logs/`, `tob/`.

## Files

| file | role |
|------|------|
| `compute_tree_of_blob.jl` | PhyloNetworks `treeofblobs` → ToB tree (trusted ground truth). |
| `tree_of_blobs.py` | Newick parsing, bridge-based ToB splits (fallback), blob-adjacency diagnostics. |
| `gen_bipartitions.py` | positives + mixed negatives → labelled TSV. |
| `common.py` | path/env resolution, Julia ToB call, run/score helpers. |
| `analyze_split.py` | single dataset+delta analysis with weakness report. |
| `run_experiment.py` | multi-dataset × multi-delta sweep. |
| `run_split_test.sh` | shell entry point (the one to use day-to-day). |

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

## Changes made to the TREE-QMC source (none alter row_sweep_test_idx logic)

*Build-only portability (Apple clang 16 / libc++):*
1. `CMakeLists.txt`: `CMAKE_CXX_STANDARD 20 → 17` (no C++20 features are used;
   C++20 + new libc++ exposed an MQLib circular-include bug).
2. `src/csvparser.hpp`: VLA-with-initializer → `std::vector<int>`.
3. `src/tree_of_blobs.cpp` (blob-search path): `std::random_shuffle` (removed in
   C++17 libc++) → `std::shuffle` + `mt19937`.

*Required init for the `--rowsweep` path (previously missing — it segfaulted):*
4. top of `run_split_experiment`: `add_r_libpaths_and_load(RINS)` (loads
   MSCquartets, needed by `pvalue()`) and `for (Tree *t : input)
   t->LCA_preprocessing();` (builds the LCA index `Tree::get_quartet` needs).
   Every other quartet path already did both.
