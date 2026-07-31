# experiments

Tree-of-blobs reconstruction experiments: run a method on a set of CAMUS
replicates, compare each estimate against the true tree of blobs, and collect
the results into one CSV per configuration.

Three scripts do all of it, and none of them is method-specific:

| script | role |
|---|---|
| `submit_tob.sh` | expand a dataset selection and launch one job per replicate, on SLURM or locally |
| `run_tob.sbatch` | one replicate, one configuration -> estimate + one-row CSV shard |
| `collect_tob.sh` | merge a run's shards into `results/<run label>.csv` |

Everything that distinguishes one experiment from another is passed in, so
adding a method means adding one line to the `case` in `run_tob.sbatch`, not
writing another pair of scripts.

## Running one

Nothing is defaulted: a result is only reproducible if its configuration was
written down, so it has to be typed out.

```sh
./submit_tob.sh --method rowsweep --runner slurm --gene-trees g_500.nwk \
    --datasets 'n15 n25 n50 n100 n150 n200' \
    --args '--oracle t1+cf|maj --rowsweep-tau 0.45,0.60 --rowsweep-heavy 1,2 --query-alpha 0.001'
```

`--runner local` runs the same jobs serially instead of submitting them, which
is how the small groups are usually done. `--datasets` takes groups, optionally
narrowed: `'n15 n50:00-09 n100:07'`. Add `--dry-run` to see the plan first.

The run lands in `results/runs/<label>/`, where the label is a slug of the
method and its flags, and the binary is copied in so that a later rebuild cannot
change what the numbers mean. Then:

```sh
./collect_tob.sh --run-dir results/runs/rowsweep__oracle-t1+cf-maj-rowsweep-tau-0.45-0.60-...
```

writes `results/<same label>.csv`. The name carries the method, the oracle and
the parameters, so a listing of `results/` says which configurations have been
measured.

Reproducing the originally published numbers is just another configuration:

```sh
./submit_tob.sh --method rowsweep --runner local --gene-trees iqtree_500.nwk \
    --datasets 'n15' --args '--delta 0.25 --query-alpha 0.001'
./submit_tob.sh --method cornerrow --runner local --gene-trees iqtree_500.nwk \
    --datasets 'n15' \
    --args '--corner-k 0 --heavy-sampling 2 --corner-tau 0.3125 --corner-seed 20250729 --query-alpha 0.001'
```

Both were checked against `results/rowsweep.csv` and `results/cornerrow.csv`
and agree exactly.

## What the methods are testing

For each four-taxon set `{x,y,rho,r}` the theoretical oracle returns an exact
labelled quarnet, and a contradiction is `Query(x,y,rho,r) != xy|rho r`. So a
4-blob answer counts, and so does a different quartet-tree topology.

The implementation has gene-tree counts rather than an oracle, so it substitutes
a test, chosen with `--oracle`; see `results/analysis/README.md` for the
available approximations and `doc/oracle_and_row_rule.tex` for why the tuned
configurations look the way they do. The tuned settings measured so far:

```sh
# row sweep
--oracle 't1+cf|maj' --rowsweep-tau 0.45,0.60 --rowsweep-heavy 1,2
# corner row (the cf guard is inert here, and t3 is worse than t1)
--oracle 't1|maj' --corner-tau 0.3125,0.60
```

## The rest

| path | what it is |
|---|---|
| `compute_tree_of_blob.jl` | ground truth: PhyloNetworks `treeofblobs` on `true_net.nwk`. Called by `run_tob.sbatch` and cached in `results/true_tob/`, since it depends only on the dataset |
| `dump_rowsweep_evidence.sh` | records the qCF counts and T1 p-value of every 4-set the sweep queries, so decision rules can be scored offline without rerunning anything |
| `split_test/` | a *different* experiment: does `row_sweep_test_idx` classify a given bipartition correctly, scored against labelled positives and negatives. Driven by `--rowsweep FILE`; entry point `split_test/run_split_test.sh` |

`results/` here holds scratch output from `split_test/`; the benchmark CSVs live
in the top-level `results/`.
