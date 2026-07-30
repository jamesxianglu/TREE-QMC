using PhyloNetworks
using Printf


function fail(message)
    println(stderr, "ERROR: ", message)
    exit(2)
end


function nontrivial_split_count(tree)
    tree.numhybrids == 0 || fail("expected a tree, found $(tree.numhybrids) hybrid nodes")

    # After degree-2 nodes are suppressed, every non-pendant edge of an
    # unrooted phylogenetic tree represents one non-trivial split.
    return length(tree.edge) - tree.numtaxa
end


function prepare_true_tob(true_network_path, true_tob_output)
    true_network = readnewick(true_network_path)
    true_tob = treeofblobs(true_network)
    mkpath(dirname(abspath(true_tob_output)))
    writenewick(true_tob, true_tob_output)
    removedegree2nodes!(true_tob)
    return true_tob
end


function compare_tobs(true_tob, inferred_tree_path)
    inferred_tob = readnewick(inferred_tree_path)
    true_taxa = Set(tiplabels(true_tob))
    inferred_taxa = Set(tiplabels(inferred_tob))
    if true_taxa != inferred_taxa
        missing = sort!(collect(setdiff(true_taxa, inferred_taxa)))
        extra = sort!(collect(setdiff(inferred_taxa, true_taxa)))
        fail("taxon mismatch: missing=$(missing), extra=$(extra)")
    end

    # PhyloNetworks' RF routine counts clusters with multiplicity. Suppressing
    # degree-2 nodes first gives the usual split-set RF distance and handles
    # degree-2 nodes retained by treeofblobs.
    removedegree2nodes!(inferred_tob)

    rf = hardwiredclusterdistance(true_tob, inferred_tob, false)
    true_split_count = nontrivial_split_count(true_tob)
    inferred_split_count = nontrivial_split_count(inferred_tob)

    # RF = FP + FN, while inferred_count - true_count = FP - FN.
    twice_fp = rf + inferred_split_count - true_split_count
    twice_fn = rf + true_split_count - inferred_split_count
    if twice_fp < 0 || twice_fn < 0 || !iseven(twice_fp) || !iseven(twice_fn)
        fail("inconsistent RF and split counts")
    end

    fp = div(twice_fp, 2)
    fn = div(twice_fn, 2)
    taxon_count = length(true_taxa)
    maximum_rf = 2 * (taxon_count - 3)
    maximum_rf > 0 || fail("normalized RF requires at least four taxa")
    normalized_rf = rf / maximum_rf

    # Rates, each normalized by its own tree's internal branch count rather than
    # by the shared 2(n-3) denominator of the normalized RF:
    #   fnr = true internal branches missed      / true internal branches
    #   fpr = estimated internal branches spurious / estimated internal branches
    # A fully contracted estimate has no internal branches and hence no
    # opportunity to be wrong in the FP direction, so its fpr is reported as 0
    # (fp is necessarily 0 there too); the same convention covers a star true
    # tree of blobs. mean_rate is the unweighted average of the two, which is
    # the usual way to keep over- and under-contraction from being traded off
    # against each other by tree size.
    fnr = true_split_count > 0 ? fn / true_split_count : 0.0
    fpr = inferred_split_count > 0 ? fp / inferred_split_count : 0.0
    mean_rate = (fnr + fpr) / 2

    return (;
        taxon_count,
        fn,
        fp,
        normalized_rf,
        true_internal = true_split_count,
        estimated_internal = inferred_split_count,
        fnr,
        fpr,
        mean_rate,
    )
end


function print_comparison(result)
    println("\n-- constructed tree split comparison --")
    println("  false negatives      : ", result.fn)
    println("  false positives      : ", result.fp)
    println("  true internal edges  : ", result.true_internal)
    println("  est. internal edges  : ", result.estimated_internal)
    @printf("  normalized RF        : %.10f\n", result.normalized_rf)
    @printf("  false negative rate  : %.10f\n", result.fnr)
    @printf("  false positive rate  : %.10f\n", result.fpr)
    @printf("  mean of the two      : %.10f\n", result.mean_rate)
end


function run_batch(manifest_path, output_tsv; rates::Bool=false)
    true_tob_cache = Dict{String, Any}()
    mkpath(dirname(abspath(output_tsv)))

    open(output_tsv, "w") do output
        # Without --rates the columns are exactly what the row-sweep harness has
        # always parsed; the extra columns are opt-in so those callers keep
        # working unchanged.
        if rates
            println(output, "dataset\tdelta\tquery_alpha\tfn\tfp\tnormalized_rf" *
                "\ttrue_internal\testimated_internal\tfnr\tfpr\tmean_rate")
        else
            println(output, "dataset\tdelta\tquery_alpha\tfn\tfp\tnormalized_rf")
        end
        for (line_number, line) in enumerate(eachline(manifest_path))
            line_number == 1 && continue
            isempty(strip(line)) && continue
            fields = split(line, '\t'; keepempty=true)
            length(fields) == 6 || fail(
                "$(manifest_path):$(line_number): expected 6 tab-separated fields"
            )
            dataset, delta, query_alpha, true_network_path,
                inferred_tree_path, true_tob_output = fields

            true_tob = get!(true_tob_cache, true_network_path) do
                prepare_true_tob(true_network_path, true_tob_output)
            end
            result = compare_tobs(true_tob, inferred_tree_path)
            if rates
                @printf(
                    output,
                    "%s\t%s\t%s\t%d\t%d\t%.10f\t%d\t%d\t%.10f\t%.10f\t%.10f\n",
                    dataset,
                    delta,
                    query_alpha,
                    result.fn,
                    result.fp,
                    result.normalized_rf,
                    result.true_internal,
                    result.estimated_internal,
                    result.fnr,
                    result.fpr,
                    result.mean_rate,
                )
            else
                @printf(
                    output,
                    "%s\t%s\t%s\t%d\t%d\t%.10f\n",
                    dataset,
                    delta,
                    query_alpha,
                    result.fn,
                    result.fp,
                    result.normalized_rf,
                )
            end
        end
    end
end


function main()
    if length(ARGS) == 3 && ARGS[1] != "--batch"
        true_tob = prepare_true_tob(ARGS[1], ARGS[3])
        print_comparison(compare_tobs(true_tob, ARGS[2]))
        return
    end
    if length(ARGS) == 3 && ARGS[1] == "--batch"
        run_batch(ARGS[2], ARGS[3])
        return
    end
    if length(ARGS) == 4 && ARGS[1] == "--batch" && ARGS[2] == "--rates"
        run_batch(ARGS[3], ARGS[4]; rates=true)
        return
    end
    fail(
        "usage:\n" *
        "  julia compare_inferred_tob.jl <true-network> <inferred-tree> <true-tob-output>\n" *
        "  julia compare_inferred_tob.jl --batch [--rates] <manifest.tsv> <results.tsv>"
    )
end


main()
