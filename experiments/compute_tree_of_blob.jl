using PhyloNetworks


function main()
    if length(ARGS) != 2
        println("Usage: julia compute_tree_of_blob.jl <input_network_file> <output_tree_of_blobs_file>")
        return
    end

    input_network_file = ARGS[1]
    output_tree_of_blobs_file = ARGS[2]

    # Read the network from the input file
    net = readnewick(input_network_file)
    # println("Input Network:")
    # println(net)

    # Compute the tree of blobs
    tob = treeofblobs(net)
    # println("Computed Tree of Blobs:")
    # println(tob)

    # Write the tree of blobs to the output file
    writenewick(tob, output_tree_of_blobs_file)
    println("Tree of blobs written to: $output_tree_of_blobs_file")
end

main()
