#include "instance.hpp"


extern bool DEBUG_MODE;

Instance::Instance(int argc, char **argv) {
    DEBUG_MODE = false;

    input_file = "";
    output_file = "";
    mapping_file = "";
    stree_file = "";
    table_file = "";
    output_qcfs_table_file = "";
    rowsweep_file = "";
    rowsweep_out_file = "";
    rowsweep_dump_prefix = "";
    rowsweep_dump_all_targets = false;
    rowsweep_dump_all_anchors = false;
    rowsweep_delta = 0.25;
    rowsweep_query_alpha = 0.001;
    corner_row_k = 0;        // 0 selects the largest legal row sample size, n-3
    corner_row_heavy = 1;    // each row draws min(heavy*k, row population)
    corner_row_tau_spec = "";  // empty selects the default (1+delta)/4
    corner_row_seed = 20250729;
    corner_row_cross = false;
    corner_row_corroborate_spec = "";
    rowsweep_row_mode_spec = "fixed";
    rowsweep_seed = 20250729;
    // Defaults reproduce the original row sweep exactly. An empty tau spec
    // means "derive it from delta", as (1+delta)/4 always did.
    oracle_spec = "t1";
    rowsweep_tau_spec = "";
    rowsweep_tau2_spec = "";
    rowsweep_score_out = "";
    // tau well below the theorem's 1/2: the theorem's bar assumes a *perfect*
    // oracle on a fully-bad certificate, but the surrogate detects a known
    // four-cycle only ~0.673 of the time and a cluster that only approximates
    // the reticulation branch dilutes that further. Measured over 70 replicates
    // of n15+n25, J rises monotonically as tau falls: 0.670 at 0.50, 0.739 at
    // 0.25, 0.756 at 0.15, 0.763 at 0.10 (doc/NEW_ALGORITHMS.md section 3a).
    branch_cut_tau_spec = "0.15";
    branch_cut_resolution_z = 0.0;        // deprecated, not scale-free
    branch_cut_resolution_margin = 0.0;   // off by default
    branch_cut_propagate_tau = 0.0;       // off by default (no phase 2)
    branch_cut_cycles = 0;                // 0 = the older uniform random sampling
    branch_cut_corroborate = 0;           // 0/1 = no split-sample corroboration
    branch_cut_corroborate_frac = 1.0;    // all groups must fire
    branch_cut_corroborate_bar = 0.0;     // 0 = the inner bar is tau itself
    branch_cut_cmin = 0;                  // 0 = no absolute contradiction floor
    branch_cut_min_depth = 0;             // 0 = no targeted depth lift
    branch_cut_reuse_min_anchors = 1;     // 1 = lift on every edge
    branch_cut_trim = 0;                  // 0 = untrimmed, the published rule
    branch_cut_res_min_pooled = 0;        // 0 = no support guard on the margin
    branch_cut_min_support = 4;
    branch_cut_samples = 8;
    branch_cut_seed = 20250729;
    branch_cut_contrast_alpha = 0.0;      // 0 = off, absolute-tau rule only
    branch_cut_margin_q10 = 0.0;          // margin-dispersion rules, all off
    branch_cut_margin_sd = 0.0;
    branch_cut_margin_min = 0.0;
    branch_cut_margin_min_sets = 8;
    branch_cut_cluster_margin = 0.0;
    branch_cut_tau_low = 0.0;
    branch_cut_tau_high = 0.0;
    branch_cut_mlc_d = 0.0;
    branch_cut_mlc_t = 0.0;
    branch_cut_mlc_min_group = 5;
    branch_cut_fixed_streams = false;     // shared stream, as earlier runs used
    // The three defaults below changed on 2026-08-11. Each was measured free
    // or better and the three together are the best configuration tested
    // (n150, 50 reps: err 11.92 -> 11.58*, fp_blob 1.82 -> 1.50*, J +0.018*).
    // Use --branchcut-no-{cycle-reuse,walecki,anchor-spread} for the old paths.
    branch_cut_cycle_reuse = true;        // lift the h <= t_A cap (Lem. cycle-cover)
    branch_cut_walecki = true;            // real edge-disjoint cycles, not random orders
    branch_cut_anchor_spread = true;      // sample anchor pairs, do not enumerate
    branch_cut_anchor_quality = 0;        // uniform anchor sample, unweighted
    branch_cut_anchor_power = 1.0;
    branch_cut_anchor_quality_out = "";
    branch_cut_anchor_corroborate = 0;   // no anchor-level count floor
    branch_cut_m2_full = false;          // keep the conservative h <= t_A/2 at m=2
    branch_cut_anchor_rotate = false;    // one anchor pair per cycle, as before
    branch_cut_shared_coords = false;     // one sample per track, as before
    branch_cut_mode_spec = "";            // every track uses the cluster scan
    branch_cut_score_out = "";
    branch_cut_quad_out = "";
    // Corner-row resolution test: off by default, so every earlier corner-row
    // run reproduces. The sample budget is a multiple of the corner-size sum;
    // 0 would be exhaustive, which is |A1||A2||B1||B2| ~ n^4/256 per edge.
    corner_row_resolution_margin = 0.0;
    corner_row_resolution_samples = 20;
    rowsweep_heavy_spec = "1";
    rowsweep_anchors = 1;
    oracle_cf_max = 0.80;      // only consulted by a "cf" term
    oracle_margin = 0.0;       // only consulted by a "maj" term
    root_str = "";
    annotation_tree_file = "";
    //pvalue_file = "";

    normal_mode = "2";   // use best algorithm for normalizing based on artificial taxa
    execute_mode = "0";  // use fast algorithm
    taxa_mode = "0";     // don't use shared taxon data structure
    weight_mode = "n";   // don't use weighting
    score_mode = "0";    // don't score final tree
    data_mode = "t";     // input data are trees i.e. newick strings
    brln_mode = "g";     // estimate branch lengths under MSC for gene trees
    contract = false;
    char2tree = false;
    rootonly = false;
    pcsonly = false;
    quartet_format = "((___,___),(___,___));___";

    blob = false;
    store_pvalue = false;
    load_pvalue = false;
    enable_split_test = false;
    three_fix_one_alter = false;
    two_fix_two_alter = false;
    row_sweep_blob = false;
    corner_row_blob = false;
    branch_cut_blob = false;
    quard = false;
    network = false;

    override_file = false;

    support_low = 0.0;
    support_high = 1.0;
    support_default = 1.0;
    support_threshold = 0.0;

    alpha = 1e-7;
    beta = 0.95;

    refine_seed = 12345;
    cut_seed = 1;
    iter_limit = 10;
    iter_limit_blob = std::numeric_limits<unsigned long int>::max();

    dict = NULL;
    output = NULL;

    int outcome = parse(argc, argv);

    if (outcome == 1) {
        std::cout << help_info;
        exit(0);
    }

    if (outcome == 2) {
        std::cout << "use -h for more information" << std::endl;
        exit(1);
    }

    if (outcome > 2) {
        exit(1);
    }

    std::cout << "LOGGING:" << std::endl;

    // Read mapping file
    if (mapping_file != "") prepare_indiv2taxon_map();

    // Process input

    // First try to figure out data type and whether the file exists...

    dict = new Dict;
    if (data_mode == "t") input_trees();
    else if (data_mode == "q") input_quartets();
    else input_matrix();
    
    if (data_mode == "q") {
        if (quartets.size() == 0 && qCFs_table.size() == 0) {
            std::cout << "\nERROR: Nothing read from input" << std::endl;
            exit(1);
        }
    }
    else {
        if (input.size() == 0) {
            std::cout << "\nERROR: Nothing read from input" << std::endl;
            exit(1);
        }
    }

    dict->update_singletons();

    // if (pvalue_file != "") input_pvalues();

    if (char2tree) {
        std::cout << "Writing characters as trees" << std::endl;
        std::ofstream fout(output_file);
        if (fout.fail()) {
            for (Tree *t : input) std::cout << t->to_string_basic() << std::endl;
        }
        else {
            for (Tree *t : input) fout << t->to_string_basic() << std::endl;
            fout.close();
        }
        exit(0);
    }

    // TODO: Add suppress uniforcations!
    refine_trees();
    prepare_trees();

    if (verbose > "0") {
        subproblem_csv.open(output_file + "_subproblems.csv");
        subproblem_csv << "ID,PARENT,DEPTH,SIZE,ARTIFICIAL,SUBSET";
        if (verbose > "1") {
            subproblem_csv << ",ENTRY,PRUNED,ZERO";
        }
        subproblem_csv << std::endl;
    }

    // DONE SETTING UP!
}


Instance::~Instance() {
    for (Tree *t : input) delete t;
    delete output;
    delete dict;
}



#if ENABLE_TOB
// Split a comma-separated list of numbers, tolerating a single value that is
// then broadcast to every oracle track.
static std::vector<std::string> split_csv(const std::string &text) {
    std::vector<std::string> out;
    std::size_t pos = 0;
    while (pos <= text.size()) {
        const std::size_t comma = text.find(',', pos);
        out.push_back(text.substr(pos, comma == std::string::npos
                                      ? std::string::npos : comma - pos));
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    return out;
}

// Assemble the row sweep's configuration from the command line. With the
// defaults this is the original algorithm: oracle "t1", tau = (1+delta)/4,
// heavy = 1, anchors = 1.
RowSweepParams Instance::build_row_sweep_params(std::string *error) const {
    RowSweepParams params;
    params.anchors = rowsweep_anchors;
    params.oracle = OracleSpec::parse(oracle_spec, rowsweep_query_alpha,
                                      oracle_cf_max, oracle_margin, error);
    if (!error->empty()) return params;
    const std::size_t n_tracks = params.oracle.tracks.size();

    // An empty tau spec keeps the historical derivation from delta.
    std::vector<std::string> taus = rowsweep_tau_spec.empty()
        ? std::vector<std::string>(1, std::to_string((1.0 + rowsweep_delta) / 4.0))
        : split_csv(rowsweep_tau_spec);
    std::vector<std::string> heavies = split_csv(rowsweep_heavy_spec);

    if (taus.size() != 1 && taus.size() != n_tracks) {
        *error = "--rowsweep-tau needs one value or one per oracle track";
        return params;
    }
    if (heavies.size() != 1 && heavies.size() != n_tracks) {
        *error = "--rowsweep-heavy needs one value or one per oracle track";
        return params;
    }
    for (std::size_t t = 0; t < n_tracks; ++t) {
        weight_t tau;
        unsigned long int heavy;
        const std::string &tau_text = taus[taus.size() == 1 ? 0 : t];
        const std::string &heavy_text = heavies[heavies.size() == 1 ? 0 : t];
        if (!s2d(tau_text, &tau) || !(tau > 0.0 && tau < 1.0)) {
            *error = "row-sweep tau must lie strictly in (0, 1): " + tau_text;
            return params;
        }
        if (!s2ul(heavy_text, &heavy)) {
            *error = "row-sweep heavy must be a non-negative integer: " + heavy_text;
            return params;
        }
        params.tau.push_back(tau);
        params.heavy.push_back(heavy);
    }
    // Second-stage threshold; empty leaves the rule exactly one-stage.
    if (!rowsweep_tau2_spec.empty()) {
        std::vector<std::string> tau2s = split_csv(rowsweep_tau2_spec);
        if (tau2s.size() != 1 && tau2s.size() != n_tracks) {
            *error = "--rowsweep-tau2 needs one value or one per oracle track";
            return params;
        }
        for (std::size_t t = 0; t < n_tracks; ++t) {
            weight_t tau2;
            const std::string &text = tau2s[tau2s.size() == 1 ? 0 : t];
            if (!s2d(text, &tau2) || !(tau2 > 0.0 && tau2 < 1.0)) {
                *error = "row-sweep tau2 must lie strictly in (0, 1): " + text;
                return params;
            }
            params.tau2.push_back(tau2);
        }
    }
    params.score_out = rowsweep_score_out;
    params.seed = rowsweep_seed;
    if (rowsweep_row_mode_spec == "fixed" || rowsweep_row_mode_spec.empty()) {
        params.row_mode = RowSweepRowMode::Fixed;
    } else if (rowsweep_row_mode_spec == "random") {
        params.row_mode = RowSweepRowMode::Random;
    } else if (rowsweep_row_mode_spec == "pooled") {
        params.row_mode = RowSweepRowMode::Pooled;
    } else {
        *error = "--rowsweep-row-mode must be fixed, random or pooled: "
                 + rowsweep_row_mode_spec;
        return params;
    }
    return params;
}
#endif  // ENABLE_TOB

long long Instance::solve() {
    srand(cut_seed);

    std::string mode = normal_mode + execute_mode + taxa_mode + weight_mode;

    auto start = std::chrono::high_resolution_clock::now();

    // Build or read species tree
    if (load_pvalue) {
        #if ENABLE_TOB
            std::cout << "Loading species tree with p-values" << root_str << std::endl;
            output = new SpeciesTree(input[0], dict, alpha, beta, enable_split_test);
        #else
            std::cout << "TREE-QMC was not compiled with tree of blob options!" << std::endl;
            exit(1);
        #endif  // ENABLE_TOB
    } else if (network) {

        #if ENABLE_TOB
        


        std::cout << "Loading annotation tree/TOB from " << annotation_tree_file << std::endl;

        get_annotation_tree();
        
        std::cout << "Loading species tree with p-values" << root_str << std::endl;
        

        
        
        if (iter_limit_blob == std::numeric_limits<unsigned long int>::max()) {
                    iter_limit_blob = 2 * dict->size() * dict->size();
                    std::cout << "Setting blob iteration limit to 2*ntaxa^2 = " << iter_limit_blob << std::endl;
                }

        
        if (data_mode == "q") {
            Tree *hybrid_tree_output = new SpeciesTree(annotation_tree, dict, alpha, beta, qCFs_table, iter_limit_blob);
            output_net = new Network(hybrid_tree_output, dict);
        } else {
            Tree *hybrid_tree_output = new SpeciesTree(annotation_tree, dict, alpha, beta, input, iter_limit_blob);
            output_net = new Network(hybrid_tree_output, dict);
        }

        #endif // ENABLE_TOB

    } else {
        if (stree_file != "") {
            std::cout << "Loading species tree" << root_str << std::endl;
            output = new SpeciesTree(stree_file, dict);
        } else if (data_mode == "q") {
            output = new SpeciesTree(quartets, dict, mode, iter_limit, output_file);
        } else {
            output = new SpeciesTree(input, dict, mode, iter_limit, output_file);
        }
    }
    

    if (root_str != "") {
        std::cout << "Rooting species tree at " << root_str << std::endl;
        output->root_at_clade(outgroup_taxon_set);
    }

    std::cout << "Printing output tree/network:" << std::endl;
    if (!network) {
        std::cout << output->to_string_basic() << std::endl;
    } else {
        std::cout << output_net->to_string_basic() << std::endl;
    }
    

    #if ENABLE_TOB
    if (rowsweep_file != "") {
        if (output == NULL) {
            std::cout << "ERROR: --rowsweep requires a species/gene-tree SpeciesTree object; "
                      << "it cannot be combined with --network mode." << std::endl;
            exit(1);
        }
        std::cout << "Running row-sweep split experiment" << std::endl;
        std::unordered_map<std::string, index_t> name2index;
        for (index_t idx = 0; idx < dict->size(); idx++)
            name2index[dict->index2label(idx)] = idx;
        std::string rs_error;
        RowSweepParams rs_params = build_row_sweep_params(&rs_error);
        if (!rs_error.empty()) {
            std::cout << "\nERROR: " << rs_error << std::endl;
            exit(1);
        }
        output->run_split_experiment(input, name2index, rowsweep_file,
                                     rowsweep_out_file, rs_params);
    }
    #else
    if (rowsweep_file != "") {
        std::cout << "TREE-QMC was not compiled with tree of blob options!" << std::endl;
        exit(1);
    }
    #endif  // ENABLE_TOB

    #if ENABLE_TOB
    // Instrumentation only: record the evidence the row sweep would consume on
    // this refinement and stop, without constructing a tree of blobs.
    if (rowsweep_dump_prefix != "") {
        if (output == NULL) {
            std::cout << "ERROR: --rowsweep-dump requires a species tree "
                      << "(use --blobsearchonly)." << std::endl;
            exit(1);
        }
        std::cout << "Dumping row-sweep evidence to " << rowsweep_dump_prefix << std::endl;
        std::string dump_error;
        RowSweepParams dump_params = build_row_sweep_params(&dump_error);
        if (!dump_error.empty()) {
            std::cout << "\nERROR: " << dump_error << std::endl;
            exit(1);
        }
        output->dump_row_sweep_evidence(input, rowsweep_dump_prefix,
                                        dump_params.oracle,
                                        rowsweep_dump_all_targets,
                                        rowsweep_dump_all_anchors);
        auto dump_end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            dump_end - start).count();
    }
    #endif  // ENABLE_TOB

    if (!load_pvalue && (store_pvalue || blob)) {
        #if ENABLE_TOB
            if (row_sweep_blob) {
                std::string error;
                RowSweepParams row_sweep = build_row_sweep_params(&error);
                if (!error.empty()) {
                    std::cout << "\nERROR: " << error << std::endl;
                    exit(1);
                }
                SpeciesTree* row_sweep_tree = new SpeciesTree(
                    input, dict, output, row_sweep
                );
                std::cout << "Printing output tree with pvalues:" << std::endl;
                std::cout << output->to_string_pvalue() << std::endl;
                if (!store_pvalue && blob) {
                    delete output;
                    output = row_sweep_tree;
                } else {
                    delete row_sweep_tree;
                }
            } else if (branch_cut_blob) {
                BranchCutParams branch_cut;
                branch_cut.min_support = branch_cut_min_support;
                branch_cut.resolution_z = branch_cut_resolution_z;
                branch_cut.resolution_margin = branch_cut_resolution_margin;
                branch_cut.propagate_tau = branch_cut_propagate_tau;
                branch_cut.cycles = branch_cut_cycles;
                branch_cut.trim = branch_cut_trim;
                branch_cut.resolution_min_pooled = branch_cut_res_min_pooled;
                branch_cut.samples = branch_cut_samples;
                branch_cut.query_alpha = rowsweep_query_alpha;
                branch_cut.seed = branch_cut_seed;
                branch_cut.contrast_alpha = branch_cut_contrast_alpha;
                branch_cut.margin_q10 = branch_cut_margin_q10;
                branch_cut.margin_sd = branch_cut_margin_sd;
                branch_cut.margin_min = branch_cut_margin_min;
                branch_cut.margin_min_sets = branch_cut_margin_min_sets;
                branch_cut.cluster_margin = branch_cut_cluster_margin;
                branch_cut.tau_low = branch_cut_tau_low;
                branch_cut.tau_high = branch_cut_tau_high;
                branch_cut.mlc_d = branch_cut_mlc_d;
                branch_cut.mlc_t = branch_cut_mlc_t;
                branch_cut.mlc_min_group = branch_cut_mlc_min_group;
                branch_cut.fixed_streams = branch_cut_fixed_streams;
                branch_cut.cycle_reuse = branch_cut_cycle_reuse;
                branch_cut.walecki = branch_cut_walecki;
                branch_cut.anchor_spread = branch_cut_anchor_spread;
                branch_cut.anchor_quality = branch_cut_anchor_quality;
                branch_cut.anchor_power = branch_cut_anchor_power;
                branch_cut.anchor_quality_out = branch_cut_anchor_quality_out;
                branch_cut.anchor_corroborate = branch_cut_anchor_corroborate;
                branch_cut.m2_full = branch_cut_m2_full;
                branch_cut.anchor_rotate = branch_cut_anchor_rotate;
                branch_cut.shared_coords = branch_cut_shared_coords;
                branch_cut.score_out = branch_cut_score_out;
                branch_cut.quad_out = branch_cut_quad_out;
                branch_cut.corroborate = branch_cut_corroborate;
                branch_cut.corroborate_frac = branch_cut_corroborate_frac;
                branch_cut.corroborate_bar = branch_cut_corroborate_bar;
                branch_cut.cmin = branch_cut_cmin;
                branch_cut.min_depth = branch_cut_min_depth;
                branch_cut.reuse_min_anchors = branch_cut_reuse_min_anchors;
                std::string bc_error;
                branch_cut.oracle = OracleSpec::parse(oracle_spec,
                                                      rowsweep_query_alpha,
                                                      oracle_cf_max,
                                                      oracle_margin, &bc_error);
                if (!bc_error.empty()) {
                    std::cout << "\nERROR: " << bc_error << std::endl;
                    exit(1);
                }
                {
                    std::vector<std::string> texts = split_csv(branch_cut_tau_spec);
                    const std::size_t n_tracks = branch_cut.oracle.tracks.size();
                    if (texts.size() != 1 && texts.size() != n_tracks) {
                        std::cout << "\nERROR: --branchcut-tau needs one value "
                                  << "or one per oracle track" << std::endl;
                        exit(1);
                    }
                    for (std::size_t t = 0; t < n_tracks; ++t) {
                        weight_t tau;
                        const std::string &text = texts[texts.size() == 1 ? 0 : t];
                        if (!s2d(text, &tau) || !(tau > 0.0 && tau < 1.0)) {
                            std::cout << "\nERROR: branch-cut tau must lie in "
                                      << "(0, 1): " << text << std::endl;
                            exit(1);
                        }
                        branch_cut.tau.push_back(tau);
                    }
                    // Per-track rule: `cluster` (the scan) or `global` (one
                    // rate test on the whole side). Absent = all cluster.
                    if (!branch_cut_mode_spec.empty()) {
                        std::vector<std::string> modes = split_csv(branch_cut_mode_spec);
                        if (modes.size() != 1 && modes.size() != n_tracks) {
                            std::cout << "\nERROR: --branchcut-mode needs one value "
                                      << "or one per oracle track" << std::endl;
                            exit(1);
                        }
                        for (std::size_t t = 0; t < n_tracks; ++t) {
                            const std::string &m = modes[modes.size() == 1 ? 0 : t];
                            if (m == "cluster") branch_cut.mode.push_back(0);
                            else if (m == "global") branch_cut.mode.push_back(1);
                            else {
                                std::cout << "\nERROR: --branchcut-mode must be "
                                          << "`cluster` or `global`: " << m << std::endl;
                                exit(1);
                            }
                        }
                    }
                }
                SpeciesTree* branch_cut_tree = new SpeciesTree(
                    input, dict, output, branch_cut
                );
                std::cout << "Printing output tree with pvalues:" << std::endl;
                std::cout << output->to_string_pvalue() << std::endl;
                if (!store_pvalue && blob) {
                    delete output;
                    output = branch_cut_tree;
                } else {
                    delete branch_cut_tree;
                }
            } else if (corner_row_blob) {
                CornerRowParams corner_row;
                corner_row.k = corner_row_k;
                corner_row.heavy = corner_row_heavy;
                corner_row.query_alpha = rowsweep_query_alpha;
                {
                    // The oracle has to be parsed first: it fixes how many
                    // tracks there are, and so how many taus are wanted.
                    std::string cr_error;
                    corner_row.oracle = OracleSpec::parse(oracle_spec,
                                                          rowsweep_query_alpha,
                                                          oracle_cf_max,
                                                          oracle_margin, &cr_error);
                    if (!cr_error.empty()) {
                        std::cout << "\nERROR: " << cr_error << std::endl;
                        exit(1);
                    }
                }
                {
                    std::vector<std::string> texts = corner_row_tau_spec.empty()
                        ? std::vector<std::string>(1,
                              std::to_string((1.0 + rowsweep_delta) / 4.0))
                        : split_csv(corner_row_tau_spec);
                    const std::size_t n_tracks = corner_row.oracle.tracks.size();
                    if (texts.size() != 1 && texts.size() != n_tracks) {
                        std::cout << "\nERROR: --corner-tau needs one value or "
                                  << "one per oracle track" << std::endl;
                        exit(1);
                    }
                    for (std::size_t t = 0; t < n_tracks; ++t) {
                        weight_t tau;
                        const std::string &text = texts[texts.size() == 1 ? 0 : t];
                        if (!s2d(text, &tau) || !(tau > 0.0 && tau < 1.0)) {
                            std::cout << "\nERROR: corner-row tau must lie in "
                                      << "(0, 1): " << text << std::endl;
                            exit(1);
                        }
                        corner_row.tau.push_back(tau);
                    }
                }
                {
                    // One count per oracle track, like --corner-tau: a noisy
                    // indicator can be made to demand several bad rows while a
                    // reliable one still fires on the first.
                    const std::size_t n_tracks = corner_row.oracle.tracks.size();
                    std::vector<std::string> texts =
                        corner_row_corroborate_spec.empty()
                        ? std::vector<std::string>(1, "1")
                        : split_csv(corner_row_corroborate_spec);
                    if (texts.size() != 1 && texts.size() != n_tracks) {
                        std::cout << "\nERROR: --corner-corroborate needs one "
                                  << "value or one per oracle track" << std::endl;
                        exit(1);
                    }
                    for (std::size_t t = 0; t < n_tracks; ++t) {
                        unsigned long int value;
                        const std::string &text = texts[texts.size() == 1 ? 0 : t];
                        if (!s2ul(text, &value) || value == 0) {
                            std::cout << "\nERROR: --corner-corroborate must be "
                                      << "a positive integer: " << text << std::endl;
                            exit(1);
                        }
                        corner_row.corroborate.push_back(value);
                    }
                }
                corner_row.seed = corner_row_seed;
                corner_row.cross = corner_row_cross;
                corner_row.resolution_margin = corner_row_resolution_margin;
                corner_row.resolution_samples = corner_row_resolution_samples;
                SpeciesTree* corner_row_tree = new SpeciesTree(
                    input, dict, output, corner_row
                );
                std::cout << "Printing output tree with pvalues:" << std::endl;
                std::cout << output->to_string_pvalue() << std::endl;
                if (!store_pvalue && blob) {
                    delete output;
                    output = corner_row_tree;
                } else {
                    delete corner_row_tree;
                }
            } else {
                if (!three_fix_one_alter && !two_fix_two_alter) {
                    if (iter_limit_blob == std::numeric_limits<unsigned long int>::max()) {
                        iter_limit_blob = 2 * dict->size() * dict->size();
                        std::cout << "Setting blob iteration limit to 2*ntaxa^2 = " << iter_limit_blob << std::endl;
                    }
                }
                SpeciesTree* display = new SpeciesTree(input, dict, output, iter_limit_blob, three_fix_one_alter, two_fix_two_alter, quard, output_qcfs_table_file);
                delete display;
                std::cout << "Printing output tree with pvalues:" << std::endl;
                std::cout << output->to_string_pvalue() << std::endl;
            }
        #else
            std::cout << "TREE-QMC was not compiled with tree of blob options!" << std::endl;
            exit(1);
        #endif  // ENABLE_TOB
    }

    #if ENABLE_TOB
    if (!load_pvalue && !store_pvalue && blob && !row_sweep_blob && !corner_row_blob && !branch_cut_blob) {
        SpeciesTree* tmp = new SpeciesTree(output, dict, alpha, beta, enable_split_test);
        delete output;
        output = tmp;
    }
    #endif  // ENABLE_TOB
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    return duration.count();
}


SpeciesTree *Instance::get_solution() {
    return output;
}

Network *Instance::get_network_solution() {
    return output_net;
}

void Instance::output_solution() {
    if (execute_mode == "2" || execute_mode == "3") return;

    // compute pcs for specified branch, then exit
    if (pcsonly) {
        std::cout << "Writing PCS" << std::endl;

        if (output_file != "") {
            std::ifstream fin(output_file);
            if (!fin.fail()) {
                //std::cout << override_file << std::endl;
                if (!override_file) {
                    std::cout << "  WARNING: " << output_file << " already exists, writing to stdout" << std::endl;
                    output_file = "";
                } else {
                    std::cout << "  WARNING: " << output_file << " already exists, overriding" << std::endl; 
                }
            }
            fin.close();
        }

        if (output_file != "") {
            std::ofstream fout(output_file);
            if (!fout.fail()) {
                output->write_pcs_table(input, positions, weight_mode, fout);
                fout.close();
            }
            else {
                std::cout << "  WARNING: Unable to write to " << output_file << ", writing to stdout" << std::endl; 
                output_file = "";
            }
        }

        if (output_file == "") output->write_pcs_table(input, positions, weight_mode, std::cout);

        return;
    }

    // if not pcs, write species tree
    if (score_mode == "1") {
        std::cout << "Computing branch info for species tree" << std::endl;
        output->refine();
        if (data_mode == "q") 
            output->annotate(quartets, weight_mode);
        else 
            output->annotate(input, weight_mode);
    }

    if (table_file != "") {
        std::cout << "Writing support to table" << std::endl;
        std::ofstream fout(table_file);
        if (fout.fail()) {
            std::cout << "  WARNING: Unable to write to " << table_file << std::endl; 
        } else {
            output->write_support_table(fout, brln_mode);
            fout.close();
        }
    }

    std::cout << "Writing species tree" << std::endl;
    std::ifstream fin(output_file);
    if (!fin.fail()) {
        if (!override_file) {
            std::cout << "  WARNING: " << output_file << " already exists, writing to stdout" << std::endl;
            output_file = "";
        } else {
            std::cout << "  WARNING: " << output_file << " already exists, overriding" << std::endl;
        }
    }
    fin.close();

    if (output_file != "") {
        std::ofstream fout(output_file);
        if (!fout.fail()) {
            if (score_mode == "1")
                fout << output->to_string_annotated(brln_mode) << std::endl;
            #if ENABLE_TOB
            else if (store_pvalue) 
                fout << output->to_string_pvalue() << std::endl;
            else if (network)
                fout << output_net->to_string_basic() << std::endl;
            #endif  // ENABLE_TOB
            else 
                fout << output->to_string_basic() << std::endl;
            fout.close();
            return;
        }
        std::cout << "  WARNING: Unable to write to " << output_file << ", writing to stdout" << std::endl;
    }

    if (score_mode == "1")
        std::cout << output->to_string_annotated(brln_mode) << std::endl;
    #if ENABLE_TOB
    else if (store_pvalue) 
        std::cout << output->to_string_pvalue() << std::endl;
    else if (network)
        std::cout << output_net->to_string_basic() << std::endl;
    #endif  // ENABLE_TOB
    else
        std::cout << output->to_string_basic() << std::endl;
}


int Instance::parse(int argc, char **argv) {
    std::cout << "TREE-QMC version " << VERSION << std::endl;

    std::cout << "COMMAND: ";
    for (int j = 0; j < argc; j ++) 
        std::cout << argv[j] << " ";
    std::cout << std::endl;

    int nweightparam = 0;
    int nminparam = 0;
    int nmaxparam = 0;
    int ndefaultparam = 0;
    bool help = false;

    index_t i = 1;
    while (i < argc) {
        std::string opt(argv[i]);
        if (opt == "-h" || opt == "--help") {
            return 1;
        }
        else if (opt == "-i" || opt == "--input") {
            if (i < argc - 1) input_file = argv[++ i];
        }
        else if (opt == "--chars") {
            data_mode = "c";  // input data are multi-state characters
            brln_mode = "n";  // don't estimate branch lengths!
        }
        else if (opt == "--bp") {
            data_mode = "b";  // input data are 2-state characters
            brln_mode = "b";  // estimate branch lengths under nWF+IS model
        }
        else if (opt == "--quartets") {
            data_mode = "q";
        }
        else if (opt == "--quartetformat") {
            if (i < argc - 1) {
                quartet_format = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No quartet format specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "-a" || opt == "--mapping") {
            if (i < argc - 1) {
                mapping_file = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No mapping file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "-u" || opt == "--support") {
            score_mode = "1";
        }
        else if (opt == "--override") {
            override_file = true;
        }
        else if (opt == "--blob") {
            blob = true;
        }
        else if (opt == "--network") {
            network = true;
        }
        else if (opt == "--3f1a") {
            three_fix_one_alter = true;
        }
        else if (opt == "--2f2a") {
            two_fix_two_alter = true;
        }
        else if (opt == "--rowsweep-blob") {
            row_sweep_blob = true;
        }
        else if (opt == "--oracle") {
            if (i < argc - 1) {
                oracle_spec = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No oracle spec specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--oracle-cf-max") {
            if (i < argc - 1) {
                if (!s2d(argv[++ i], &oracle_cf_max) || oracle_cf_max <= 0.0) {
                    std::cout << "\nERROR: oracle cf-max must be positive" << std::endl;
                    return 2;
                }
            }
            else {
                std::cout << "\nERROR: No cf-max specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--oracle-margin") {
            if (i < argc - 1) {
                if (!s2d(argv[++ i], &oracle_margin) || oracle_margin < 0.0
                        || oracle_margin >= 1.0) {
                    std::cout << "\nERROR: oracle margin must lie in [0, 1)" << std::endl;
                    return 2;
                }
            }
            else {
                std::cout << "\nERROR: No margin specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-tau") {
            if (i < argc - 1) {
                rowsweep_tau_spec = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No tau specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-tau2") {
            if (i < argc - 1) {
                rowsweep_tau2_spec = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No tau specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-score-out") {
            if (i < argc - 1) {
                rowsweep_score_out = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No score output file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-heavy") {
            if (i < argc - 1) {
                rowsweep_heavy_spec = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No heavy specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-anchors") {
            if (i < argc - 1) {
                if (!s2ul(argv[++ i], &rowsweep_anchors)) {
                    std::cout << "\nERROR: invalid anchor count" << std::endl;
                    return 2;
                }
            }
            else {
                std::cout << "\nERROR: No anchor count specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--cornerrow-blob") {
            corner_row_blob = true;
        }
        else if (opt == "--branchcut-blob") {
            branch_cut_blob = true;
        }
        else if (opt == "--branchcut-tau") {
            if (i < argc - 1) {
                branch_cut_tau_spec = argv[++ i];
            } else {
                std::cout << "\nERROR: --branchcut-tau needs a value" << std::endl;
                return 2;
            }
        }
        else if (opt == "--branchcut-resolution-margin") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_resolution_margin)) {
                std::cout << "\nERROR: --branchcut-resolution-margin needs a value" << std::endl;
                return 2;
            }
        }
        else if (opt == "--branchcut-resolution-z") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_resolution_z)) {
                std::cout << "\nERROR: --branchcut-resolution-z needs a value" << std::endl;
                return 2;
            }
        }
        else if (opt == "--branchcut-min-support") {
            if (i >= argc - 1 || !s2ul(argv[++ i], &branch_cut_min_support)) {
                std::cout << "\nERROR: --branchcut-min-support needs an integer" << std::endl;
                return 2;
            }
        }
        else if (opt == "--branchcut-samples") {
            if (i >= argc - 1 || !s2ul(argv[++ i], &branch_cut_samples)) {
                std::cout << "\nERROR: --branchcut-samples needs an integer" << std::endl;
                return 2;
            }
        }
        else if (opt == "--branchcut-seed") {
            if (i >= argc - 1 || !s2ul(argv[++ i], &branch_cut_seed)) {
                std::cout << "\nERROR: --branchcut-seed needs an integer" << std::endl;
                return 2;
            }
        }
        else if (opt == "--quard") {
            quard = true;
        }
        else if (opt == "--write_qcfs") {
            if (i < argc - 1) {
                output_qcfs_table_file = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No output qCF table file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--alpha") {
            alpha = std::stod(argv[++ i]);
        }
        else if (opt == "--beta") {
            beta = std::stod(argv[++ i]);
        }
        else if (opt == "--smt") {
            enable_split_test = true;
        }
        /*else if (opt == "--pvalue") {
            if (i < argc - 1) {
                pvalue_file = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No pvalue file specified" << std::endl;
                return 2;
            }
        }*/
        else if (opt == "--store_pvalue") {
            store_pvalue = true;
        }
        else if (opt == "--load_pvalue") {
            load_pvalue = true;
        }
        else if (opt == "--iter_limit_blob") {
            if (i < argc - 1) {
                iter_limit_blob = std::stoi(argv[++ i]);
            }
        }
        else if (opt == "--at") {
            if (i < argc - 1) {
                annotation_tree_file = argv[++ i];
            }
            std::cout << "Loading annotation tree/TOB from " << annotation_tree_file << std::endl;
            
        } else if (opt == "-q" || opt == "--supportonly") {
            score_mode = "1";
            if (i < argc - 1) {
                stree_file = argv[++ i];
            }
            else  {
                std::cout << "\nERROR: No species tree file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--blobsearchonly") {
            if (i < argc - 1) {
                stree_file = argv[++ i];
                std::cout << "Need to read " << stree_file << std::endl;
            }
            else  {
                std::cout << "\nERROR: No species tree file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "-r" || opt == "--rootonly") {
            rootonly = true;
            if (i < argc - 1) {
                stree_file = argv[++ i];
            }
            else  {
                std::cout << "\nERROR: No species tree file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--pcsonly") {
            pcsonly = true;
            if (i < argc - 1) {
                stree_file = argv[++ i];
            }
            else  {
                std::cout << "\nERROR: No species tree file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "-o" || opt == "--output") {
            if (i < argc - 1) output_file = argv[++ i];
        }
        else if (opt == "--root") {
            if (i < argc - 1) {
                root_str = argv[++ i];
            }
            else  {
                std::cout << "\nERROR: No root specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--char2tree") {
            data_mode = "c";  // input data are multi-state characters
            char2tree = true;
        }
        else if (opt == "--writetable") {
            if (i < argc - 1) {
                table_file = argv[++ i];
            }
            else  {
                std::cout << "\nERROR: No table file specified" << std::endl;
                return 2;
            }
        }

        // Handle weighting options
        else if (opt == "--hybrid") {
            weight_mode = 'h';
            nweightparam += 1;
        }
        else if (opt == "--fast") {
            weight_mode = 'f';
            nweightparam += 1;
        }
        else if (opt == "-w" || opt == "--weight") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param == "1" || param == "h") weight_mode = "h";
            else if (param == "2" || param == "s") weight_mode = "s";
            else if (param == "3" || param == "l") weight_mode = "l";
            else if (param == "4" || param == "n") weight_mode = "n";
            else if (param == "5" || param == "f") weight_mode = "f";
            else {
                std::cout << "\nERROR: invalid (weight) mode parameter: " << param << std::endl;
                return 2;
            }
            nweightparam += 1;
        }

        // Handle support range options
        else if (opt == "-B" || opt == "--bayes") {
            support_default = 0.333;
            support_low = 0.333;
            support_high = 1.0;
            ndefaultparam += 1;
            nminparam += 1;
            nmaxparam += 1;
        }
        else if (opt == "-L" || opt == "--lrt") {
            support_default = 0.0;
            support_low = 0.0;
            support_high = 1.0;
            ndefaultparam += 1;
            nminparam += 1;
            nmaxparam += 1;
        }
        else if (opt == "-S" || opt == "--bootstrap") {
            support_default = 0.0;
            support_low = 0.0;
            support_high = 100.0;
            ndefaultparam += 1;
            nminparam += 1;
            nmaxparam += 1;
        }
        else if (opt == "-d" || opt == "--default") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_default)) {
                std::cout << "\nERROR: invalid default support value: " << param << std::endl;
                return 2;
            }
            ndefaultparam += 1;
        }
        else if (opt == "-n" || opt == "--min") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_low)) {
                std::cout << "\nERROR: invalid min support value: " << param << std::endl;
                return 2;
            }
            nminparam += 1;
        }
        else if (opt == "-x" || opt == "--max") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_high)) {
                std::cout << "\nERROR: invalid max support value: " << param << std::endl;
                return 2;
            }
            nmaxparam += 1;
        }

        // Handle contraction option
        else if (opt == "--contract") {
            contract = true;
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2d(param, &support_threshold)) {
                std::cout << "\nERROR: invalid parameter: " << param << std::endl;
                return 2;
            }
        }

        // Handle artificial taxon normalization options
        else if (opt == "--norm_atax") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "\nERROR: invalid graph normalization parameter: " << param << std::endl;
                return 2;
            }
            normal_mode = param;
        }
        else if (opt == "--shared") {
            taxa_mode = "1";
        }

        // Handle other options
        else if (opt == "-v" || opt == "--verbose") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2") {
                std::cout << "\nERROR: invalid verbose parameter " << param << std::endl;
                return 2;
            }
            verbose = param;
        }
        else if (opt == "-e" || opt == "--execute") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (param != "0" && param != "1" && param != "2" && param != "3" && param != "4") {
                std::cout << "\nERROR: invalid execution parameter " << param << std::endl;
                return 2;
            }
            execute_mode = param;
        }
        else if (opt == "--polyseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &refine_seed)) {
                std::cout << "\nERROR: invalid polyseed parameter: " << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--cutseed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &cut_seed)) {
                std::cout << "\nERROR: invalid cutseed parameter: " << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--iterlimit") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &iter_limit)) {
                std::cout << "\nERROR: invalid parameter: " << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep") {
            if (i < argc - 1) {
                rowsweep_file = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No bipartition file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep_out") {
            if (i < argc - 1) {
                rowsweep_out_file = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No rowsweep output file specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-dump") {
            if (i < argc - 1) {
                rowsweep_dump_prefix = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No rowsweep dump prefix specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-dump-all-targets") {
            rowsweep_dump_all_targets = true;
        }
        else if (opt == "--rowsweep-dump-all-anchors") {
            rowsweep_dump_all_anchors = true;
        }
        else if (opt == "--delta") {
            if (i < argc - 1) {
                rowsweep_delta = std::stod(argv[++ i]);
            }
            else {
                std::cout << "\nERROR: No delta specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--corner-k") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &corner_row_k)) {
                std::cout << "\nERROR: invalid parameter: " << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--heavy-sampling") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &corner_row_heavy) || corner_row_heavy == 0) {
                std::cout << "\nERROR: --heavy-sampling must be a positive integer: "
                          << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--corner-tau") {
            if (i < argc - 1) {
                corner_row_tau_spec = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No corner-row tau specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--corner-corroborate") {
            if (i < argc - 1) {
                corner_row_corroborate_spec = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No corner-row corroboration specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-row-mode") {
            if (i < argc - 1) {
                rowsweep_row_mode_spec = argv[++ i];
            }
            else {
                std::cout << "\nERROR: No row-sweep row mode specified" << std::endl;
                return 2;
            }
        }
        else if (opt == "--rowsweep-seed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &rowsweep_seed)) {
                std::cout << "\nERROR: invalid parameter: " << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--corner-cross") {
            corner_row_cross = true;
        }
        else if (opt == "--branchcut-trim") {
            if (i + 1 < argc) {
                if (!s2ul(argv[++ i], &branch_cut_trim)) {
                    std::cout << "\nERROR: invalid --branchcut-trim" << std::endl; return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-trim needs a value" << std::endl; return 1; }
        }
        else if (opt == "--resolution-min-pooled") {
            if (i + 1 < argc) {
                if (!s2ul(argv[++ i], &branch_cut_res_min_pooled)) {
                    std::cout << "\nERROR: invalid --resolution-min-pooled" << std::endl; return 1;
                }
            } else { std::cout << "\nERROR: --resolution-min-pooled needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-cycles") {
            if (i + 1 < argc) {
                if (!s2ul(argv[++ i], &branch_cut_cycles)) {
                    std::cout << "\nERROR: invalid --branchcut-cycles" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-cycles needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-contrast") {
            if (i < argc - 1) {
                if (!s2d(argv[++ i], &branch_cut_contrast_alpha)
                    || branch_cut_contrast_alpha < 0.0
                    || branch_cut_contrast_alpha > 1.0) {
                    std::cout << "\nERROR: --branchcut-contrast must lie in [0, 1]" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-contrast needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-margin-q10") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_margin_q10)) {
                std::cout << "\nERROR: --branchcut-margin-q10 needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-margin-sd") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_margin_sd)) {
                std::cout << "\nERROR: --branchcut-margin-sd needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-margin-min") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_margin_min)) {
                std::cout << "\nERROR: --branchcut-margin-min needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-margin-min-sets") {
            if (i >= argc - 1 || !s2ul(argv[++ i], &branch_cut_margin_min_sets)) {
                std::cout << "\nERROR: --branchcut-margin-min-sets needs an integer" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-cluster-margin") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_cluster_margin)) {
                std::cout << "\nERROR: --branchcut-cluster-margin needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-tau-high") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_tau_high)) {
                std::cout << "\nERROR: --branchcut-tau-high needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-tau-low") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_tau_low)) {
                std::cout << "\nERROR: --branchcut-tau-low needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-mlc-d") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_mlc_d)) {
                std::cout << "\nERROR: --branchcut-mlc-d needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-mlc-t") {
            if (i >= argc - 1 || !s2d(argv[++ i], &branch_cut_mlc_t)) {
                std::cout << "\nERROR: --branchcut-mlc-t needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-mlc-min-group") {
            if (i >= argc - 1 || !s2ul(argv[++ i], &branch_cut_mlc_min_group)) {
                std::cout << "\nERROR: --branchcut-mlc-min-group needs an integer" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-mode") {
            if (i < argc - 1) branch_cut_mode_spec = argv[++ i];
            else { std::cout << "\nERROR: --branchcut-mode needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-fixed-streams") {
            branch_cut_fixed_streams = true;
        }
        else if (opt == "--branchcut-cycle-reuse") {
            branch_cut_cycle_reuse = true;
        }
        else if (opt == "--branchcut-no-cycle-reuse") {
            branch_cut_cycle_reuse = false;
        }
        else if (opt == "--branchcut-walecki") {
            branch_cut_walecki = true;
        }
        else if (opt == "--branchcut-no-walecki") {
            branch_cut_walecki = false;
        }
        else if (opt == "--branchcut-anchor-spread") {
            branch_cut_anchor_spread = true;
        }
        else if (opt == "--branchcut-anchor-power") {
            if (i < argc - 1) branch_cut_anchor_power = std::stod(argv[++ i]);
        }
        else if (opt == "--branchcut-anchor-quality-out") {
            if (i < argc - 1) branch_cut_anchor_quality_out = argv[++ i];
        }
        else if (opt == "--branchcut-anchor-rotate") {
            branch_cut_anchor_rotate = true;
            branch_cut_anchor_spread = true;
        }
        else if (opt == "--branchcut-m2-full") {
            branch_cut_m2_full = true;
        }
        else if (opt == "--branchcut-anchor-corroborate") {
            branch_cut_anchor_corroborate = 2;
            if (i < argc - 1 && argv[i + 1][0] != '-')
                branch_cut_anchor_corroborate = std::stoul(argv[++ i]);
        }
        else if (opt == "--branchcut-no-anchor-spread") {
            branch_cut_anchor_spread = false;
        }
        else if (opt == "--branchcut-anchor-quality") {
            branch_cut_anchor_quality = 4000;
            if (i < argc - 1 && argv[i + 1][0] != '-')
                branch_cut_anchor_quality = std::stoul(argv[++ i]);
            branch_cut_anchor_spread = true;
        }
        else if (opt == "--branchcut-shared-coords") {
            branch_cut_shared_coords = true;
            branch_cut_fixed_streams = true;
        }
        else if (opt == "--branchcut-quad-out") {
            if (i < argc - 1) branch_cut_quad_out = argv[++ i];
            else { std::cout << "\nERROR: --branchcut-quad-out needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-reuse-min-anchors") {
            if (i < argc - 1) {
                if (!s2ul(argv[++ i], &branch_cut_reuse_min_anchors)) {
                    std::cout << "\nERROR: invalid --branchcut-reuse-min-anchors" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-reuse-min-anchors needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-min-depth") {
            if (i < argc - 1) {
                if (!s2ul(argv[++ i], &branch_cut_min_depth)) {
                    std::cout << "\nERROR: invalid --branchcut-min-depth" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-min-depth needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-cmin") {
            if (i < argc - 1) {
                if (!s2ul(argv[++ i], &branch_cut_cmin)) {
                    std::cout << "\nERROR: invalid --branchcut-cmin" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-cmin needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-corroborate-bar") {
            if (i < argc - 1) {
                if (!s2d(argv[++ i], &branch_cut_corroborate_bar)
                    || branch_cut_corroborate_bar < 0.0
                    || branch_cut_corroborate_bar >= 1.0) {
                    std::cout << "\nERROR: --branchcut-corroborate-bar must lie in [0, 1)" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-corroborate-bar needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-corroborate-frac") {
            if (i < argc - 1) {
                if (!s2d(argv[++ i], &branch_cut_corroborate_frac)
                    || branch_cut_corroborate_frac <= 0.0
                    || branch_cut_corroborate_frac > 1.0) {
                    std::cout << "\nERROR: --branchcut-corroborate-frac must lie in (0, 1]" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-corroborate-frac needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-corroborate") {
            if (i < argc - 1) {
                if (!s2ul(argv[++ i], &branch_cut_corroborate)) {
                    std::cout << "\nERROR: invalid --branchcut-corroborate" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-corroborate needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-score-out") {
            if (i < argc - 1) branch_cut_score_out = argv[++ i];
            else { std::cout << "\nERROR: --branchcut-score-out needs a value" << std::endl; return 1; }
        }
        else if (opt == "--branchcut-propagate-tau") {
            if (i + 1 < argc) {
                if (!s2d(argv[++ i], &branch_cut_propagate_tau)
                    || branch_cut_propagate_tau < 0.0
                    || branch_cut_propagate_tau >= 1.0) {
                    std::cout << "\nERROR: --branchcut-propagate-tau must lie in [0, 1)" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --branchcut-propagate-tau needs a value" << std::endl; return 1; }
        }
        else if (opt == "--cornerrow-resolution-margin") {
            if (i + 1 < argc) {
                if (!s2d(argv[++ i], &corner_row_resolution_margin)) {
                    std::cout << "\nERROR: invalid --cornerrow-resolution-margin" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --cornerrow-resolution-margin needs a value" << std::endl; return 1; }
        }
        else if (opt == "--resolution-samples") {
            if (i + 1 < argc) {
                if (!s2ul(argv[++ i], &corner_row_resolution_samples)) {
                    std::cout << "\nERROR: invalid --resolution-samples" << std::endl;
                    return 1;
                }
            } else { std::cout << "\nERROR: --resolution-samples needs a value" << std::endl; return 1; }
        }
        else if (opt == "--corner-seed") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (! s2ul(param, &corner_row_seed)) {
                std::cout << "\nERROR: invalid parameter: " << param << std::endl;
                return 2;
            }
        }
        else if (opt == "--query-alpha") {
            std::string param = "";
            if (i < argc - 1) param = argv[++ i];
            if (!s2d(param, &rowsweep_query_alpha) || !std::isfinite(rowsweep_query_alpha) ||
                    rowsweep_query_alpha < 0.0 || rowsweep_query_alpha > 1.0) {
                std::cout << "\nERROR: --query-alpha must be between 0 and 1: "
                          << param << std::endl;
                return 2;
            }
        }
        else {
            std::cout << "ERROR: Unrecognized option: " << opt << std::endl;
            exit(1);
        }

        i ++;
    }

    if (input_file == "") {
        std::cout << "\nERROR: No input file specified" << std::endl;
        return 2;
    }

    if (row_sweep_blob && (three_fix_one_alter || two_fix_two_alter || quard)) {
        std::cout << "\nERROR: --rowsweep-blob cannot be combined with "
                  << "--3f1a, --2f2a, or --quard" << std::endl;
        return 2;
    }

    if (corner_row_blob && (three_fix_one_alter || two_fix_two_alter || quard)) {
        std::cout << "\nERROR: --cornerrow-blob cannot be combined with "
                  << "--3f1a, --2f2a, or --quard" << std::endl;
        return 2;
    }

    if (corner_row_blob && row_sweep_blob) {
        std::cout << "\nERROR: --cornerrow-blob cannot be combined with "
                  << "--rowsweep-blob" << std::endl;
        return 2;
    }

    if (branch_cut_blob && (row_sweep_blob || corner_row_blob
                            || three_fix_one_alter || two_fix_two_alter || quard)) {
        std::cout << "\nERROR: --branchcut-blob cannot be combined with "
                  << "--rowsweep-blob, --cornerrow-blob, --3f1a, --2f2a "
                  << "or --quard" << std::endl;
        return 2;
    }

    // Report settings and check that they make sense
    std::cout << std::endl << "SETTINGS:" << std::endl;

    // Input files
    std::cout << "input file: " << input_file << std::endl;
    if (mapping_file != "") std::cout << "mapping file: " << mapping_file << std::endl;
    if (stree_file != "") std::cout << "species tree file: " << stree_file << std::endl;
    if (table_file != "") std::cout << "table file: " << table_file << std::endl;
    if (output_qcfs_table_file != "")
        std::cout << "qCF table file: " << output_qcfs_table_file << std::endl;
    if (rowsweep_file != "" || row_sweep_blob) {
        std::cout << "row-sweep query alpha: " << rowsweep_query_alpha << std::endl;
    }
    if (corner_row_blob) {
        std::cout << "corner-row query alpha: " << rowsweep_query_alpha << std::endl;
        std::cout << "corner-row tau: ";
        if (corner_row_tau_spec.empty())
            std::cout << (1.0 + rowsweep_delta) / 4.0 << " (from delta = "
                      << rowsweep_delta << ")" << std::endl;
        else
            std::cout << corner_row_tau_spec << std::endl;
        std::cout << "corner-row k: ";
        if (corner_row_k == 0) std::cout << "n-3" << std::endl;
        else std::cout << corner_row_k << std::endl;
        std::cout << "corner-row heavy sampling: " << corner_row_heavy
                  << " (each row draws min(" << corner_row_heavy
                  << "k, row population))" << std::endl;
        std::cout << "corner-row seed: " << corner_row_seed << std::endl;
        std::cout << "corner-row row partners: "
                  << (corner_row_cross ? "opposite corner only (--corner-cross)"
                                       : "the whole near side")
                  << std::endl;
    }

    // Output file
    if (output_file == "")
         std::cout << "output file: std" << std::endl;
    else
        std::cout << "output file: " << output_file << std::endl;

    // Get outgroup taxon set for rooting
    if (root_str != "") {
        prepare_root_taxa();
        std::cout << "target root placement:";
        for (auto &taxon : outgroup_taxon_set) std::cout << " " << taxon;
        std::cout << std::endl;
    }

    if (data_mode == "t") {
        // Process weighting options
        if (nweightparam > 1) {
            std::cout << "\nERROR: Conflicting weighting options specified" << std::endl;
            return 2;
        }
        if (weight_mode == "s")
            std::cout << "weighting mode: support only" << std::endl;
        else if (weight_mode == "h")
            std::cout << "weighting mode: hybrid" << std::endl;
        else if (weight_mode == "l")
            std::cout << "weighting mode: length only" << std::endl;
        else if (weight_mode == "f")
            std::cout << "weighting mode: none (fast)" << std::endl;
        else if (weight_mode == "n")
            std::cout << "weighting mode: none" << std::endl;
        else {
            std::cout << "\nERROR: Unrecognized weight option!" << std::endl;
            return 2;
        }
        if (weight_mode != "h")
            std::cout << "  WARNING: --hybrid option is recommended" << std::endl;

        // Process support branch options
        if (! load_pvalue && (weight_mode == "s" || weight_mode == "h" || contract)) {
            // Check support options make sense
            if (nminparam  == 0 || nmaxparam == 0 || ndefaultparam == 0) {
                std::cout << "\nERROR: Must specify min, max, and default support values or use preset option" << std::endl;
                return 2;
            }
            else if (nminparam > 1 || nmaxparam > 1 || ndefaultparam > 1) {
                std::cout << "\nERROR: Multiple support parameters specified" << std::endl;
                return 2;
            }

            std::cout << "support min: " << (double)support_low << std::endl;
            std::cout << "support max: " << (double)support_high << std::endl;
            std::cout << "support default: " << (double)support_default << std::endl;

            // Check support parameters make sense 
            if (support_low >= support_high || support_high <= support_low) {
                std::cout << "\nERROR: Conflicting support values specified" << std::endl;
                return 3;
            }
            if (support_low < 0) {
                std::cout << "\nERROR: Support min value must be non-negative" << std::endl;
                return 3;
            }
            if (support_default < support_low || support_default > support_high) {
                std::cout << "\nERROR: Default support value is outside of range" << std::endl;
                return 3;
            }
        } else {
            // Ignore any support values that were specified and set back to defaults
            if (nminparam > 0 || nmaxparam > 0 || ndefaultparam > 0)
                std::cout << "  WARNING: Running in unweighted mode (no contraction) so ignoring support options" << std::endl;
            support_low = 0.0;
            support_high = 1.0;
            support_default = 1.0;
        }

        // Process contraction options
        if (contract) {
            if (weight_mode == "f") {
                std::cout << "\nERROR: Cannot contract low support branches in --fast mode" << std::endl;
                return 2;
            }

            std::cout << "contract support threshold: " << (double)support_threshold << " ("
                      << (double)(support_low + ((support_high - support_low) * support_threshold)) 
                      << " pre-transform)" << std::endl;

            if (support_threshold < 0 || support_threshold > 1) {
                std::cout << "\nERROR: Support threshold must be between 0 and 1 (it's applied after support values are mapped to this interval)" << std::endl;
                return 1;
            }
        }
    }
    else {
        if (nweightparam > 0 || nminparam > 0 || nmaxparam > 0 || ndefaultparam > 0) {
            if (data_mode == "b")
                std::cout << "  WARNING: Running in bipartition mode so ignoring any weight or support options" << std::endl;
            else
                std::cout << "  WARNING: Running in character mode so ignoring any weight or support options" << std::endl;
        }
        weight_mode = "n";
        support_low = 0.0;
        support_high = 1.0;
        support_default = 1.0;
        support_threshold = 0.0;
        contract = false;
    }

    // Print options for tree building
    if (stree_file == "") {
        // Process normalization mode
        std::cout << "normalization for artifical taxa mode: n" + normal_mode;
        if (taxa_mode == "1") std::cout << " (shared)";
        std::cout << std::endl;
        if (normal_mode != "2")
            std::cout << "  WARNING: --norm_atax 2 is recommended" << std::endl;

        // Process execution mode
        if (execute_mode == "0") {
            std::cout << "execution mode: efficient" << std::endl;
        }
        else if (execute_mode == "1") {
            std::cout << "execution mode: brute force" << std::endl;
        }
        else if (execute_mode == "2") {
            std::cout << "execution mode: compute weighted quartets, then exit" << std::endl;
            std::cout << "good edges will be saved in: " << output_file << "_quartets.txt" << std::endl;
        }
        else if (execute_mode == "3") {
            std::cout << "execution mode: compute good and bad edges, then exit" << std::endl;
            std::cout << "good edges will be saved in: " << output_file << "_good_edges.txt" << std::endl;
            std::cout << "bad edges will be saved in: " << output_file << "_bad_edges.txt" << std::endl;
        }
        else {
            DEBUG_MODE = true;
            execute_mode = "0";
            std::cout << "execution mode: efficient with brute force validation" << std::endl;
        }

        std::cout << "random seed for refinement: " << refine_seed << std::endl;
        std::cout << "random seed for max-cut: " << cut_seed << std::endl;
        std::cout << "max-cut heuristic iteration limit: " << iter_limit << std::endl;
    }

    std::cout << std::endl;

    return 0;
}


std::string Instance::get_execution_mode() {
    return execute_mode;
}


void Instance::input_trees() {
    std::ifstream fin(input_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << input_file << std::endl;
        exit(1);
    }

    std::cout << "Reading input" << std::endl;
    std::string newick;
    int mintax = INDEX_WIDTH;
    int maxtax = 0;
    std::size_t pos = 1;
    while (std::getline(fin, newick)) {
        // TODO: change to function that checks if newick string is valid, before proceeding
        if (newick.find(";") != std::string::npos) {
            Tree *t = new Tree(newick, dict, indiv2taxon, support_low, support_default);
            if (t->size() > maxtax) maxtax = t->size();
            if (t->size() < mintax) mintax = t->size();
            if (t->size() > 3) {
                input.push_back(t);
                if (pcsonly) positions.push_back(pos);
            }
            else {
                std::cout << "  WARNING: Input tree on line " << pos << " has fewer than 4 species so ignoring" << std::endl;
            }
            pos++;
        }
    }

    std::cout << "Found" << std::endl;
    std::cout << "    " << input.size() << " trees\n";
    std::cout << "    " << dict->size() << " taxa\n";

    if (mintax != maxtax && taxa_mode == "1") {
        std::cout << "    some input trees are missing taxa" << std::endl;
        std::cout << "WARNING: --shared option should NOT be used with missing taxa!" << std::endl;
    }

    fin.close();
}



void Instance::get_annotation_tree() {
    std::ifstream fin(annotation_tree_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << annotation_tree_file << std::endl;
        exit(1);
    }
    std::string line;
    std::getline(fin, line);
    fin.close();

    annotation_tree = new Tree(line, dict, indiv2taxon, support_low, support_default);
    // annotation_tree->refine();
    // annotation_tree->prepare(weight_mode, support_low, support_high, contract, support_threshold);
    std::cout << "Annotation tree/TOB has " << annotation_tree->size() << " taxa." << std::endl;
}


void Instance::input_matrix() {
    std::ifstream fin(input_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << input_file << std::endl;
        exit(1);
    }
    std::string line;
    std::getline(fin, line);
    fin.close();

    if (line[0] == '(') {
        input_trees();
        return;
    }

    CharMat *cmat = new CharMat(input_file);

    int mintax = INDEX_WIDTH;
    int maxtax = 0;
    std::size_t pos = 1;
    while (cmat->size() > 0) {
        std::string newick = cmat->pop_newick();
        if (newick != "") {
            Tree *t = new Tree(newick, dict, indiv2taxon, support_low, support_default);
            if (t->size() > maxtax) maxtax = t->size();
            if (t->size() < mintax) mintax = t->size();
            input.push_back(t);
            if (pcsonly) positions.push_back(pos);
        }
        pos++;
    }


    std::cout << "Found" << std::endl;
    std::cout << "    " << input.size() << " informative characters out of " << pos - 1 << "\n";
    std::cout << "    " << dict->size() << " taxa\n";

    if (mintax != maxtax && taxa_mode == "1") {
        std::cout << "    some input characters are missing taxa" << std::endl;
        std::cout << "WARNING: --shared option should NOT be used with missing taxa!" << std::endl;
    }
}


void Instance::refine_trees() {
    srand(refine_seed);

    std::size_t total = 0;
    for (Tree *t : input)
        total += t->refine();

    if (data_mode == "t")
        std::cout << "    " << total << " polytomies across all trees\n";

    if (total > 0 && weight_mode == "f") {
        std::cout << "  WARNING: polytomies were refined arbitrarily!" << std::endl;

        // Save refined trees for inspection!
        std::ofstream fout(input_file + ".refined");
        if (!fout.fail()) {
            for (Tree *t : input)
                fout << t->to_string() << std::endl;
            fout.close();
        }
    }
}


void Instance::prepare_trees() {
    for (Tree *t : input) {
        t->prepare(weight_mode, support_low, support_high, contract, support_threshold);
        // std::cout << t->to_string() << std::endl;
    }
}


void Instance::prepare_root_taxa() {
    std::string taxon;

    // Try to read outgroup taxon set from file (1 taxon per line)
    std::ifstream fin(root_str);
    if (!fin.fail()) {
        while (std::getline(fin, taxon)) {
            if (taxon != "")
                outgroup_taxon_set.insert(taxon);
        }
        fin.close();
        return;
    }

    // Read outgroup from comma separated string
    // TODO: use string stream with getline and delim=','
    std::vector<int> splits;

    splits.push_back(-1);
    for (int j = 0; j < root_str.size(); j++)
        if (root_str[j] == ',')
            splits.push_back(j);
    splits.push_back(root_str.size() + 1);

    for (int j = 0; j < splits.size() - 1; j++) {
        taxon = root_str.substr(splits[j]+1, splits[j+1] - splits[j] - 1);
        outgroup_taxon_set.insert(taxon);
    }
}


void Instance::prepare_indiv2taxon_map() {
    std::string taxon, indiv, line, delim;
    std::size_t sep;

    std::ifstream fin(mapping_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open mapping file " << mapping_file << std::endl;
        exit(1);
    }

    std::cout << "Reading mapping file" << std::endl;

    // Process first line and figure out the delimitator
    std::getline(fin, line);
    delim = " ";
    sep = line.find(delim);
    if (sep == std::string::npos) {
        delim = "\t";
        sep = line.find(delim);
        if (sep == std::string::npos) {
            std::cout << "\nERROR: Unable to process mapping file" << std::endl;
            exit(1);
        }
    }
    if (sep == 0 || sep + 1 == std::string::npos) {
        std::cout << "  WARNING: Unable to process mapping file line " << line << std::endl;
    }
    else {
        indiv = line.substr(0, sep);
        taxon = line.substr(sep + 1, std::string::npos);
        indiv2taxon.insert({indiv, taxon});
    }

    // Process remainder of file
    while (std::getline(fin, line)) {
        sep = line.find(delim);
        if (sep == std::string::npos || sep == 0 || sep + 1 == std::string::npos) {
            std::cout << "  WARNING: Unable to process mapping file line " << line << std::endl;
        }
        else {
            indiv = line.substr(0, sep);
            taxon = line.substr(sep + 1, std::string::npos);
            indiv2taxon.insert({indiv, taxon});
        }
    }

    //for (auto itr = indiv2taxon.begin(); itr != indiv2taxon.end(); ++itr)
    //    std::cout << itr->first << " : " << itr->second << std::endl;

    fin.close();
}


void Instance::input_quartets() {
    std::ifstream fin(input_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << input_file << std::endl;
        exit(1);
    }

    // Read the first line
    std::string firstLine;
    if (!(std::getline(fin, firstLine))) {
        std::cout << "\nERROR: Unable to read " << input_file << std::endl;
        exit(1);
    }

    if (firstLine == "t1,t2,t3,t4,CF12_34,CF13_24,CF14_23,ngenes")
        input_quartets_phylonetworks();
    else if (firstLine == "qind,t1,t2,t3,t4,CF12_34,CF13_24,CF14_23")
        input_qcfs();
    else
        input_quartets_basic();

    fin.close();
}


void Instance::input_quartets_basic() {
    std::ifstream fin(input_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << input_file << std::endl;
        exit(1);
    }

    std::vector<std::string> tokens;
    std::string delimiter = "___";
    size_t pos = 0, count = 0;
    std::string temp = quartet_format;
    while ((pos = quartet_format.find(delimiter)) != std::string::npos) {
        count ++;
        tokens.push_back(quartet_format.substr(0, pos));
        quartet_format.erase(0, pos + delimiter.length());
    }
    if (quartet_format != "") tokens.push_back(quartet_format);
    if (count != 5) {
        std::cout << "\nERROR: Invalid quartet format " << temp << std::endl;
        exit(1);
    }
    
    std::string line;
    size_t j = 0;
    while (std::getline(fin, line)) {
        j ++;

        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

        index_t indices[4];
        for (size_t i = 0; i < count; i ++) {
            size_t pos = line.find(tokens[i]);
            if (pos == std::string::npos) {
                std::cout << "\nWARNING: Invalid quartet; input truncated at line " << j << std::endl;
                j = -1;
                break;
            }
            if (i > 0)
                indices[i - 1] = dict->label2index(line.substr(0, pos));
            line.erase(0, pos + tokens[i].length());
        }
        if (j == -1) break;
        
        quartet_t quartet = join(indices);
        weight_t weight = 1.0;
        if (line != "") weight = std::stod(line);
        if (quartets.find(quartet) == quartets.end()) 
            quartets[quartet] = 0;
        quartets[quartet] += weight;
    }

    fin.close();
}


void Instance::input_quartets_phylonetworks() {
    std::ifstream fin(input_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << input_file << std::endl;
        exit(1);
    }

    std::string line;

    // Read header
    std::getline(fin, line);

    size_t j = 0;
    while (std::getline(fin, line)) {
        j ++;
        std::string temp = line;

        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

        index_t indices_12v34[4];
        index_t indices_13v24[4];
        index_t indices_14v23[4];
        weight_t weights[4];

        for (size_t i = 0; i < 7; i ++) {
            size_t pos = line.find(',');
            if (pos == std::string::npos) {
                std::cout << "\nWARNING: Invalid qCF; input truncated at line " << j << std::endl;
                j = -1;
                break;
            }

            if (i < 4)
                indices_12v34[i] = dict->label2index(line.substr(0, pos));
            else
                weights[i - 4] = std::stod(line.substr(0, pos));

            line.erase(0, pos + 1);
        }
        if (j == -1) break;

        weights[3] = std::stod(line);  // number of genes

        // Add quartet 1,2|3,4
        quartet_t quartet_12v34 = join(indices_12v34);
        quartets[quartet_12v34] = weights[0] * weights[3];

        // Add quartet 1,3|2,4
        indices_13v24[0] = indices_12v34[0];
        indices_13v24[1] = indices_12v34[2];
        indices_13v24[2] = indices_12v34[1];
        indices_13v24[3] = indices_12v34[3];
        quartet_t quartet_13v24 = join(indices_13v24);
        quartets[quartet_13v24] = weights[1] * weights[3];

        // Add 1,4|2,3
        indices_14v23[0] = indices_12v34[0];
        indices_14v23[1] = indices_12v34[3];
        indices_14v23[2] = indices_12v34[1];
        indices_14v23[3] = indices_12v34[2];
        quartet_t quartet_14v23 = join(indices_14v23);
        quartets[quartet_14v23] = weights[2] * weights[3];
    }

    fin.close();
}


void Instance::input_qcfs() {
    DataFrame df; 
    df.read_csv(input_file);
    std::cout << "Reading QCFS from " << input_file << " with " << df.data.size() << " rows." << std::endl;
    df.print_info();

    auto get_string = [](const Cell &cell) -> std::string {
        return std::get<std::string>(cell);
    };

    auto get_weight = [](const Cell &cell) -> weight_t {
        if (std::holds_alternative<float>(cell))
            return static_cast<weight_t>(std::get<float>(cell));
        if (std::holds_alternative<int>(cell))
            return static_cast<weight_t>(std::get<int>(cell));
        return std::stod(std::get<std::string>(cell)); // fallback
    };

    for (auto &row : df.data) {
        index_t indices[4];
        indices[0] = dict->label2index(get_string(row[1]));
        indices[1] = dict->label2index(get_string(row[2]));
        indices[2] = dict->label2index(get_string(row[3]));
        indices[3] = dict->label2index(get_string(row[4]));
        quartet_t quartet = join(indices);
        auto weights = std::array<weight_t, 3>{get_weight(row[5]),get_weight(row[6]),get_weight(row[7])};
        // std::cout << "Read quartet " << dict->index2label(indices[0]) << "," << dict->index2label(indices[1]) << "|"
        //           << dict->index2label(indices[2]) << "," << dict->index2label(indices[3]) 
        //           << " with qCFs: " << weights[0] << ", " << weights[1] << ", " << weights[2] << std::endl;
        auto weights_integers = std::array<size_t, 3>{
            static_cast<size_t>(weights[0] * 1e9),
            static_cast<size_t>(weights[1] * 1e9),
            static_cast<size_t>(weights[2] * 1e9)
        };
        weights[0] = weights_integers[0];
        weights[1] = weights_integers[1];
        weights[2] = weights_integers[2];
        qCFs_table[quartet] = weights;
        // std::cout << "Added quartet " << dict->index2label(indices[0]) << "," << dict->index2label(indices[1]) << "|"
        //           << dict->index2label(indices[2]) << "," << dict->index2label(indices[3]) 
        //           << " with qCFs: " << weights[0] << ", " << weights[1] << ", " << weights[2] << std::endl;
    }

}

/*
// Reads results of TINNiK tests from file 
void Instance::input_pvalues() {
    std::ifstream fin(pvalue_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open pvalue file " << pvalue_file << std::endl;
        exit(1);
    }

    std::cout << "Reading pvalue file" << std::endl;
    index_t N = dict->size();
    size_t Q = N * (N - 1) * (N - 2) * (N - 3) / 24;
    std::vector<std::string> labels;

    for (index_t i = 0; i < N; i ++)
        labels.push_back(dict->index2label(i));

    std::sort(labels.begin(), labels.end());
    std::vector<std::vector<index_t>> quartets;

    for (index_t i = 0; i < N; i ++) {
        std::vector<index_t> line;
        for (size_t j = 0; j < Q; j ++) {
            index_t k;
            fin >> k;
            line.push_back(k);
            //std::cout << k << " ";
        }
        quartets.push_back(line);
        // std::cout << i << std::endl;
    }

    std::vector<std::vector<weight_t>> pvalues;
    for (index_t i = 0; i < 5; i ++) {
        std::vector<weight_t> line;
        for (size_t j = 0; j < Q; j ++) {
            weight_t k;
            fin >> k;
            line.push_back(k);
            //std::cout << k << " ";
        }
        pvalues.push_back(line);
        // std::cout << i << std::endl;
    }

    for (size_t i = 0; i < Q; i ++) {
        index_t indices[4], k = 0;
        for (index_t j = 0; j < N; j ++) {
            if (quartets[j][i] == 1) 
                indices[k ++] = j;
        }
        for (index_t k = 0; k < 4; k ++) {
            std::string label = labels[indices[k]];
            index_t j = dict->label2index(label);
            indices[k] = j;
        }
        std::sort(indices, indices + 4);
        quartet_t q = join(indices);
        std::vector<weight_t> values;
        for (index_t j = 0; j < 5; j ++) 
            values.push_back(pvalues[j][i]);
        quartet2pvalue[q] = values;
    }

    fin.close();
}
*/
