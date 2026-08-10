#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include "utility.hpp"
#include "tree.hpp"
#include "charmat.hpp"
#if ENABLE_TOB
#include "network.hpp"
#include "csvparser.hpp"
#endif

class Instance {
    public:
        Instance(int argc, char **argv);
        ~Instance();
        long long solve();
        SpeciesTree *get_solution();
        Network *get_network_solution();
        void output_solution();
        std::string get_execution_mode();
    private:
        std::vector<Tree *> input;
        Tree * annotation_tree;
        std::unordered_map<quartet_t, weight_t> quartets;
        std::unordered_map<quartet_t, std::array<weight_t, 3>> qCFs_table;
        std::vector<std::size_t> positions;
        Dict *dict;
        SpeciesTree *output;
        Network *output_net;
        std::unordered_set<std::string> outgroup_taxon_set;
        std::unordered_map<std::string, std::string> indiv2taxon;
        std::string root_str, quartet_format;
        std::string input_file, output_file, mapping_file, stree_file, table_file; //, pvalue_file;
        std::string output_qcfs_table_file;
        std::string rowsweep_file, rowsweep_out_file, rowsweep_dump_prefix;
        bool rowsweep_dump_all_targets, rowsweep_dump_all_anchors;
        weight_t rowsweep_delta, rowsweep_query_alpha;
        unsigned long int corner_row_k, corner_row_heavy, corner_row_seed;
        unsigned long int corner_row_resolution_samples;
        double corner_row_resolution_margin;
        double branch_cut_propagate_tau;
        unsigned long int branch_cut_cycles, branch_cut_trim, branch_cut_res_min_pooled;
        unsigned long int branch_cut_corroborate;
        weight_t branch_cut_corroborate_frac, branch_cut_corroborate_bar;
        unsigned long int branch_cut_cmin, branch_cut_min_depth;
        bool corner_row_cross;
        std::string rowsweep_row_mode_spec;
        unsigned long int rowsweep_seed;
        std::string corner_row_corroborate_spec;
        std::string oracle_spec, rowsweep_tau_spec, rowsweep_heavy_spec;
        std::string rowsweep_tau2_spec;
        std::string rowsweep_score_out;
        std::string branch_cut_tau_spec;
        weight_t branch_cut_resolution_z, branch_cut_resolution_margin;
        unsigned long int branch_cut_min_support, branch_cut_samples, branch_cut_seed;
        double branch_cut_contrast_alpha;
        double branch_cut_margin_q10, branch_cut_margin_sd, branch_cut_margin_min;
        unsigned long int branch_cut_margin_min_sets, branch_cut_mlc_min_group;
        double branch_cut_mlc_d, branch_cut_mlc_t;
        double branch_cut_cluster_margin, branch_cut_tau_low, branch_cut_tau_high;
        bool branch_cut_fixed_streams, branch_cut_cycle_reuse, branch_cut_shared_coords;
        std::string branch_cut_mode_spec, branch_cut_score_out, branch_cut_quad_out;
        unsigned long int rowsweep_anchors;
        weight_t oracle_cf_max, oracle_margin;
        std::string corner_row_tau_spec;
        std::string annotation_tree_file;
        std::string normal_mode, weight_mode, execute_mode, taxa_mode, score_mode, data_mode, brln_mode;
        unsigned long int refine_seed, cut_seed, iter_limit, iter_limit_blob;
        weight_t support_low, support_high, support_default, support_threshold, blob_threshold, alpha, beta;
        bool contract, char2tree, rootonly, pcsonly, blob, store_pvalue, load_pvalue, enable_split_test, override_file, three_fix_one_alter, two_fix_two_alter, row_sweep_blob, corner_row_blob, branch_cut_blob, quard, network;
        int parse(int argc, char **argv);
        #if ENABLE_TOB
        RowSweepParams build_row_sweep_params(std::string *error) const;
        #endif  // ENABLE_TOB
        void prepare_root_taxa();
        void prepare_indiv2taxon_map();
        void input_trees();
        void get_annotation_tree();
        void input_matrix();
        void prepare_trees();
        void refine_trees();
        void input_quartets();
        void input_quartets_basic();
        void input_quartets_phylonetworks();
        void input_qcfs();
        //void input_pvalues();
};

extern std::ofstream subproblem_csv;
extern std::string verbose;
extern unsigned long long count[8];
extern std::unordered_map<quartet_t, std::vector<weight_t>> quartet2pvalue;
#endif
