#ifndef TREE_HPP
#define TREE_HPP

#include "utility.hpp"
#include "dict.hpp"
#include "taxa.hpp"
#include <tuple>
#include <string>
#include <vector>
#include <cstdint>
#include <random>
#include <array>

class QCFWriter;

class Node {
    friend class Network;
    friend class NetworkNode;
    friend class Tree;
    friend class SpeciesTree;
    public:
        Node(index_t index);
        Node(index_t index, bool isfake);
        ~Node();
        void new_states(index_t size);
        void delete_states();
        bool is_leaf();
        void set_parent(Node *parent);
        Node* get_parent();
        Node* get_sibling();
        void add_child(Node *child);
        bool remove_child(Node *child);
        void print_leaves_below_index();
        // std::vector<Node *> children;
        // index_t index;
    private:
        Node *parent;
        std::vector<Node *> children, ancestors;
        index_t index, size, depth;
        // std::vector<Node *> ancestors;
        // index_t size, depth;
        weight_t s1, s2, support, length;
        
        bool isfake;

        // Below is data only needed for gene trees, would be good to have GTNode class
        weight_t /* **doublet, */ *singlet;
        // std::map<index_t, weight_t> doublet;
        std::vector<std::pair<index_t, weight_t>> *doublet;
        static weight_t get_doublet(weight_t *singlet, weight_t s1, weight_t s2, index_t x, index_t y);
        weight_t get_doublet(index_t a, index_t b);
        void add_doublet(index_t a, index_t b, weight_t c);

        // for weighted quartets:
        weight_t length_, support_[2], plength, tdoublet[2], tdoublet_[2];
        weight_t *ssinglet, *ssinglet_, *pdoublet[2], *ptriplet[2], *mdoublet[2], *mdoublet_[2];
        weight_t **sdoublet[2], **sdoublet_[2], **striplet[2];

        // Below is data only needed for species trees, would be good to have STNode class
        weight_t f[3];
        #if ENABLE_TOB
        unsigned long int blob_id;
        weight_t min_f[3]; //, max_f[3];
        weight_t min_pvalue, max_pvalue;
        index_t minimizer[4]; // indices of the 4-tuple giving min p-value
        size_t split_match_count;  // number of sampled 4 taxon most frequent topology mathcing with split defined by this edge
        size_t split_mismatch_count; // number of sampled 4 taxon most frequent topology not mathcing with split defined by this edge
        std::vector<index_t> minimizers;
        std::vector<std::vector<index_t>> multi_partitions; 
        std::unordered_map<index_t, index_t> taxon2partition_id_mapping;
        index_t hybrid_index; // index of hybrid bucket in the multi_partitions
        index_t pivots[2];
        std::vector<index_t> circle_ordering; 
        #endif  // ENABLE_TOB
};

class Tree {
    friend class SpeciesTree;
    friend class Network;
    public:
        Tree();
        Tree(const std::string &newick,
             Dict *dict, 
             const std::unordered_map<std::string, std::string> &indiv2taxon, 
             weight_t support_low, weight_t support_default);
        virtual ~Tree();
        std::string to_string();
        std::string to_string_basic();
        size_t refine();
        void prepare(std::string weight_mode, weight_t low, weight_t high, bool contract, weight_t threshold);
        index_t size();
        std::unordered_map<index_t, index_t> &get_indices();
        weight_t ***build_graph(Taxa &subset);
        weight_t ***build_wgraph(Taxa &subset);
        void get_quartets(std::unordered_map<quartet_t, weight_t> *quartets);
        void get_wquartets(std::unordered_map<quartet_t, weight_t> *quartets);
        void get_wquartets_(std::unordered_map<quartet_t, weight_t> *quartets);
        std::string to_string(std::unordered_map<quartet_t, weight_t> &quartets);
        void test(Taxa &subset);
        Node* find_node(index_t index);
        Node* get_root();
        weight_t total_weight();
        void get_bipartitions(std::vector<Node *> *internal, std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> *bips);
        static std::string display_bipartition(std::vector<Node *> &A, std::vector<Node *> &B);
        void get_quardpartitions(std::vector<Node *> *internal, std::vector<std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>>> *quads, Dict *dict);
        std::string display_quardpartitions(std::vector<Node *> &A, std::vector<Node *> &B, std::vector<Node *> &C, std::vector<Node *> &D, Dict *dict);
        index_t get_quartet(index_t *indices);
    protected:
        Node *root;
        Node *pcs_node;
        std::unordered_map<index_t, Node*> index2node;
        Dict *dict;
        index_t pseudonym();
        std::string display_tree(Node *root);
        std::string display_tree_basic(Node *root);
        std::string display_tree_index(Node *root);
        Node* find_node_for_split(std::unordered_set<index_t> &clade);
        void reroot_on_edge_above_node(Node *node);
    private:
        weight_t support_low, support_default;
        weight_t total_quartet_weight;
        index_t pseudonyms;
        std::unordered_map<index_t, index_t> indices;
        std::vector<Node *> lca_euler;
        std::vector<index_t> lca_euler_depth;
        std::unordered_map<Node *, std::size_t> lca_first_occurrence;
        std::vector<std::size_t> lca_log2;
        std::vector<std::vector<std::size_t>> lca_sparse_table;
        void LCA_euler_dfs(Node *node, index_t depth);
        void clear_states(Node *root);
        void build_states(Node *root, Taxa &subset);
        void depth(Node *root, index_t depth);
        weight_t get_doublet(Node *subtree, index_t x, index_t y, bool complement);
        void sa_doublet(Node *root, weight_t sum, index_t x);
        weight_t aa_doublet(Node *root, index_t x, index_t y);
        static bool cmp(const std::pair<index_t, weight_t> &a, const std::pair<index_t, weight_t> &b);
        void sort_doublet(Node *root);
        std::unordered_set<index_t> bad_edges(Node *root, Taxa &subset, weight_t ***graph);
        void good_edges(Node *root, Taxa &subset, weight_t ***graph);
        Node *build_tree(const std::string &newick,
                         const std::unordered_map<std::string, std::string> &indiv2taxon);
        Node *build_subtree_from(Node *root);
        size_t refine_tree(Node *root);
        void prepare_tree(Node *root, std::string weight_mode, weight_t low, weight_t high, bool contract, weight_t threshold);
        //void resolve_support(Node *root);
        void add_indices(Node *root, std::vector<index_t> &indices);
        void get_leaves(Node *root, std::vector<Node *> *leaves);
        void get_bipartition(Node *root, std::vector<Node *> *A, std::vector<Node *> *B);
        void get_bipartitions(Node *root, std::vector<Node *> *internal, std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> *bips);
        void get_quardpartition(Node *root, std::vector<Node *> *A, std::vector<Node *> *B, std::vector<Node *> *C, std::vector<Node *> *D, Dict *dict);
        void get_quardpartitions(Node *root, std::vector<Node *> *internal, std::vector<std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>>> *quads, std::unordered_set<uint64_t>* seen_edges, Dict *dict);
        void get_leaf_set(Node *root, std::unordered_set<Node *> *leaf_set);
        void get_depth(Node *root, index_t depth);
        //for weighted quartets:
        void clear_wstates(Node *root);
        void build_wstates(Node *root, Taxa &subset);
        void build_ssinglet(Node *root, Taxa &subset);
        void build_ssinglet_(Node *root);
        void build_sdoublet(Node *root);
        void build_sdoublet_(Node *root);
        void build_striplet(Node *root);
        void build_striplet_(Node *root);
        std::unordered_set<index_t> wg_edges(Node *root, Taxa &subset, weight_t ***graph);
        std::unordered_set<index_t> wb_edges(Node *root, Taxa &subset, weight_t ***graph);
        template <typename function1, typename function2>
        weight_t squartet(function1 f1, function2 f2, index_t size, index_t x, index_t y);
        void test_ssinglet(Node *root, Taxa &subset);
        void test_ssinglet_(Node *root, Taxa &subset);
        void test_sdoublet(Node *root, Taxa &subset);
        void test_sdoublet_(Node *root, Taxa &subset);
        void test_striplet(Node *root, Taxa &subset);
        void test_striplet_(Node *root, Taxa &subset);
        void test_pxlet(Node *root, std::unordered_set<index_t> &subtree, Taxa &subset);
        void test_pxlet_(Node *root, std::unordered_set<index_t> &subtree, Taxa &subset);
        void test_graph(Node *root, Taxa &subset, weight_t ***graph);
        void build_wstates(Node *root);
        void build_ssinglet(Node *root, std::unordered_map<index_t, index_t> quad);
        weight_t get_qcount(std::unordered_map<index_t, index_t> quad);
        weight_t get_qfreq(std::unordered_map<index_t, index_t> quad);
        weight_t freq_(Node *root);
        void clear_wstates_(Node *root);
        void build_wstates_s(Node *root);
        void build_ssinglet_s(Node *root);
        weight_t freq_s(Node *root);
        weight_t total_weight_bf();
        void LCA_preprocessing();
        void LCA_depth_first_search(Node *root, std::vector<Node *> &stack);
        void LCA_preprocessing_with_ett_rmq_sparse_table();
        Node *LCA_via_rmq(Node *x, Node *y) const;
        Node *LCA_fast(Node *x, Node *y);
        Node *LCA_naive(Node *a, Node *b);
};

#if ENABLE_TOB
// Parameters of CornerRowContraction (Algorithm "refinement-edge-sweep"): the
// row sample size k, the heavy-sampling multiplier, the contradiction threshold
// tau in (0,1), the per-query level of the T1 test standing in for the oracle,
// and the private random seed used to draw the row samples in Phase 1.
// How the quarnet oracle is approximated on observed gene trees.
//
// The perfect oracle answers "is the quarnet on {x,y,rho,r} the tree quartet
// xy|rho r". Every stand-in answers something narrower, and which one is used
// matters more than the aggregation rule built on top. A spec is parsed from a
// string so that every method (row sweep, corner row, ...) can be pointed at the
// same set of approximations:
//
//   "t1"          T1 rejects the target at query_alpha. What has always shipped.
//   "t1+cf"       ... and the target's concordance factor is below cf_max. With
//                 1000 gene trees T1 rejects on a minor-frequency asymmetry of a
//                 few percent even when the target holds 90% of the trees; those
//                 rejections are reticulation signal from elsewhere in the 4-set,
//                 not evidence against the split. Measured: contradictions on
//                 true splits have median cf 0.83-0.89, genuine ones ~0.50.
//   "t3"          T3 rejects, i.e. the qCF vector fits no tree quartet at all.
//                 This is the closest thing MSCquartets offers to a hypothesis
//                 test whose alternative is a 4-blob, so it is the natural
//                 primitive for corner row, which is looking for nothing else.
//                 Unlike t1 it is target-free: it says "not a tree", not "not
//                 this tree", so on its own it misses a contradiction that comes
//                 from a different tree quartet.
//   "maj"         the observed majority topology is not the target, by more than
//                 `margin`. No significance test, so it is far more sensitive and
//                 far noisier; it is the only one of these that sees an ASTRAL
//                 resolution error, because those sit on branches so short that
//                 T1 has no power.
//
// Terms joined by '|' are separate *tracks*: each gets its own row threshold and
// its own corroboration requirement, and an edge is rejected if any track fires.
// That is what lets a noisy indicator be admitted under a stricter rule than a
// reliable one.
struct OracleTerm {
    bool use_t1;        // require T1 to reject the target
    bool use_t3;        // require T3 to reject every tree quartet
    bool use_cf;        // require the target concordance factor < cf_max
    bool use_majority;  // require some other topology to lead by > margin
    bool use_sym = false;  // require the two MINOR counts to be asymmetric
                        // (chi-square, 1 df).  If e is a cut edge then, once
                        // the near-side lineages fail to coalesce, they cross e
                        // as an interchangeable pair, so the two discordant
                        // topologies are exactly exchangeable whatever blobs
                        // sit elsewhere -- and asymmetry between them is the
                        // reticulation signature itself.  Measured on n15+n25
                        // this carries all of T1's discriminating power (J
                        // within 0.005 of t1 at every alpha), while T1's other
                        // half -- "the target is not the plurality" -- carries
                        // essentially none (J 0.22).  Unlike T1 it is closed
                        // form, so it needs no R call at all.
};

struct OracleSpec {
    std::vector<OracleTerm> tracks;
    weight_t query_alpha;
    weight_t cf_max;
    weight_t margin;

    // "t1" (the default) reproduces the original surrogate exactly.
    static OracleSpec parse(const std::string &text, weight_t query_alpha,
                            weight_t cf_max, weight_t margin, std::string *error);
    std::string to_string() const;
};

// Row sweep configuration. Defaults reproduce the original algorithm exactly:
// one track "t1", tau = (1+delta)/4, heavy = 1, anchors = 1.
//
// `tau` is exposed directly because deriving it as (1+delta)/4 confines delta in
// [0, 1/3) to tau in [0.25, 0.333), and almost no row lands in that band -- which
// is why the published delta sweep is nearly flat. tau = 0.45 is a small but
// consistent improvement and is unreachable through delta.
//
// `heavy[t]` demands that track t find a bad row at that many *distinct* partners
// r, capped at the number of partners the edge has (so an edge whose smaller side
// holds two taxa is unaffected). Note this is a different quantity from
// CornerRowParams::heavy, which scales a row's sample size.
//
// `anchors` > 1 repeats the sweep from several taxa of R rather than always R[0].
// It helped at n15 and not at n25; see results/analysis/ROWSWEEP_FINDINGS.md.
// How a row's far-side pair is chosen. The published rule fixes one anchor
// rho = R[0] for the whole edge, so every row (x, r) is conditioned on the pair
// (rho, r) and a single anomalous pair contaminates the entire row. See
// doc/HANDOFF.md section 5.
enum class RowSweepRowMode {
    Fixed,   // rho = R[0] for every row (published behaviour)
    Random,  // each partner column draws its own rho from R \ {r}
    Pooled   // as Random, and one row per x pooling every column, so a single
             // pair contributes 1/|R| of the statistic instead of all of it
};

struct RowSweepParams {
    OracleSpec oracle;
    std::vector<weight_t> tau;            // one per track
    std::vector<unsigned long int> heavy; // one per track
    unsigned long int anchors;            // 0 means every taxon of R
    RowSweepRowMode row_mode = RowSweepRowMode::Fixed;  // Fixed = published rule
    unsigned long int seed = 20250729;    // anchor draws for Random/Pooled

    // Second-stage threshold, one per track, applied only to edges adjacent to
    // an edge the first stage rejected.  A blob of degree k resolved binarily
    // contributes exactly k-3 spurious edges of T' and they form a *connected
    // subtree*, so every spurious edge of a degree>=5 blob has a spurious
    // neighbour by construction.  Measured over the shipped runs at n50-n200:
    // among the edges the sweep KEEPS, 84.8% of the blob-internal ones sit next
    // to a rejection against 19.9% of the true ones -- a likelihood ratio of
    // 4.25 that no per-edge test can see.  Empty, or tau2 >= tau, leaves the
    // rule exactly one-stage.
    std::vector<weight_t> tau2;

    // Optional path for a per-edge dump of the continuous statistic, so a
    // whole tau curve can be swept offline from a single run.
    std::string score_out;
};

// BranchCutContraction (theory.tex Algorithm "branch-cut-contraction").
//
// For a refinement-only edge the reticulation branch at the blob vertex is an
// edge-induced cluster of T' contained in one corner, and *every* query whose
// near-side pair crosses that cluster is a perfect contradiction. Enumerating
// the clusters therefore gives a witness family of density 1 rather than the
// density 1/2 of the corner-row lemma, which is what widens the certified noise
// range from delta < 1/3 to delta < 1/2.
//
// Two deliberate departures from the algorithm as stated, both measured on the
// recorded evidence (doc/NEW_ALGORITHMS.md section 3):
//
//   tau  The theorem rejects when a cluster's crossing coordinates contradict
//        more than half the time, because under a *perfect* oracle that family
//        is fully bad. The statistical surrogate detects a known four-cycle only
//        about 0.673 of the time, and a cluster that is only approximately the
//        reticulation branch dilutes that further, so 1/2 is far too high: it
//        halves trr_blob against the corner row. Sweeping tau down to ~0.1-0.3
//        recovers it at no cost in sensitivity. tau is therefore exposed and
//        defaults well below 1/2.
//
//   min_support  The rule maximises over ~2n clusters, so a cluster with a
//        single crossing coordinate that happens to contradict satisfies
//        C > tau*M outright. A perfect oracle's union bound absorbs this; a real
//        surrogate does not. Clusters with fewer than min_support crossing
//        coordinates are skipped.
struct BranchCutParams {
    std::vector<weight_t> tau;     // reject when C_U > tau[t] * M_U, one per
                                   // oracle track (so a noisy indicator such as
                                   // `maj` can be admitted under a stricter bar
                                   // than the reticulation certificate needs)
    unsigned long int min_support; // ignore clusters with M_U < this
    unsigned long int samples;     // per side, as a multiple of |side|; 0 = all
    // Number of edge-disjoint spanning cycles per side, i.e. the `h` of
    // Lem. cycle-cover and Thm. cycle-cover-global. 0 keeps the older uniform
    // random sampling, which is NOT the published algorithm.
    //
    // Why this matters. Tests_A(e) is the multigraph t_A K_|A| with
    // t_A = |B1||B2| parallel classes, one per anchor pair (b1,b2). Lem.
    // cycle-cover picks h edge-disjoint spanning cycles from it, and a spanning
    // cycle meets every vertex, so it crosses every nonempty proper subset an
    // even positive number of times: **M_U >= 2h for EVERY cluster U,
    // deterministically**. That bound is what the proof of
    // Thm. cycle-cover-global uses for its depth, and it is exactly what
    // uniform random sampling fails to provide -- which is why the random path
    // needs `min_support`, a guard the published algorithm does not have.
    //
    // Construction used here: take one spanning cycle per distinct anchor pair.
    // Cycles in different parallel classes are automatically edge-disjoint, so
    // no Walecki decomposition is needed whenever h <= t_A; that covers every
    // edge in this benchmark. The per-cycle vertex order is a seeded
    // permutation drawn before any query, as Ass. protocol requires.
    //
    // Budget is h*|side| coordinates per side, so `--branchcut-cycles h` is
    // budget-matched to `--branchcut-samples h`.
    unsigned long int cycles;
    // Lift the h <= t_A cap of the one-cycle-per-anchor-pair construction to
    // Lem. cycle-cover's own bound h <= t_A * floor((|A|-1)/2), taking the extra
    // cycles per parallel class as fresh random vertex orders with repeated
    // coordinates dropped. See the note in the constructor: the cap, not
    // saturation, is what made h = 64 and h = 200 indistinguishable from h = 32,
    // because for a far cherry t_A = 1.
    bool cycle_reuse;
    // Draw the extra within-class cycles as a genuine Walecki decomposition of
    // K_m instead of as fresh random vertex orders. Random orders are
    // edge-disjoint only in expectation, so `cycle_reuse` de-duplicates the
    // collisions -- and a dropped coordinate is a lost crossing, which is why
    // the shipped reuse path does NOT deliver Lem. cycle-cover's M_U >= 2h
    // (measured: 48.7% of clusters short of it at h = 32, minimum 24 of 64).
    // Walecki's construction is exactly edge-disjoint, so no coordinate is ever
    // dropped and the bound is restored by construction.
    bool walecki;
    // Draw the h anchor pairs as a uniform sample of B1 x B2 rather than as the
    // first h of the enumeration f1i = c % |B1|, f2i = (c/|B1|) % |B2|. That
    // enumeration is degenerate whenever |B1| >= h: every class then shares the
    // SAME b2 = B2[0], so all coordinates of the edge lean on one anchor taxon
    // and their errors are perfectly correlated -- which is exactly the
    // independence ACROSS parallel classes that Ass. independence assumes.
    bool anchor_spread;
    // Weight the anchor sample by a per-taxon QUALITY score instead of drawing
    // it uniformly. Anchor identity is irrelevant under the pure MSC -- quartet
    // frequencies depend only on the internal branch -- so any anchor effect
    // comes from GENE TREE ESTIMATION error: a taxon that is placed unreliably
    // in the estimated gene trees drags down the resolution margin of every
    // 4-set that contains it. The score is therefore the mean per-4-set margin
    // over a global sample, which is reference-free (no true topology needed)
    // and drawn BEFORE any edge is tested, so Ass. protocol still holds.
    // 0 disables; otherwise this is the number of global 4-sets to sample.
    unsigned long int anchor_quality;
    // Exponent on the quality score. The raw per-taxon mean margin spans only
    // about [0.8, 1.05] at n25, so gamma = 1 is nearly uniform; the exponent is
    // what decides whether the weighting bites at all.
    double anchor_power;
    // Write the per-taxon quality vector here for offline validation.
    std::string anchor_quality_out;
    // A cluster may fire only if its CONTRADICTIONS span at least this many
    // distinct anchor pairs. Ass. blocked allows arbitrary dependence WITHIN an
    // anchor pair, so the effective independent sample of a side is its number
    // of anchor pairs `a`, not M_U -- and Prop. count-floor's first-order
    // regime is then a < 1/tau, i.e. a <= 6 at tau = 0.15. Measured at n150:
    // 68% of testable error sits on edges with min(t_A,t_B) <= 4, where that
    // condition always holds. The bar is min(K, a), so a side with a single
    // anchor pair is left exactly as before; this can only REMOVE rejections,
    // which is the conjunctive-guard shape that has worked in this project.
    unsigned long int anchor_corroborate;
    // Vary the anchor pair PER COORDINATE instead of per cycle. Lem.
    // cycle-cover's crossing guarantee is a property of the NEAR-side pairs
    // only, so nothing requires a cycle's m coordinates to share one (b1,b2) --
    // and today they do, which means a thin cluster's 2h crossings carry only h
    // distinct anchors, two per anchor. Rotating gives up to 2h, doubling the
    // effective independent sample exactly where Prop. count-floor bites.
    // Needs only t_A >= 2, so unlike anchor_spread it reaches the t_min <= 4
    // population that holds 68% of the remaining error.
    bool anchor_rotate;
    // Systematic rather than uniform-random anchor assignment: walk a random
    // permutation of the anchor pool, one entry per coordinate, reshuffling on
    // each pass. Uniform sampling has multinomial variance, so a thin cluster's
    // two crossings can land on the SAME anchor by chance -- exactly what
    // rotation is meant to prevent. A systematic walk guarantees even coverage
    // and strictly lower variance at identical cost.
    bool anchor_balanced;
    // Rotate ONLY on sides with t_A <= this. Measured at n150: rotation cuts
    // false rejections where anchors are few (fn_test -0.06/rep at t_min 2-4,
    // the regime Ass. blocked calls correlated) and costs blob detection where
    // they are many (fp_blob +0.04/rep at t_min > 32, where a fixed anchor
    // concentrates the four-cycle signal). Restricting it to the first regime
    // keeps the benefit and drops the cost. 0 = no restriction; the principled
    // value is ceil(1/tau), the anchor count below which Prop. count-floor is
    // first order.
    unsigned long int anchor_rotate_max;
    // A side whose near set has only m = 2 taxa has exactly ONE near-side pair,
    // so its only variation is the anchor. The code emits ONE query per cycle
    // there (a 2-cycle would repeat the same 4-set), yet the availability
    // accounting keeps the spanning-cycle bound 2h <= t_A and so caps h at
    // t_A/2 -- HALF the available anchor pairs are never queried. Since one
    // "cycle" on two vertices contributes one crossing, the honest bound is
    // M_U >= h and h <= t_A. This flag uses it, doubling the evidence on
    // exactly the cherry sides that carry 32% of the remaining error.
    bool m2_full;
    // Minimum number of ANCHOR PAIRS (t_A = |B1||B2|) an edge must have before
    // cycle_reuse is allowed to lift its depth. 1 = lift everywhere.
    //
    // Measured at n100 (arm C vs arm A, 50 replicates, doc/FINDINGS §24): reuse
    // rescues net +11 true edges on edges with t_A in [2,32] and net **-4** on
    // edges with t_A = 1, where it is measurably harmful. The mechanism is the
    // one Ass. blocked describes: with a single anchor pair the extra cycles are
    // not new anchors, they re-query the same (b1,b2) with different near-side
    // pairs, so the depth buys CORRELATED evidence and amplifies whatever is
    // idiosyncratic about that anchor instead of averaging it away. Depth pays
    // only where there is anchor diversity to spread it over.
    unsigned long int reuse_min_anchors;
    // Minimum depth to lift a THIN edge to, using the same extra-cycle
    // mechanism as `cycle_reuse` but stopping as soon as the count floor is
    // cleared instead of going all the way to `cycles`.
    //
    // Prop. count-floor: a cluster crossed by m coordinates fires on a single
    // erroneous answer exactly when m < 1/tau, and Lem. cycle-cover gives
    // m >= 2h, so h >= 1/(2 tau) leaves the first-order regime -- at tau = 0.15
    // that is h >= 4. The obstruction is availability, not budget: taking one
    // cycle per anchor pair caps h at t_A = |B1||B2|, and a far cherry gives
    // t_A = 1. `cycle_reuse` lifts the cap to Lem. cycle-cover's own bound and
    // then draws the full `cycles`, which costs 2.4x wall clock; this lifts it
    // only to `min_depth`, which is all the proposition asks for.
    // 0 = off.
    unsigned long int min_depth;
    // Additive slack on the cluster rule: reject only when C_U > tau*M_U + trim.
    // For a BINARY indicator this is exactly the "drop the top-b contradictors"
    // trimmed statistic proposed in the earliest notes and never built.
    //
    // It is nearly free for branch cut specifically, because branch cut's
    // certificate is FULL density: Lem. reticulation-cut makes every coordinate
    // crossing the reticulation branch a contradiction, so on a collapsed edge
    // C_L ~ delta_1 * M_L = 0.673 * 2h, which at h = 32, tau = 0.15 is 43
    // against a threshold of 9.6 -- headroom of 33. Corner row's half-full row
    // (Lem. half-full-corner-row) only reaches 0.5 * delta_1 * k, headroom 12,
    // so the same trim is far riskier there.
    //
    // Meanwhile a falsely rejected true edge sits *just* above the threshold --
    // that is what makes it marginal -- so a small trim removes it. This targets
    // fn_test, the error term that grows fastest with n.
    unsigned long int trim;
    // Minimum pooled cross-corner 4-sets before the resolution margin may fire.
    // The margin is a ratio of pooled gene-tree counts; an edge whose corners are
    // tiny (a corner can be a single taxon) pools very few 4-sets and gets a
    // noisy margin, which costs sensitivity for no detection. 0 = no guard.
    unsigned long int resolution_min_pooled;
    // Edge-level test aimed at the OTHER half of N. A blob-internal edge is
    // compatible with T and needs four-cycle signal; an ASTRAL-error edge
    // *conflicts* with T, so some 4-set's perfect answer is a resolved quartet
    // with a different topology -- a first-order difference in concordance
    // factors rather than a second-order asymmetry between minors.
    //
    // For a 4-set with one taxon per corner the three quartet topologies are in
    // bijection with the three ways of pairing the corners, so summing qCFs
    // over the cross-corner family scores the three local resolutions:
    //   w0 = (A1A2)(B1B2)  the resolution T' chose,  w1, w2 = the two NNIs.
    // Reject when the pooled effect size
    //   margin = (w0 - max(w1,w2)) / (w0 + w1 + w2)
    // falls below `resolution_margin`. Measured medians: 0.26 on retained
    // edges against 0.009 on ASTRAL-error edges, and every one of the 46
    // ASTRAL-error edges observed lies in [-0.05, 0.05].
    //
    // `margin` is dimensionless. The earlier form
    //   z = (w0 - max(w1,w2)) / sqrt(w0 + max(w1,w2))
    // has identical discriminating power (held-out AUC 0.978 for both) but
    // scales as sqrt(M L) in the number of pooled 4-sets M and loci L, so its
    // threshold drifts with n: the best cut moves 4.2 -> 5.7 -> 7.3 from n15 to
    // n50 while median M grows 24 -> 135. The margin's best cut stays inside
    // [0.018, 0.044] and, more usefully, sensitivity at a *fixed* cut is stable
    // (0.988 / 0.979 / 0.985 at 0.01). Prefer `resolution_margin`;
    // `resolution_z` is kept only to reproduce earlier runs.
    //
    // A per-4-set counter is much weaker: the fraction of 4-sets whose target
    // pairing wins scores AUC 0.772 against the margin's 0.974, so pool the
    // counts rather than counting wins. Combining the margin with the
    // refinement's own support and branch length does not help either (0.979
    // against 0.978 held out), so the signal is one-dimensional.
    weight_t resolution_z;
    weight_t resolution_margin;
    // DISPERSION of the per-4-set cross-corner margins, as opposed to the
    // pooled location statistic above. Predicted by Lem. branch-restriction:
    //   retained    -- every cross-corner 4-set has w0 as its perfect answer,
    //                  so the margins are uniformly HIGH;
    //   ASTRAL error-- the split conflicts with T throughout, uniformly LOW;
    //   blob-internal - only the 4-sets meeting four distinct branches at the
    //                  blob vertex are affected, so the margins are a MIXTURE:
    //                  high variance and a low lower tail, with a mean that can
    //                  sit anywhere. Pooling averages the components together
    //                  and is blind to it, which is why the contradiction test
    //                  and the pooled margin both miss the "silent" blob edges.
    // Off by default (0). Reject when q10(mu) < margin_q10, when sd(mu) >
    // margin_sd, or when min(mu) < margin_min, whichever are set.
    weight_t margin_q10, margin_sd, margin_min;
    unsigned long int margin_min_sets;   // guard: need this many 4-sets
    // Localised margin contrast (see the note above the helper of that name):
    // the largest gap, over clusters of T' lying properly inside one corner,
    // between the mean per-4-set margin outside the cluster and inside it.
    // `mlc_d` thresholds the raw gap, `mlc_t` its Welch t. 0 = off.
    weight_t mlc_d, mlc_t;
    // Cluster-level corroboration: a cluster whose localised margin gap exceeds
    // `cluster_margin` may fire at the lowered contradiction bar `tau_low`
    // instead of `tau`. Both 0 = off, which is the published rule exactly.
    weight_t cluster_margin, tau_low, tau_high;
    unsigned long int mlc_min_group;     // 4-sets needed on each side of U
    // Seed-and-propagate (theory.tex Alg. seed-and-propagate), 0 = off.
    // Phase 1 scores an edge by the MAXIMUM over ~2n clusters, so its null is
    // inflated by the maximisation. Lem. propagation says a retained edge gives
    // a null family against ANY fixed cluster, so once the reticulation branch
    // is known no maximisation is needed and the bar can drop. Measured on
    // n15+n25: the 99th percentile of the statistic on retained edges is 0.080
    // for the max against 0.0030 for a fixed cluster, 27x tighter. Phase 2
    // re-tests the surviving edges against the clusters that fired in phase 1,
    // at this threshold. Because a fixed cluster's rate never exceeds the max,
    // it can only fire where propagate_tau < rate <= tau -- a strict addition.
    weight_t propagate_tau;
    // ---- added 2026-08-08, every one default-inert ---------------------
    //
    // CONTRAST test. The published rule compares the crossing contradiction
    // rate C_U/M_U with an absolute bar tau. That bar has to serve two jobs at
    // once: it must sit above the oracle's error rate delta_0 and below the
    // certificate's density delta_1. The trouble is that delta_0 is not a
    // constant of the run -- on a SHORT branch the local violation rate is far
    // above the global 0.0005, which is exactly where branch cut's false
    // rejections live (the h = 8 -> 32 experiment showed they are noise-limited,
    // not threshold-limited, so no choice of tau and no trim reaches them).
    //
    // But a short retained branch raises the contradiction rate UNIFORMLY over
    // the side's coordinates, while Lem. reticulation-cut says a collapsed edge
    // raises it only on the coordinates that CROSS the reticulation branch:
    // a coordinate whose near pair lies wholly inside U, or wholly outside it,
    // never meets four distinct branches at the blob vertex and so has a tree
    // quarnet. The informative quantity is therefore not the level of C_U/M_U
    // but its CONTRAST against the same edge's off-cluster rate.
    //
    // Conditioning on the side's total contradiction count removes the unknown
    // local rate entirely: under the null "contradiction is independent of
    // crossing U", C_U ~ Hypergeometric(M_all, C_all, M_U), which is Fisher's
    // exact test and needs no estimate of delta_0. Reject when that upper tail,
    // Bonferroni-adjusted over the eligible clusters (the rule maximises over
    // ~2n of them), falls below this level. 0 = off.
    //
    // It also makes the bar size-adaptive for free: at equal rate a cluster with
    // M_U = 2h is far less significant than one with M_U = h|A|/2, which the
    // fixed tau treats identically.
    weight_t contrast_alpha;
    // Per-track rule. 0 = the cluster scan above; 1 = a single GLOBAL rate test
    // on the whole side, with no maximisation over clusters.
    //
    // Why a track would want this: the cluster scan is the right shape for a
    // certificate track (t1/sym), whose signal is localised on the reticulation
    // branch. It is the wrong shape for `maj`, which fires when the target
    // quartet loses its plurality -- a property of the edge, not of any cluster.
    // Running `maj` through a max over ~2n clusters of M_U ~ 2h coordinates
    // gives it the maximisation penalty with none of the localisation benefit.
    // Empty = every track uses the cluster scan, i.e. the original behaviour.
    std::vector<int> mode;
    // One independent RNG stream per (edge, track, side) instead of a single
    // stream consumed in scan order. The early exits mean a decision taken by
    // track 0 skips track 1's draws, so under the shared stream ANY flag that
    // changes a decision also reshuffles every later sample; two configurations
    // are then compared on different query sets. With fixed streams the
    // coordinates an edge is tested on depend only on the seed, so paired
    // comparisons isolate the rule being changed. Off reproduces earlier runs.
    bool fixed_streams;
    // Draw the cycle-cover coordinates once per (edge, side) and evaluate EVERY
    // oracle track on the same 4-sets, instead of giving each track its own
    // sample. The qCF cache is keyed by 4-set, so the second track's queries
    // become cache hits and the wall clock of a two-track spec drops towards
    // that of one track. Requires fixed_streams (there is no other way to make
    // two tracks draw the same stream when the first may exit early).
    bool shared_coords;
    // Per-edge sufficient statistics for offline threshold sweeps: every
    // eligible cluster's (M_U, C_U) together with the side totals, so tau,
    // min_support, trim, the contrast level and the resolution margin can all
    // be swept from one run per replicate. Costs no queries.
    std::string score_out;
    // Every sampled cross-corner 4-set with its margin, as
    // edge_id, a1, a2, b1, b2, mu. The four corner taxa locate the 4-set inside
    // the corner structure, so ANY margin statistic -- dispersion, a
    // cluster-localised contrast, a per-taxon analysis of variance -- can be
    // built and compared offline from one run instead of one run per statistic.
    std::string quad_out;
    // SPLIT-SAMPLE CORROBORATION. Partition the h spanning cycles into g groups
    // by cycle index mod g and require the SAME cluster to fire in EVERY group
    // before the edge is rejected. 0 or 1 = off (the published rule).
    //
    // Why the same cluster and not the same edge. Lem. reticulation-cut says a
    // collapsed edge has ONE distinguished cluster -- the reticulation branch --
    // on which the contradiction density is delta_1 for every coordinate that
    // crosses it. That density does not depend on which cycles supplied the
    // coordinates, so a true certificate fires in every group. A false rejection
    // instead comes from the maximisation over ~2n clusters finding a cluster
    // whose few crossing coordinates happen to contradict; the argmax cluster is
    // then a property of the draw, and there is no reason for the SAME cluster to
    // be extreme in an independent group of cycles.
    //
    // It costs no queries: the h cycles are drawn and queried exactly as before,
    // only the counting is stratified. Since a spanning cycle crosses every
    // nonempty proper subset at least twice, every group still has
    // M_U^(j) >= 2*floor(h/g) deterministically, so Thm. cycle-cover-global
    // applies inside each group at depth h/g and completeness survives with a
    // union bound over the g groups; the false-rejection probability is a
    // product over g nearly independent draws.
    unsigned long int corroborate;
    // Fraction of the g groups that must fire, in (0, 1]. 1 = all of them (the
    // conjunction). Values below 1 give the two-level majority test that
    // Assumption anchor-blocked dependence calls for: each group is one anchor
    // stratum, a stratum votes by its own threshold rule at `tau`, and the edge
    // is rejected when more than `corroborate_frac` of the strata vote to.
    // Under blocked dependence the effective sample size is the number of
    // STRATA, not the number of coordinates, so the outer vote is the one whose
    // Chernoff bound is valid.
    weight_t corroborate_frac;
    // Inner bar a single group must clear, when corroborating. 0 = use `tau`.
    //
    // Requiring every group to clear the SAME tau as the pooled count is
    // strictly stricter than the pooled rule, so it can only cost detections.
    // The loss is concentrated on collapsed edges whose certificate is genuine
    // but thin: with M_U ~ 2h/g per group, a reticulation branch at the
    // measured delta_1 = 0.673 clears tau = 0.15 comfortably in expectation but
    // fluctuates below it in one group often enough to matter. A LOWER inner
    // bar keeps the corroboration requirement -- the same cluster must look bad
    // in both halves -- while asking less of each half, which is the right
    // shape: what distinguishes a draw-driven cluster from a certificate is
    // whether the second half sees ANY excess, not whether it independently
    // reaches the pooled threshold.
    weight_t corroborate_bar;
    // Absolute floor on a cluster's contradiction count: reject only when
    // C_U >= cmin as well as C_U > tau*M_U. 0 = off.
    //
    // The cycle cover guarantees M_U >= 2h for every cluster, so at h = 32 the
    // bar tau*M_U = 9.6 already demands ten contradictions and `cmin` is inert.
    // It binds exactly where the guarantee does not: h is capped at
    // t_A = |B1||B2|, so an edge whose far corners are a cherry has t_A = 1,
    // h collapses to 1, and M_U falls to 2 -- where C_U > 0.15*M_U is satisfied
    // by a SINGLE contradicting coordinate. Those clusters are the multiplicity
    // leak: the rule maximises over ~2n of them and the smallest ones have a
    // null tail that no fixed RATE threshold controls, only a COUNT floor does.
    unsigned long int cmin;
    weight_t query_alpha;
    unsigned long int seed;
    OracleSpec oracle;
};

struct CornerRowParams {
    unsigned long int k;        // 0 requests the largest legal value, n - 3
    unsigned long int heavy;    // row m draws min(heavy * k, |Omega_x(e)|) tuples
    std::vector<weight_t> tau;  // one per oracle track
    weight_t query_alpha;
    unsigned long int seed;
    bool cross;                 // draw the row partner y from the opposite near
                                // corner only, so every queried 4-set spans all
                                // four corners; false reproduces Omega_x(e)
    std::vector<unsigned long int> corroborate;  // bad rows demanded, per track;
                                // 1 is the original "flag on the first bad row"
    // Edge-level corner-resolution test, shared with branch cut; see the long
    // note on BranchCutParams::resolution_margin. It targets the ASTRAL-error
    // half of N, which the contradiction test cannot see because those edges
    // are not blobs at all: they *conflict* with T, so some 4-set's perfect
    // answer is a resolved quartet with a different topology. 0 disables it,
    // which reproduces the original behaviour exactly.
    weight_t resolution_margin;
    unsigned long int resolution_samples;  // multiple of the corner-size sum
    OracleSpec oracle;          // "t1" reproduces the original behaviour exactly
};
#endif  // ENABLE_TOB

class SpeciesTree : public Tree {
    public:
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, std::string mode, unsigned long int iter_limit, std::string output_file);
        SpeciesTree(std::unordered_map<quartet_t, weight_t> &input_quartets, Dict *dict, std::string mode, unsigned long int iter_limit, std::string output_file);
        SpeciesTree(std::string stree_file, Dict *dict);
        #if ENABLE_TOB
        SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, bool enable_split_test);
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, SpeciesTree* display, const RowSweepParams &row_sweep);
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, SpeciesTree* display, const CornerRowParams &corner_row);
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, SpeciesTree* display, const BranchCutParams &branch_cut);
        // Every edge-induced cluster of T' as a leaf-index bitmask (both sides
        // of every edge). The branch-cut certificate is one such cluster.
        bool quartet_counts(std::vector<Tree *> &input, index_t a, index_t b, index_t c, index_t d, weight_t out[3]);
        // Pooled cross-corner resolution score for one edge: sums raw gene-tree
        // counts over 4-sets drawn from A1 x A2 x B1 x B2 into the three
        // corner-pairing scores w0 (the pairing T' asserts), w1, w2 (the NNIs),
        // and reports the scale-free margin (w0 - max(w1,w2)) / sum(w) and the
        // older scale-dependent z. See BranchCutParams for what these measure.
        //
        // `rng` is taken by reference and consumed in place: the branch-cut
        // caller shares one stream between cluster sampling and this test, so
        // the number of draws taken here is part of that run's identity.
        // `samples` is a multiple of the corner-size SUM (0 = exhaustive),
        // which differs from the cluster scan's multiple of |side|.
        //
        // Returns false, leaving the outputs untouched, when a corner is empty
        // or no sampled 4-set resolved.
        bool corner_resolution_score(std::vector<Tree *> &input, const std::vector<index_t> *corners, unsigned long int samples, std::mt19937_64 &rng, double *margin, double *z, std::size_t *pooled, std::vector<double> *per_set = nullptr, std::vector<std::array<index_t, 4>> *per_quad = nullptr);
        void branch_cut_clusters(Node *root, std::size_t nwords, std::vector<std::vector<std::uint64_t>> &out, std::vector<std::uint64_t> &universe);
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, SpeciesTree* display, QCFWriter *qcf_writer = nullptr);
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, SpeciesTree* display, unsigned long int iter_limit_blob, QCFWriter *qcf_writer = nullptr);
        SpeciesTree(std::vector<Tree *> &input, Dict *dict, SpeciesTree* display, unsigned long int iter_limit_blob, bool three_fix_one_alter, bool two_fix_two_alter, bool is_quard, const std::string &output_qcfs_table_file = "");
        SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, std::vector<Tree *> &gene_trees, unsigned long int iter_limit_blob);
        SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit_blob);
        void hybrid_voting(std::vector<Tree *> &gene_trees,Dict *dict, Node * hybrid_blob, unsigned long int iter_limit, std::vector<std::unordered_set<index_t>> &banned_buckets);
        void hybrid_voting(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Dict *dict, Node* blob_node, unsigned long int iter_limit, std::vector<std::unordered_set<index_t>> &banned_buckets);
        void pivot_scan(std::vector<Tree *> &gene_trees, Dict *dict, Node *hybrid_blob, unsigned long int iter_limit);
        void circle_sorting(std::vector<Tree *> &gene_trees, unsigned long int iter_limit, Node* hybrid_blob_node);
        void circle_sorting_enmuerate_pivots(std::vector<Tree *> &gene_trees, unsigned long int iter_limit, Node * blob_node);
        void circle_sorting_enmuerate_pivots(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit, Node * blob_node);
        #endif  // ENABLE_TOB
        ~SpeciesTree();
        void print_leaves(std::vector<Node *> &leaves, std::ostream &os);
        void print_leaf_set(std::unordered_set<Node *> &leaf_set, std::ostream &os);
        void annotate(std::vector<Tree *> &input, std::string & qfreq_mode);
        void annotate(std::unordered_map<quartet_t, weight_t> &quartets, std::string & qfreq_mode);
        void root_at_clade(std::unordered_set<std::string> &clade_taxon_set);
        void put_back_root();
        std::string to_string_annotated(std::string brln_mode);
        void write_support_table(std::ostream &os, std::string brln_mode);
        void write_pcs_table(std::vector<Tree *> &input, std::vector<std::size_t> &positions, std::string &qfreq_mode, std::ostream &os);
        #if ENABLE_TOB
        std::string to_string_pvalue();
        std::string display_tree_pvalue(Node *root);
        void run_split_experiment(std::vector<Tree *> &input, const std::unordered_map<std::string, index_t> &name2index, const std::string &bipartition_file, const std::string &output_file, const RowSweepParams &params);
        void dump_row_sweep_evidence(std::vector<Tree *> &input, const std::string &out_prefix, const OracleSpec &oracle, bool all_targets, bool all_anchors);
        #endif  // ENABLE_TOB
    private:
        index_t artifinyms;
        std::string mode;
        unsigned long int iter_limit;
        #if ENABLE_TOB
        std::unordered_map<quartet_t,weight_t> pvalues, pvalues_star;
        std::unordered_map<quartet_t, std::array<weight_t, 3>> qCFs_cache;
        // Empirical oracle caches are isolated from the legacy T3/blob searches.
        std::unordered_map<quartet_t, std::array<weight_t, 3>> row_sweep_qcfs_cache;
        std::unordered_map<quartet_t, std::array<weight_t, 3>> row_sweep_t1_pvalues_cache;
        std::unordered_map<quartet_t, weight_t> row_sweep_t3_pvalues_cache;
        std::unordered_map<quartet_t, std::array<weight_t, 3>> qCFs_average_cache;
        #endif  // ENABLE_TOB
        index_t artifinym();
        Node *construct_stree(std::vector<Tree *> &input, Taxa &subset, index_t parent_pid, index_t depth);
        Node *construct_stree(std::unordered_map<quartet_t, weight_t> &input, Taxa &subset, index_t parent_pid, index_t depth);
        Node *reroot(Node *root, std::unordered_set<index_t> &visited);
        Node *reroot_stree(Node *root, index_t artificial);
        Node *artificial2node(Node *root, index_t artificial);
        void get_qfreq_around_branch(Node *root, std::vector<Tree *> &input, std::string &qfreq_mode);
        void get_qfreq_around_branch(Node *root, std::unordered_map<quartet_t, weight_t> &quartets, std::string &qfreq_mode);
        std::string display_tree_annotated(Node *root, std::string brln_mode);
        void write_support_table_row(Node *root, std::ostream &os, std::string brln_mode);
        #if ENABLE_TOB
        weight_t search(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, size_t iter_limit, index_t *minimizer, index_t branch_id = 0, QCFWriter *qcf_writer = nullptr);
        weight_t search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, size_t iter_limit);
        weight_t neighbor_search(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *current, weight_t *min, index_t branch_id = 0, QCFWriter *qcf_writer = nullptr);
        weight_t neighbor_search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *current, weight_t *min);
        weight_t search(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *minimizer, index_t branch_id = 0, QCFWriter *qcf_writer = nullptr);
        weight_t search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B);
        weight_t search_3f1a(std::vector<Tree *> &input, std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>> *quad, index_t* minimizer, index_t branch_id = 0, QCFWriter *qcf_writer = nullptr);
        weight_t search_quard(std::vector<Tree *> &input, std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>> *quad, index_t* minimizer, index_t branch_id = 0, QCFWriter *qcf_writer = nullptr);
        weight_t search_quard(std::vector<Tree *> &input, std::vector<std::vector<index_t>> &quad, index_t* minimizer);
        
        size_t neighbor_search_quard(std::vector<Tree *> &input,
                                         std::vector<std::vector<index_t>> &quad,
                                         index_t *current,
                                         weight_t *min,
                                         index_t branch_id = 0,
                                         QCFWriter *qcf_writer = nullptr);
        size_t neighbor_search_quard(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table,
                                         std::vector<std::vector<index_t>> &quad,
                                         index_t *current,
                                         weight_t *min);
        weight_t search_quard_heuristic(std::vector<Tree *> &input,
                                            std::vector<std::vector<index_t>> &quad,
                                            size_t iter_limit,
                                            index_t *minimizer,
                                            index_t branch_id = 0,
                                            QCFWriter *qcf_writer = nullptr);
        weight_t search_quard_heuristic(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table,
                                            std::vector<std::vector<index_t>> &quad,
                                            unsigned long int iter_limit,
                                            index_t *minimizer);
        void generate_minimizers(std::vector<Tree *> &input, Node *root, Dict *dict, unsigned long int iter_limit, weight_t alpha);
        void generate_minimizers(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Node *root, Dict *dict, unsigned long int iter_limit, weight_t alpha);
        void postorder_nodes(Node* root, std::vector<Node*>& out);
        std::array<index_t, 2> hybrid_siblings_from_top2_qcfs(const std::array<weight_t, 3>& qCFs, const std::array<index_t,4>& quad_ids, index_t hybrid_partition_id);
        weight_t F_bucket_topology(const std::vector<Tree*>& input,
                         const std::array<std::vector<index_t>,4>& buckets,
                         index_t num_taxa,
                         const std::array<int,4>& roles);
        std::array<weight_t,3> freq_three_toplogies(const std::vector<Tree*>& input,
                                       Node* blob_node,
                                       const std::array<index_t,4>& quad_ids,
                                       index_t num_taxa);
                                       

        bool is_bucket_i_less_than_bucket_j(index_t partition_i, index_t partition_j, index_t pivot, Node* blob_node, std::vector<Tree *> gene_trees, size_t iter_limit, size_t &failed_counts);
        bool is_bucket_i_less_than_bucket_j(index_t partition_i, index_t partition_j, index_t pivot_index, Node* blob_node, std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit, size_t &failed_counts);
        Node *build_refinement(Node *root, std::unordered_set<Node *> false_positive);
        weight_t get_pvalue(std::vector<Tree *> &input, index_t *indices);
        weight_t get_pvalue_star(std::vector<Tree *> &input, index_t *indices);
        std::pair<Node *, std::vector<index_t>> hybrid_info_tree(Node *root, Dict *dict, std::unordered_set<Node *> &false_positive_alpha, std::unordered_set<Node *> &false_positive_beta);
        std::vector<index_t > compute_taxon2parition_mapping(std::vector<Tree *> &input,Node *root, Dict *dict, std::vector<Node *> &hybrid_blob_nodes, std::unordered_set<Node *> & full_leaf_indices, unsigned long int iter_limit_blob, weight_t alpha);
        std::vector<index_t> compute_taxon2parition_mapping(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Node *root, Dict *dict, std::vector<Node *> &hybrid_blob_nodes, std::unordered_set<Node *> & full_leaf_nodes, unsigned long int iter_limit_blob, weight_t alpha);
        // for network circle sorting
        std::pair<weight_t, std::array<weight_t, 3>> get_pvalue_and_qCFs(std::vector<Tree *> &input, index_t *indices, index_t branch_id = 0, QCFWriter *qcf_writer = nullptr);
        std::pair<weight_t, std::array<weight_t, 3>> get_pvalue_and_qCFs(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, index_t *indices);
        std::array<std::array<index_t, 4>, 2> computed_displayed_quartet_toplogy(index_t *indices);
        std::array<std::array<index_t, 4>, 2> computed_displayed_quartet_toplogy(index_t *indices,
                                                const std::array<weight_t,3>& qcf);
        std::array<index_t, 2> siblings_in_two_best_topologies(const std::array<std::array<index_t,4>,2> &best2,
                                             index_t taxon);
        bool is_match_with_split(const std::array<weight_t,3>& qcf, index_t node_a1_id, index_t node_a2_id, index_t *indices);
        bool query_pairs_together(std::vector<Tree *> &input, index_t x, index_t y, index_t rho, index_t r, const OracleSpec &oracle, std::size_t track = 0);
        bool row_sweep_test_idx(std::vector<Tree *> &input, const std::vector<index_t> &A, const std::vector<index_t> &B, const RowSweepParams &params);
        // The continuous statistic behind row_sweep_test_idx: for each track,
        // the value the rule compares against tau[t], so that
        //   REJECT  <=>  score[t] > tau[t] for some t.
        // Only the Fixed row mode is scored; the others return an empty vector.
        std::vector<weight_t> row_sweep_scores_idx(std::vector<Tree *> &input, const std::vector<index_t> &A, const std::vector<index_t> &B, const RowSweepParams &params);
        bool corner_sets_for_edge(Tree *refinement, Node *edge, std::vector<index_t> corners[4]);
        weight_t search_2f2a(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t* minimizer, size_t &split_match_count, size_t &split_mismatch_count, index_t branch_id = 0, QCFWriter *qcf_writer = nullptr);
        #endif  // ENABLE_TOB
};


extern std::ofstream subproblem_csv, quartets_txt, good_edges_txt, bad_edges_txt;
extern std::string verbose;
extern unsigned long long count[8];

#if ENABLE_TOB
extern RInside RINS;
#endif  // ENABLE_TOB

#endif
