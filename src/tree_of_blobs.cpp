#if ENABLE_TOB
#include "tree.hpp"
#include "qcf_writer.hpp"
#include "rlib_dirs.hpp"
#include <random>
#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>

void get_qCFs(std::vector<Tree *> &input, index_t *indices, weight_t *qCFs) {
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);

    weight_t f[3] = {0, 0, 0};
    for (Tree *t : input) {
        index_t topology = t->get_quartet(temp);
        if (topology >= 0) f[topology] += 1;
    }
    qCFs[0] = f[0];
    qCFs[1] = f[1];
    qCFs[2] = f[2];
}




std::pair<Node *, std::vector<index_t>> SpeciesTree::hybrid_info_tree(Node *root, Dict *dict, std::unordered_set<Node *> &false_positive_alpha, std::unordered_set<Node *> &false_positive_beta) {
    std::vector<index_t> minimizers = {};
    if (root->children.size() == 0) {
        // std::cout << "Processing leaf node " << root->index << std::endl;
        Node *new_root = new Node(root->index, false);
        index2node[new_root->index] = new_root;
        return {new_root, minimizers};
    } else {
        
        // std::cout << "Processing node " << root->index << " with " << root->children.size() << " children" << std::endl;
        bool in_false_positive_alpha = false_positive_alpha.find(root) != false_positive_alpha.end();
        
        bool in_false_positive_beta = false_positive_beta.find(root) != false_positive_beta.end();
        
        Node *new_root = new Node(pseudonym(), in_false_positive_alpha || in_false_positive_beta);
        //  Node *new_root = new Node(pseudonym(), in_false_positive_alpha);
        if (in_false_positive_alpha) {
            minimizers.insert(minimizers.end(), std::begin(root->minimizer), std::end(root->minimizer));
        }

        for (Node *child : root->children) {
            // main the child and minimizers that no place to go 
            auto [new_child, child_minimizers] = hybrid_info_tree(child, dict, false_positive_alpha, false_positive_beta);
            
            
            if (!new_child->isfake) {
                new_root->children.push_back(new_child);
                // then the minimizer has a place to go 
                new_child->minimizers = std::move(child_minimizers);
                index2node[new_child->index] = new_child;
            }
            else {
                for (Node *grand_child : new_child->children) 
                    new_root->children.push_back(grand_child);
                
                minimizers.insert(minimizers.end(), child_minimizers.begin(), child_minimizers.end());
                
                new_child->children.clear();
                delete new_child;
            }
        }
        for (Node *new_child : new_root->children) 
            new_child->parent = new_root;
        return {new_root, minimizers};

    }
}

std::vector<std::array<std::vector<index_t>,4 >> sample_quads(const std::vector<std::vector<index_t>>& groups, size_t R, unsigned seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, 3);
    std::vector<std::array<std::vector<index_t>, 4>> results;
    results.reserve(2*R);

    std::vector<size_t> idx(groups.size());
    std::iota(idx.begin(), idx.end(), 0);

    for (size_t r = 0; r < 2*R; ++r) {
        std::shuffle(idx.begin(), idx.end(), rng);

        std::array<std::vector<index_t>, 4> buckets;
        for (size_t i = 0; i < idx.size(); ++i) {
            int b = static_cast<int>(i % 4);
            const auto& g = groups[idx[i]];
            buckets[b].insert(buckets[b].end(), g.begin(), g.end());
        }
        results.push_back(std::move(buckets));
    }
    return results;
}
// void sample_tree_quads(const std::vector<std::vector<index_t>>& cluster_a, 
//     const std::vector<std::vector<index_t>>& cluster_b, 
//     const std::vector<std::vector<index_t>>&cluster_c,  
//     std::vector<std::array<std::vector<index_t>,4 >>&result,
//     unsigned seed = 42) {

//     std::mt19937 rng(seed);
//     // case 1: root 
//     if (cluster_b.empty() && cluster_c.empty()) {
//         std::vector<size_t> idx(cluster_a.size());
//         std::iota(idx.begin(), idx.end(), 0);
//         std::shuffle(idx.begin(), idx.end(), rng);
//         std::array<std::vector<index_t>, 4> buckets;
//         for (size_t i = 0; i < idx.size(); ++i) {
//             int b = static_cast<int>(i % 4);
//             const auto& g = cluster_a[idx[i]];
//             buckets[b].insert(buckets[b].end(), g.begin(), g.end());
//         }
//         result.push_back(std::move(buckets));
//         std::vector<index_t> above_left = buckets[2];
//         above_left.insert(above_left.end(), buckets[3].begin(), buckets[3].end());
//         std::vector<index_t> above_right = buckets[0];
//         above_right.insert(above_right.end(), buckets[1].begin(), buckets[1].end());
//         sample_tree_quads(buckets[0], buckets[1], above_left);
//         sample_tree_quads(buckets[1], buckets[0], above_right);
//         sample_tree_quads(buckets[2], buckets[3], above_right);
//         sample_tree_quads(buckets[3], buckets[2], above_left);
//     } else if (cluster_a.size() == 1) {
//         return;
//     } else {
//         // shuffle cluster a
//         std::vector<size_t> idx(cluster_a.size());
//         std::iota(idx.begin(), idx.end(), 0);
//         std::shuffle(idx.begin(), idx.end(), rng);
//         std::array<std::vector<index_t>, 2> buckets;
//         for (size_t i = 0; i < idx.size(); ++i) {
//             int b = static_cast<int>(i % 2);
//             const auto& g = cluster_a[idx[i]];
//             buckets[b].insert(buckets[b].end(), g.begin(), g.end());
//         }
//         std::array<std::vector<index_t>, 4> new_quard_partitions;
//         new_quard_partitions[0] = buckets[0];
//         new_quard_partitions[1] = buckets[1];
        
//         for (auto &part : cluster_b) {
//             new_quard_partitions[2].insert(new_quard_partitions[2].end(), part.begin(), part.end());
//         }
        
//         for (auto &part : cluster_c) {
//             new_quard_partitions[3].insert(new_quard_partitions[3].end(), part.begin(), part.end());
//         }
//         // combine cluster b and c
//         result.push_back(std::move(new_quard_partitions));
//     }
//     std::vector<std::array<std::vector<index_t>, 4>> results;
//     results.reserve(2*R);

//     std::vector<size_t> idx(groups.size());
//     std::iota(idx.begin(), idx.end(), 0);

//     for (size_t r = 0; r < 2*R; ++r) {
//         std::shuffle(idx.begin(), idx.end(), rng);

//         std::array<std::vector<index_t>, 4> buckets;
//         for (size_t i = 0; i < idx.size(); ++i) {
//             int b = static_cast<int>(i % 4);
//             const auto& g = groups[idx[i]];
//             buckets[b].insert(buckets[b].end(), g.begin(), g.end());
//         }
//         results.push_back(std::move(buckets));
//     }
//     return results;
// }



// randomly generate minmizers for a blbo node in a tob with expected qcfs
void SpeciesTree::generate_minimizers(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Node *root, Dict *dict, unsigned long int iter_limit, weight_t alpha) {
    size_t blob_size = root->multi_partitions.size();
    std::vector<std::array<std::vector<index_t>,4 >> quads = sample_quads(root->multi_partitions, blob_size - 3, 42 + root->index);
    size_t count = 0;
    // weight_t current_pvalue = search_quard(gene_trees, quard_vec, current_minimizer);
    
    auto cmp = [](const std::pair<std::array<index_t, 4>, weight_t> &a,
                      const std::pair<std::array<index_t, 4>, weight_t> &b) {
        return a.second > b.second; // min-heap based on weight_t
    };

    std::priority_queue<std::pair<std::array<index_t, 4>, weight_t>,
                        std::vector<std::pair<std::array<index_t, 4>, weight_t>>,
                        decltype(cmp)> pq(cmp);
    for (auto vec : quads) {

        if (count >= blob_size - 3) break;

        index_t current_minimizer[4];
        std::vector<std::vector<index_t>>  quard_vec = {vec[0], vec[1], vec[2], vec[3]};

        weight_t current_pvalue = search_quard_heuristic(qCFs_table, quard_vec, iter_limit, current_minimizer);
        
        if (current_pvalue < alpha) {
            std::cout << "Generated minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(current_minimizer[i]) << " ";
            }
            std::cout << " with p-value " << current_pvalue << std::endl;

            root->minimizers.push_back(current_minimizer[0]);
            root->minimizers.push_back(current_minimizer[1]);
            root->minimizers.push_back(current_minimizer[2]);
            root->minimizers.push_back(current_minimizer[3]);
            count += 1;
        } else {
            std::cout << "Discarded minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(current_minimizer[i]) << " ";
            }
            std::cout << " with p-value " << current_pvalue << std::endl;
        }

        pq.push({{current_minimizer[0], current_minimizer[1], current_minimizer[2], current_minimizer[3]}, current_pvalue});
        

    }

    if (count < blob_size - 3) {
        std::cout << "Not enough minimizers found under alpha = " << alpha << ", filling up with best minimizers found." << std::endl;
        while (count < blob_size - 3 && !pq.empty()) {
            auto [minimizer, pvalue] = pq.top();
            pq.pop();

            root->minimizers.push_back(minimizer[0]);
            root->minimizers.push_back(minimizer[1]);
            root->minimizers.push_back(minimizer[2]);
            root->minimizers.push_back(minimizer[3]);
            count += 1;

            std::cout << "Added minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(minimizer[i]) << " ";
            }
            std::cout << " with p-value " << pvalue << std::endl;
        }
    }
   
}


// randomly generate minmizers for a blbo node in a tob 
void SpeciesTree::generate_minimizers(std::vector<Tree *> &gene_trees, Node *root, Dict *dict, unsigned long int iter_limit, weight_t alpha) {
    size_t blob_size = root->multi_partitions.size();
    std::vector<std::array<std::vector<index_t>,4 >> quads = sample_quads(root->multi_partitions, blob_size - 3, 42 + root->index);
    size_t count = 0;
    // weight_t current_pvalue = search_quard(gene_trees, quard_vec, current_minimizer);
    
    auto cmp = [](const std::pair<std::array<index_t, 4>, weight_t> &a,
                      const std::pair<std::array<index_t, 4>, weight_t> &b) {
        return a.second > b.second; // min-heap based on weight_t
    };

    std::priority_queue<std::pair<std::array<index_t, 4>, weight_t>,
                        std::vector<std::pair<std::array<index_t, 4>, weight_t>>,
                        decltype(cmp)> pq(cmp);
    for (auto vec : quads) {

        if (count >= blob_size - 3) break;

        index_t current_minimizer[4];
        std::vector<std::vector<index_t>>  quard_vec = {vec[0], vec[1], vec[2], vec[3]};
        
        auto get_nodes_from_indices_vec = [&](const std::vector<index_t> &indices) {
            std::vector<Node *> nodes;
            nodes.reserve(indices.size());
            for (index_t idx : indices) {
                nodes.push_back(index2node[idx]);
            }
            return nodes;

        }; 

        std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>> quard_nodes = {
            get_nodes_from_indices_vec(vec[0]),
            get_nodes_from_indices_vec(vec[1]),
            get_nodes_from_indices_vec(vec[2]),
            get_nodes_from_indices_vec(vec[3])
        };
             
        // weight_t current_pvalue = search_quard_heuristic(gene_trees, quard_vec, iter_limit, current_minimizer);
        weight_t current_pvalue = search_3f1a(gene_trees, &quard_nodes, current_minimizer);
        
        auto [_, qcfs] = get_pvalue_and_qCFs(gene_trees, current_minimizer);

        if (current_pvalue < alpha) {
            std::cout << "Generated minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(current_minimizer[i]) << " ";
            }
            std::cout << " with p-value " << current_pvalue << " with qcfs: ";
            for (const auto& qcf : qcfs) {
                std::cout << qcf << " ";
            }
            std::cout << std::endl;

            root->minimizers.push_back(current_minimizer[0]);
            root->minimizers.push_back(current_minimizer[1]);
            root->minimizers.push_back(current_minimizer[2]);
            root->minimizers.push_back(current_minimizer[3]);
            count += 1;
        } else {
            std::cout << "Discarded minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(current_minimizer[i]) << " ";
            }
            std::cout << " with p-value " << current_pvalue << " with qcfs: ";
            for (const auto& qcf : qcfs) {
                std::cout << qcf << " ";
            }
            std::cout << std::endl;
        }

        pq.push({{current_minimizer[0], current_minimizer[1], current_minimizer[2], current_minimizer[3]}, current_pvalue});
        

    }

    if (count < blob_size - 3) {
        std::cout << "Not enough minimizers found under alpha = " << alpha << ", filling up with best minimizers found." << std::endl;
        while (count < blob_size - 3 && !pq.empty()) {
            auto [minimizer, pvalue] = pq.top();
            pq.pop();

            root->minimizers.push_back(minimizer[0]);
            root->minimizers.push_back(minimizer[1]);
            root->minimizers.push_back(minimizer[2]);
            root->minimizers.push_back(minimizer[3]);
            count += 1;

            std::cout << "Added minimizers ";
            for (int i = 0; i < 4; i ++) {
                std::cout << dict->index2label(minimizer[i]) << " ";
            }
            std::cout << " with p-value " << pvalue << std::endl;
        }
    }
   
}

// compute taxon to partition mapping for a blob node in a tob using expected qCFS
std::vector<index_t> SpeciesTree::compute_taxon2parition_mapping(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Node *root, Dict *dict, std::vector<Node *> &hybrid_blob_nodes, std::unordered_set<Node *> & full_leaf_nodes, unsigned long int iter_limit_blob, weight_t alpha) {
    std::vector<index_t> taxon_below;
    if (root->children.size() == 0) {
        taxon_below.push_back(root->index);
    } else {
        index_t partition_id = 0;
        std::unordered_set<index_t > seen_taxa;
        for (Node *child : root->children) {
            std::vector<index_t> child_taxon_below = compute_taxon2parition_mapping(qCFs_table, child, dict, hybrid_blob_nodes, full_leaf_nodes, iter_limit_blob, alpha);
            root->multi_partitions.push_back(child_taxon_below);
            for (index_t taxon : child_taxon_below) {
                taxon_below.push_back(taxon);
                root->taxon2partition_id_mapping[taxon] = partition_id;
            }
            
            seen_taxa.insert(child_taxon_below.begin(), child_taxon_below.end());
            
            partition_id++;
        }

        
        bool flag = true;
        for (Node *taxon : full_leaf_nodes) {
            if (seen_taxa.find(taxon->index) == seen_taxa.end()) {
                if (flag) {
                    root->multi_partitions.push_back({}); // outside taxa partition
                    flag = false;
                }
                root->taxon2partition_id_mapping[taxon->index] = partition_id; // outside taxa
                root->multi_partitions.back().push_back(taxon->index);
            }
        }


        if (root->multi_partitions.size() > 4) {
            hybrid_blob_nodes.push_back(root);

            // for the case that take tob as input so we need to generate minimizers first
            size_t blob_size = root->multi_partitions.size();

            if (root->minimizers.empty()) {
                std::cout << "Generating minimizers for blob id " << root->index << " with " << blob_size << " partitions." << std::endl;
                generate_minimizers(qCFs_table, root, dict, iter_limit_blob, alpha);
            }

            // printing out the minimizers
            std::cout << "Hybrid blob id " << root->index << " with partitions: ";
            for (index_t i = 0; i < root->multi_partitions.size(); i ++) {
                std::cout << "[";
                for (index_t taxon : root->multi_partitions[i]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
                std::cout << "] ";
            }
            std::cout << " and minimizers: ";
            for (index_t taxon : root->minimizers) {
                std::cout << dict->index2label(taxon) << " ";
            }
            std::cout << std::endl;
            std::cout << "the total size of minimizers is " << root->minimizers.size() << std::endl;
        
        }
        
    }

    return taxon_below;
}


std::vector<index_t> SpeciesTree::compute_taxon2parition_mapping(std::vector<Tree *> &gene_trees, Node *root, Dict *dict, std::vector<Node *> &hybrid_blob_nodes, std::unordered_set<Node *> & full_leaf_nodes, unsigned long int iter_limit_blob, weight_t alpha) {
    std::vector<index_t> taxon_below;
    if (root->children.size() == 0) {
        taxon_below.push_back(root->index);
    } else {
        index_t partition_id = 0;
        std::unordered_set<index_t > seen_taxa;
        for (Node *child : root->children) {
            std::vector<index_t> child_taxon_below = compute_taxon2parition_mapping(gene_trees, child, dict, hybrid_blob_nodes, full_leaf_nodes, iter_limit_blob, alpha);
            root->multi_partitions.push_back(child_taxon_below);
            for (index_t taxon : child_taxon_below) {
                taxon_below.push_back(taxon);
                root->taxon2partition_id_mapping[taxon] = partition_id;
            }
            
            seen_taxa.insert(child_taxon_below.begin(), child_taxon_below.end());
            
            partition_id++;
        }

        
        bool flag = true;
        for (Node *taxon : full_leaf_nodes) {
            if (seen_taxa.find(taxon->index) == seen_taxa.end()) {
                if (flag) {
                    root->multi_partitions.push_back({}); // outside taxa partition
                    flag = false;
                }
                root->taxon2partition_id_mapping[taxon->index] = partition_id; // outside taxa
                root->multi_partitions.back().push_back(taxon->index);
            }
        }


        if (root->multi_partitions.size() > 4) {
            hybrid_blob_nodes.push_back(root);

            // for the case that take tob as input so we need to generate minimizers first
            size_t blob_size = root->multi_partitions.size();

            if (root->minimizers.empty()) {
                std::cout << "Generating minimizers for blob id " << root->index << " with " << blob_size << " partitions." << std::endl;
                generate_minimizers(gene_trees, root, dict, iter_limit_blob, alpha);
            }

            // printing out the minimizers
            std::cout << "Hybrid blob id " << root->index << " with partitions: ";
            for (index_t i = 0; i < root->multi_partitions.size(); i ++) {
                std::cout << "[";
                for (index_t taxon : root->multi_partitions[i]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
                std::cout << "] ";
            }
            std::cout << " and minimizers: ";
            for (index_t taxon : root->minimizers) {
                std::cout << dict->index2label(taxon) << " ";
            }
            std::cout << std::endl;
            std::cout << "the total size of minimizers is " << root->minimizers.size() << std::endl;
        
        }
        
    }

    return taxon_below;
}

void SpeciesTree::hybrid_voting(std::vector<Tree *> &gene_trees, Dict *dict, Node* blob_node, unsigned long int iter_limit, std::vector<std::unordered_set<index_t>> &banned_buckets) {

    // add_r_libpaths_and_load(RINS);
    // for (Tree *t : gene_trees) t->LCA_preprocessing();

    std::unordered_map<index_t, int> parititons_votes; 
    parititons_votes.reserve(blob_node->multi_partitions.size());
    std::unordered_set<index_t> banned_buckets_index_of_current_node;
    
    std::vector<std::unordered_set<index_t>> partition_sets;
    partition_sets.reserve(blob_node->multi_partitions.size());

    for (auto &p : blob_node->multi_partitions) {
        partition_sets.emplace_back(p.begin(), p.end());
    }

    for (size_t i = 0; i < banned_buckets.size(); i++) {
        for (size_t j = 0; j < blob_node->multi_partitions.size(); j++) {
            if (std::all_of( banned_buckets[i].begin(),  banned_buckets[i].end(),
            [&](index_t a) {
                return partition_sets[j].count(a) > 0;
            })) {
                banned_buckets_index_of_current_node.insert(j);
            }
        }
    }
    
    for (index_t i = 0; i < blob_node->multi_partitions.size(); i ++) {
        parititons_votes[i] = 0ULL;
        
    }

    for (size_t k = 0; k + 3 < blob_node->minimizers.size(); k += 4) {

        std::unordered_map<index_t, size_t> seen_partitions_in_quad;
        seen_partitions_in_quad.reserve(4);
            
        for (size_t t = 0; t < 4; ++t) {
            index_t minimizer = blob_node->minimizers[k + t];
            
            auto it = blob_node->taxon2partition_id_mapping.find(minimizer);
            index_t partition_id = it->second;

            if (it != blob_node->taxon2partition_id_mapping.end()) {
                seen_partitions_in_quad[partition_id] += 1;
                // seen_partitions_in_quad.insert(partition_id);
            } else {
                seen_partitions_in_quad[partition_id] = 1;
            }
        }

        for (const auto& [partition_id, votes] : seen_partitions_in_quad) {
            if (banned_buckets_index_of_current_node.find(partition_id) != banned_buckets_index_of_current_node.end()) {
                parititons_votes[partition_id] -= 1;
            } else if (votes == 1) {
                parititons_votes[partition_id] += 1;
            } else if (votes >= 2) {
                // parititons_votes[partition_id] = -100;
                parititons_votes[partition_id] += 1;
            }
            
        }
    }

    if (blob_node->multi_partitions.size() < 5) {
        std::cout << "Blob id " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization" << std::endl;
    } else {
            
            // sort the partition votes by values 

            std::vector<std::pair<index_t, int>> sorted_partitions (parititons_votes.begin(), parititons_votes.end());
            std::sort(sorted_partitions.begin(), sorted_partitions.end(), [](auto &a, auto &b) { return a.second > b.second; });
        
            // sorted_partitions.erase(std::remove_if(sorted_partitions.begin(), sorted_partitions.end(), [](std::pair<index_t, int> a) return a.second < 0), sorted_partitions.end());

            // if (sorted_partitions.size() > 5) {
            //     sorted_partitions.resize(5); // keep top 5 partitions only
            // }

            if (sorted_partitions[0].second < 0) {
                std::cout << "Blob id " << blob_node->index << " all partitions conflicted leave it as unresolved. ";
                return;
            }


            
            std::cout << "number of parititons : " << blob_node->multi_partitions.size() << std::endl;
            // print out the top 5 partitions
            std::cout << "Blob id " << blob_node->index << " partitions votes: ";
            for (auto &p : sorted_partitions) {
                std::cout << "Partition " << p.first << " with votes " << p.second << "; ";
            }
            std::cout << std::endl;

            std::cout << "sorted partitions size: " << sorted_partitions.size() << std::endl;


            // iterates all the way of choose 4 from 5 paritions

            if (sorted_partitions[0].second > sorted_partitions[1].second) {
                blob_node->hybrid_index = sorted_partitions[0].first;
            } else {

                weight_t best_pvalue = -1.0;
                index_t best_missing = -1;

                for (index_t missing = 0; missing < 5; ++missing) {
                std::vector<std::vector<index_t>> quard_vec;
                
                
                for (index_t j = 0; j < 5; j++) {
                    if (j != missing) {
                        
                        quard_vec.push_back(blob_node->multi_partitions[sorted_partitions[j].first]);
                        
                    }
                }

                // std::cout << "quard_vec size: " << quard_vec.size() << std::endl;

                
                index_t current_minimizer[4];
                // weight_t current_pvalue = search_quard(gene_trees, quard_vec, current_minimizer);
                weight_t current_pvalue = search_quard_heuristic(gene_trees, quard_vec, iter_limit, current_minimizer);
                // std::cout << "debugging line 210" << std::endl;
                auto [_, qcfs] = get_pvalue_and_qCFs(gene_trees, current_minimizer);
                std::cout << "Trying missing partition " << missing << " with minimizers ";
                for (int i = 0; i < 4; i ++) {
                    std::cout << dict->index2label(current_minimizer[i]) << " ";
                }
                std::cout << " resulting in p-value " << current_pvalue << " with qcfs: ";
                for (const auto& qcf : qcfs) {
                    std::cout << qcf << " ";
                }
                std::cout << std::endl;
                if (best_pvalue < current_pvalue) {
                    best_pvalue = current_pvalue;
                    best_missing = missing;
                }

                
            }
            // std::cout << "best missing : " << best_missing << std::endl;
            // std::cout << "Blob id " << blob_node->index << " identified hybridization partition id: " << sorted_partitions[best_missing].first << " with best p-value: " << best_pvalue << std::endl;
        
            blob_node->hybrid_index = sorted_partitions[best_missing].first;
            
        }

        std::unordered_set<index_t> banned_bucket;
        
        for (size_t j = 0; j < blob_node->multi_partitions.size(); j++) {
            if (j != blob_node->hybrid_index) {
                for (index_t taxon : blob_node->multi_partitions[j]) {
                    banned_bucket.insert(taxon);
                }
            }
        }
        
        banned_buckets.push_back(banned_bucket);
        if (parititons_votes[blob_node->hybrid_index] > 0) {
            
            std::cout << "Hybrid blob id " << blob_node->index << " identified hybridization partition : [";
            
            for (index_t taxon : blob_node->multi_partitions[blob_node->hybrid_index]) {
                std::cout << dict->index2label(taxon) << " ";
            }
            
            std::cout << "]" << std::endl;
        } else {
            blob_node->hybrid_index = blob_node->multi_partitions.size(); // leave it as unresolved
            std::cout << "Blob id " << blob_node->index << " hybridization partition : " << blob_node->hybrid_index << " conflicts with other hybrid node; leave it as unresolved." << std::endl;
        }
        
    }
}

/// hybrid voting with expected qCFS
void SpeciesTree::hybrid_voting(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, Dict *dict, Node* blob_node, unsigned long int iter_limit, std::vector<std::unordered_set<index_t>> &banned_buckets) {

    // add_r_libpaths_and_load(RINS);
    // for (Tree *t : gene_trees) t->LCA_preprocessing();

    std::unordered_map<index_t, int> parititons_votes; 
    parititons_votes.reserve(blob_node->multi_partitions.size());
    std::unordered_set<index_t> banned_buckets_index_of_current_node;
    
    std::vector<std::unordered_set<index_t>> partition_sets;
    partition_sets.reserve(blob_node->multi_partitions.size());

    for (auto &p : blob_node->multi_partitions) {
        partition_sets.emplace_back(p.begin(), p.end());
    }

    for (size_t i = 0; i < banned_buckets.size(); i++) {
        for (size_t j = 0; j < blob_node->multi_partitions.size(); j++) {
            if (std::all_of( banned_buckets[i].begin(),  banned_buckets[i].end(),
            [&](index_t a) {
                return partition_sets[j].count(a) > 0;
            })) {
                banned_buckets_index_of_current_node.insert(j);
            }
        }
    }
    
    for (index_t i = 0; i < blob_node->multi_partitions.size(); i ++) {
        parititons_votes[i] = 0ULL;
        
    }

    for (size_t k = 0; k + 3 < blob_node->minimizers.size(); k += 4) {

        std::unordered_map<index_t, size_t> seen_partitions_in_quad;
        seen_partitions_in_quad.reserve(4);
            
        for (size_t t = 0; t < 4; ++t) {
            index_t minimizer = blob_node->minimizers[k + t];
            
            auto it = blob_node->taxon2partition_id_mapping.find(minimizer);
            index_t partition_id = it->second;

            if (it != blob_node->taxon2partition_id_mapping.end()) {
                seen_partitions_in_quad[partition_id] += 1;
                // seen_partitions_in_quad.insert(partition_id);
            } else {
                seen_partitions_in_quad[partition_id] = 1;
            }
        }

        for (const auto& [partition_id, votes] : seen_partitions_in_quad) {
            if (banned_buckets_index_of_current_node.find(partition_id) != banned_buckets_index_of_current_node.end()) {
                parititons_votes[partition_id] -= 1;
            } else if (votes == 1) {
                parititons_votes[partition_id] += 1;
            } else if (votes >= 2) {
                // parititons_votes[partition_id] = -100;
                parititons_votes[partition_id] += 1;
            }
            
        }
    }

    if (blob_node->multi_partitions.size() < 5) {
        std::cout << "Blob id " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization" << std::endl;
    } else {
            
            // sort the partition votes by values 

        std::vector<std::pair<index_t, int>> sorted_partitions (parititons_votes.begin(), parititons_votes.end());
        std::sort(sorted_partitions.begin(), sorted_partitions.end(), [](auto &a, auto &b) { return a.second > b.second; });
        
        // sorted_partitions.erase(std::remove_if(sorted_partitions.begin(), sorted_partitions.end(), [](std::pair<index_t, int> a) return a.second < 0), sorted_partitions.end());

        // if (sorted_partitions.size() > 5) {
        //     sorted_partitions.resize(5); // keep top 5 partitions only
        // }

        if (sorted_partitions[0].second < 0) {
            std::cout << "Blob id " << blob_node->index << " all partitions conflicted leave it as unresolved. ";
            return;
        }


            
        std::cout << "number of parititons : " << blob_node->multi_partitions.size() << std::endl;
        // print out the top 5 partitions
        std::cout << "Blob id " << blob_node->index << " partitions votes: ";
        for (auto &p : sorted_partitions) {
            std::cout << "Partition " << p.first << " with votes " << p.second << "; ";
        }
        std::cout << std::endl;

        std::cout << "sorted partitions size: " << sorted_partitions.size() << std::endl;


            // iterates all the way of choose 4 from 5 paritions
            
        weight_t best_pvalue = -1.0;
        index_t best_missing = -1;
            
        for (index_t missing = 0; missing < 5; ++missing) {
            std::vector<std::vector<index_t>> quard_vec;
                
                
            for (index_t j = 0; j < 5; j++) {
                if (j != missing) {
                        
                    quard_vec.push_back(blob_node->multi_partitions[sorted_partitions[j].first]);
                        
                }
            }

            // std::cout << "quard_vec size: " << quard_vec.size() << std::endl;

                
            index_t current_minimizer[4];
            // weight_t current_pvalue = search_quard(gene_trees, quard_vec, current_minimizer);
            weight_t current_pvalue = search_quard_heuristic(qCFs_table, quard_vec, iter_limit, current_minimizer);
            // std::cout << "debugging line 210" << std::endl;
            if (best_pvalue < current_pvalue) {
                best_pvalue = current_pvalue;
                best_missing = missing;
            }

                
        }
        // std::cout << "best missing : " << best_missing << std::endl;
        // std::cout << "Blob id " << blob_node->index << " identified hybridization partition id: " << sorted_partitions[best_missing].first << " with best p-value: " << best_pvalue << std::endl;
        
        blob_node->hybrid_index = sorted_partitions[best_missing].first;
        std::unordered_set<index_t> banned_bucket;
        
        for (size_t j = 0; j < blob_node->multi_partitions.size(); j++) {
            if (j != blob_node->hybrid_index) {
                for (index_t taxon : blob_node->multi_partitions[j]) {
                    banned_bucket.insert(taxon);
                }
            }
        }
        
        banned_buckets.push_back(banned_bucket);
        if (parititons_votes[blob_node->hybrid_index] > 0) {
            
            std::cout << "Hybrid blob id " << blob_node->index << " identified hybridization partition : [";
            
            for (index_t taxon : blob_node->multi_partitions[blob_node->hybrid_index]) {
                std::cout << dict->index2label(taxon) << " ";
            }
            
            std::cout << "]" << std::endl;
        } else {
            blob_node->hybrid_index = blob_node->multi_partitions.size(); // leave it as unresolved
            std::cout << "Blob id " << blob_node->index << " hybridization partition : " << blob_node->hybrid_index << " conflicts with other hybrid node; leave it as unresolved." << std::endl;
        }
        
    }
}


void SpeciesTree::pivot_scan(std::vector<Tree *> &gene_trees, Dict *dict, Node *blob_node, unsigned long int iter_limit) {
    
        
    if (blob_node->multi_partitions.size() < 5) {
        std::cout << "Blob id " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization" << std::endl;
        
    } else {
        std::unordered_map<index_t, index_t> pivot_scores;
            
        index_t hybrid_partition_index = blob_node->hybrid_index;
        std::vector<index_t> multi_partitions_index_no_pivot;
            
        for (index_t i = 0; i < blob_node->multi_partitions.size(); i++) {
            pivot_scores[i] = 0;
            if (i != hybrid_partition_index) {
                multi_partitions_index_no_pivot.push_back(i);
            }
        }
        // random shuffle
        static std::mt19937 rsw_rng(12345);
        std::shuffle(multi_partitions_index_no_pivot.begin(),
                multi_partitions_index_no_pivot.end(), rsw_rng);

        // index_t partition_dmmy = multi_partitions_index_no_pivot[0];
        
        // for (index_t fixed = 1;  fixed < multi_partitions_index_no_pivot.size(); fixed ++) {
            
        //     for (index_t k = fixed + 1; k < multi_partitions_index_no_pivot.size(); k ++) {
                
        //     index_t partition_k = multi_partitions_index_no_pivot[k];
        //     index_t partititon_fixed = multi_partitions_index_no_pivot[fixed];
        //     // below is the minimum p-value method
        //     // std::vector<std::vector<index_t>> quad = {blob_node->multi_partitions[partititon_fixed], 
        //     //     blob_node->multi_partitions[partition_dmmy],
        //     //     blob_node->multi_partitions[partition_k],
        //     //     blob_node->multi_partitions[hybrid_partition_index]
        //     // };

        //     // index_t minimizer[4];
        //     // weight_t min_pvalue = search_quard_heuristic(gene_trees, quad, iter_limit, minimizer);
        //     // index_t hybrid_taxon_index_pos = -1;
        //     // for (index_t i = 0; i < 4; i++) {
        //     //     if (blob_node->taxon2partition_id_mapping[minimizer[i]] == hybrid_partition_index) {
        //     //         hybrid_taxon_index_pos = i;
        //     //         break;
        //     //     }
        //     // }
        //     // index_t hybrid_taxon_index = minimizer[hybrid_taxon_index_pos];
        //     // std::array<std::array<index_t,4>,2> displayed_quartets = computed_displayed_quartet_toplogy(minimizer);
        //     // std::array<index_t, 2> hybrid_siblings = siblings_in_two_best_topologies(displayed_quartets, hybrid_taxon_index);
        //     //above is the minimum p-value method

        //     // comment out below is average qCFs method
        //     std::array<index_t, 4> quad_ids = {
        //         multi_partitions_index_no_pivot[partititon_fixed],
        //         partition_dmmy,
        //         partition_k,
        //         hybrid_partition_index
        //     };

        //     std::sort(quad_ids.begin(), quad_ids.end());
            

        //     std::array<weight_t, 3> qCFs = freq_three_toplogies(gene_trees, blob_node, quad_ids, dict->size());
        //     auto sib_parts = hybrid_siblings_from_top2_qcfs(qCFs, quad_ids, hybrid_partition_index);
        //     std::array<index_t, 2> hybrid_siblings = {sib_parts[0], sib_parts[1]};
        //     // comment out above is average qCFs method
            


        //     if (hybrid_siblings[0] == partititon_fixed || hybrid_siblings[1] == partititon_fixed) {
        //         pivot_scores[partititon_fixed] += 1;
        //     }
        //     }
        
        // }
        // // get the top 2 highest pivot scores
        // index_t highest_pivot = -1;
        // index_t second_highest_pivot = -1;
        // index_t highest_score = -1;
        // index_t second_highest_score = -1;
        // for (const auto& [partition_id, score] : pivot_scores) {
        //     if (score > highest_score) {
        //         second_highest_score = highest_score;
        //         second_highest_pivot = highest_pivot;
        //         highest_score = score;
        //         highest_pivot = partition_id;
        //     } else if (score > second_highest_score) {
        //         second_highest_score = score;
        //         second_highest_pivot = partition_id;
        //     }
        // }
        // blob_node->pivots[0] = highest_pivot;
        // blob_node->pivots[1] = second_highest_pivot;

        // commened out the pivot scan below 
        index_t partititon_i = multi_partitions_index_no_pivot[0];
        index_t partittion_j = multi_partitions_index_no_pivot[1];
        for (index_t k = 2; k < multi_partitions_index_no_pivot.size(); k ++) {
                
            index_t partition_k = multi_partitions_index_no_pivot[k];

            std::array<index_t, 4> quad_ids = {
                partititon_i,
                partittion_j,
                partition_k,
                hybrid_partition_index
            };
                

            std::array<weight_t, 3> qCFs = freq_three_toplogies(gene_trees, blob_node, quad_ids, dict->size());
            auto sib_parts = hybrid_siblings_from_top2_qcfs(qCFs, quad_ids, hybrid_partition_index);
            std::array<index_t, 2> hybrid_siblings = {sib_parts[0], sib_parts[1]};
            partititon_i = hybrid_siblings[0];
            partittion_j = hybrid_siblings[1];

        }

        
        blob_node->pivots[0] = partititon_i < partittion_j ? partititon_i : partittion_j;
        blob_node->pivots[1] = partititon_i < partittion_j ? partittion_j : partititon_i;
        // commened out the pivot scan above

        std::cout << "blob id : " << blob_node->index << " has hybrid bucket: [ ";
        for (index_t taxon : blob_node->multi_partitions[blob_node->hybrid_index]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] and pivot bucket 0 : [ ";
            
        for (index_t taxon : blob_node->multi_partitions[blob_node->pivots[0]]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] " << std::endl;

        std::cout << "] and pivot bucket 1 : [ ";
            
        for (index_t taxon : blob_node->multi_partitions[blob_node->pivots[1]]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] " << std::endl;


        // Print all pivot scores for all buckets
        // std::cout << "pivot scores by bucket id: ";
        // for (index_t pid = 0; pid < (index_t)blob_node->multi_partitions.size(); ++pid) {
            // std::cout << pid << "=" << pivot_scores[pid];
            // if (pid + 1 < (index_t)blob_node->multi_partitions.size()) std::cout << ", ";
        // }
        std::cout << std::endl;
        std::cout << "piviot candidate 1 = " << blob_node->pivots[0] << std::endl;
        std::cout << "piviot candidate 2 = " << blob_node->pivots[1] << std::endl;

    }

}

std::array<index_t, 2> SpeciesTree::hybrid_siblings_from_top2_qcfs(const std::array<weight_t, 3>& qCFs,
                               const std::array<index_t,4>& quad_ids, index_t hybrid_partition_id) {
    int hpos = -1;
    for (int i = 0; i < 4; ++i) {
        if (quad_ids[i] == hybrid_partition_id) { hpos = i; break; }
    }
    if (hpos < 0) {
        std::cout << "Error: hybrid_partition_id not found in quad_ids." << std::endl;
        return { (index_t)-1, (index_t)-1 };
    }

    // Helper: for a given topology t in {0,1,2}, return the sibling POSITION of hpos.
    auto sibling_pos = [&](int topo) -> int {
        // topo 0: (0,1)|(2,3)
        if (topo == 0) {
            if (hpos == 0) return 1;
            if (hpos == 1) return 0;
            if (hpos == 2) return 3;
            return 2; // hpos == 3
        }
        // topo 1: (0,2)|(1,3)
        if (topo == 1) {
            if (hpos == 0) return 2;
            if (hpos == 2) return 0;
            if (hpos == 1) return 3;
            return 1; // hpos == 3
        }
        // topo 2: (0,3)|(1,2)
        // pairs: (0,3) and (1,2)
        if (hpos == 0) return 3;
        if (hpos == 3) return 0;
        if (hpos == 1) return 2;
        return 1; // hpos == 2
    };

    // pick top-2 topologies by qCF (same tie-break)
    int best = 0, second = 1;
    auto better = [&](int a, int b) -> bool {
        if (qCFs[a] > qCFs[b]) return true;
        if (qCFs[a] < qCFs[b]) return false;
        return a < b;
    };

    if (better(1, best)) { second = best; best = 1; }
    else { second = 1; }

    if (better(2, best)) {
        second = best;
        best = 2;
    } else if (better(2, second)) {
        second = 2;
    }

    // sibling PARTITION IDS (indices into multi_partitions)
    int sibpos1 = sibling_pos(best);
    int sibpos2 = sibling_pos(second);

    return { quad_ids[sibpos1], quad_ids[sibpos2] };
}


bool SpeciesTree::is_bucket_i_less_than_bucket_j(index_t partition_i, index_t partition_j, index_t pivot_index, Node* blob_node, std::vector<Tree *> gene_trees, unsigned long int iter_limit, size_t &failed_counts) {
    index_t hybrid_index = blob_node->hybrid_index;

    // iter_limit = std::ceil(2 * iter_limit / std::log(dict->size()));
    // below is minmimum p-value method
    // std::vector<std::vector<index_t>> quad = {blob_node->multi_partitions[partition_i], 
    //                 blob_node->multi_partitions[partition_j],
    //                 blob_node->multi_partitions[pivot_index],
    //                 blob_node->multi_partitions[hybrid_index]
    //             };
    
    // index_t minimizer[4];

    // weight_t min_pvalue = search_quard_heuristic(gene_trees, quad, iter_limit, minimizer);
    
    
    // index_t hybrid_taxon_index_pos = -1;
    
    // for (index_t i = 0; i < 4; i++) {
    //     if (blob_node->taxon2partition_id_mapping[minimizer[i]] == hybrid_index) {
    //         hybrid_taxon_index_pos = i;
    //         break;
    //     }
    // }
    
    // index_t hybrid_taxon_index = minimizer[hybrid_taxon_index_pos];

    // std::array<std::array<index_t,4>,2> displayed_quartets = computed_displayed_quartet_toplogy(minimizer);

    // std::array<index_t, 2> hybrid_siblings = siblings_in_two_best_topologies(displayed_quartets, hybrid_taxon_index);
    // above is minmimum p-value method

    // comment out below is average qCFs method
    std::array<index_t, 4> quad_ids = {partition_i, partition_j, pivot_index, hybrid_index};
    // std::cout << "current 4 tuple: ";
    // std::cout << " [ "
    // for (index_t i : blob_node->multi_partitions[partition_i]) {
    //         std::cout << dict->index2label(i) << " ";
    //     }
    // std::cout << " ] | [ ";
    // for (index_t i : blob_node->multi_partitions[partition_j]) {
    //         std::cout << dict->index2label(i) << " ";
    //     }
    // std::cout << " ] | [ ";
    // for (index_t i : blob_node->multi_partitions[pivot_index]) {

    //         std::cout << dict->index2label(i) << " ";
    //     }
    // std::cout << " ] | [ ";
    // for (index_t i : blob_node->multi_partitions[hybrid_index]) {
    //         std::cout << dict->index2label(i) << " ";
    //     }
    // std::cout << " ] ";

    // std::sort(quad_ids.begin(), quad_ids.end());
    std::array<weight_t, 3> qCFs = freq_three_toplogies(gene_trees, blob_node, quad_ids, dict->size());
    // std::cout << "with qCFs: ";
    // for (index_t t = 0; t < 3; t++) {
    //     std::cout << qCFs[t] << " ";
    // }
    // std::cout << std::endl;

    auto sib_parts = hybrid_siblings_from_top2_qcfs(qCFs, quad_ids, hybrid_index);
    std::array<index_t, 2> hybrid_siblings = {sib_parts[0], sib_parts[1]};
    // comment out above is average qCFs method

    bool seen_pivot_with_hybrid = false;
    bool bucket_i_less_than_j = true;
    for (index_t i = 0; i < 2; i++) {
        if (hybrid_siblings[i] == pivot_index) {
            seen_pivot_with_hybrid = true;
        }
    }

    if (!seen_pivot_with_hybrid) {
        // std::cout << "Warnning for blob id : " << blob_node->index << " pivot : " << pivot_index <<" bucket : " << partition_i << " and bucket " << partition_j << " We did not see hybrid|pivot siblings in all displayed quartet topology" << std::endl;
        failed_counts += 1;

        std::cout << "failed 4-tuple and hybrid sibling and pvalue: ";
        for (index_t i : blob_node->multi_partitions[partition_i]) {
            std::cout << dict->index2label(i) << " ";
        }
        std::cout << " | ";
        for (index_t i : blob_node->multi_partitions[partition_j]) {
            std::cout << dict->index2label(i) << " ";
        }
        std::cout << " | ";
        for (index_t i : blob_node->multi_partitions[pivot_index]) {
            std::cout << dict->index2label(i) << " ";
        }
        std::cout << " | ";
        for (index_t i : blob_node->multi_partitions[hybrid_index]) {
            std::cout << dict->index2label(i) << " ";
        }; 
        
        std::cout << " | hybrid: [";
        for (auto taxon : blob_node->multi_partitions[hybrid_index]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] ";
        std::cout << " hybrid siblings: ";
        for (index_t t = 0; t < 2; t++) {
            std::cout << "[ ";
            for (auto taxon : blob_node->multi_partitions[hybrid_siblings[t]]) {
                std::cout << dict->index2label(taxon) << " ";
            }
            std::cout << "] ";
        }
        std::cout << " | qCFs: ";
        for (index_t t = 0; t < 3; t++) {
            std::cout << qCFs[t] << " ";
        }
        std::cout << std::endl;
    }

    for (index_t i = 0; i < 2; i++) {
        if (hybrid_siblings[i] == pivot_index) {
            continue;
        } else if (hybrid_siblings[i] == partition_i) {
            bucket_i_less_than_j = false;
        }
    }
    // std::cout << "finished comparing bucket " << partition_i << " and bucket " << partition_j << " with pivot " << pivot_index << ". Result: " << (bucket_i_less_than_j ? "true" : "false") << std::endl;
    return bucket_i_less_than_j;


}

// using expected qCFs
bool SpeciesTree::is_bucket_i_less_than_bucket_j(index_t partition_i, index_t partition_j, index_t pivot_index, Node* blob_node, std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit, size_t &failed_counts) {
    index_t hybrid_index = blob_node->hybrid_index;

    // iter_limit = std::ceil(2 * iter_limit / std::log(dict->size()));
    // below is minmimum p-value method
    std::vector<std::vector<index_t>> quad = {blob_node->multi_partitions[partition_i], 
                    blob_node->multi_partitions[partition_j],
                    blob_node->multi_partitions[pivot_index],
                    blob_node->multi_partitions[hybrid_index]
                };
    

    index_t minimizer[4] = {quad[0][0], quad[1][0], quad[2][0], quad[3][0]};

    std::sort(minimizer, minimizer + 4);
    // std::cout << "quad taxon: ";
    // for (index_t t = 0; t < 4; t++) {
    //     std::cout << dict->index2label(minimizer[t]) << " ";
    // }
    auto [min_pvalue, qcf] = get_pvalue_and_qCFs(qCFs_table, minimizer);
    // std::cout << "min pvalue: " << min_pvalue << std::endl;
    // std::cout << "qCFs: ";
    // for (index_t t = 0; t < 3; t++) {
    //     std::cout << qcf[t] << " ";
    // }
    // std::cout << std::endl;

    index_t hybrid_taxon_index_pos = -1;
    
    for (index_t i = 0; i < 4; i++) {
        if (blob_node->taxon2partition_id_mapping[minimizer[i]] == hybrid_index) {
            hybrid_taxon_index_pos = i;
            break;
        }
    }
    
    index_t hybrid_taxon_index = minimizer[hybrid_taxon_index_pos];

    std::array<std::array<index_t,4>,2> displayed_quartets = computed_displayed_quartet_toplogy(minimizer, qcf);

    std::array<index_t, 2> hybrid_siblings = siblings_in_two_best_topologies(displayed_quartets, hybrid_taxon_index);
    // above is minmimum p-value method
    // std::cout << "hybrid siblings: ";
    
    // for (index_t t = 0; t < 2; t++) {
    //     std::cout << dict->index2label(hybrid_siblings[t]) << " ";
    // }
    // std::cout << std::endl;

    bool seen_pivot_with_hybrid = false;
    bool bucket_i_less_than_j = true;
    for (index_t i = 0; i < 2; i++) {
        if (blob_node->taxon2partition_id_mapping[hybrid_siblings[i]] == pivot_index) {
            seen_pivot_with_hybrid = true;
        }
    }

    if (!seen_pivot_with_hybrid) {
        // std::cout << "Warnning for blob id : " << blob_node->index << " pivot : " << pivot_index <<" bucket : " << partition_i << " and bucket " << partition_j << " We did not see hybrid|pivot siblings in all displayed quartet topology" << std::endl;
        failed_counts += 1;
        std::cout << "failed 4-tuple and hybrid sibling and pvalue: ";
        for (index_t t = 0; t < 4; t++) {
            std::cout << dict->index2label(minimizer[t]) << " ";
        }
        std::cout << " | pvalue: " << min_pvalue;
        std::cout << " | hybrid: [";
        for (auto taxon : blob_node->multi_partitions[hybrid_index]) {
            std::cout << dict->index2label(taxon) << " ";
        }
        std::cout << "] ";
        std::cout << " hybrid siblings: ";
        for (index_t t = 0; t < 2; t++) {
            std::cout << dict->index2label(hybrid_siblings[t]) << " ";
        }
        std::cout << " | qCFs: ";
        for (index_t t = 0; t < 3; t++) {
            std::cout << qcf[t] << " ";
        }
        std::cout << std::endl;

    }

    for (index_t i = 0; i < 2; i++) {
        if (blob_node->taxon2partition_id_mapping[hybrid_siblings[i]] == pivot_index) {
            continue;
        } else if (blob_node->taxon2partition_id_mapping[hybrid_siblings[i]] == partition_i) {
            bucket_i_less_than_j = false;
        }
    }
    // std::cout << "finished comparing bucket " << partition_i << " and bucket " << partition_j << " with pivot " << pivot_index << ". Result: " << (bucket_i_less_than_j ? "true" : "false") << std::endl;
    return bucket_i_less_than_j;


}


void SpeciesTree::circle_sorting_enmuerate_pivots(std::vector<Tree *> &gene_trees, unsigned long int iter_limit, Node * blob_node) {

    if (blob_node->multi_partitions.size() == 4) {
        std::cout << "blob id : " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization " << std::endl;
        return; 
    } 

    const index_t hybrid_partition_index = blob_node->hybrid_index;

    // Build base list of candidate partitions (excluding hybrid and pivot later)
    std::vector<index_t> base_ids;
    base_ids.reserve(blob_node->multi_partitions.size());
    for (index_t pid = 0; pid < blob_node->multi_partitions.size(); ++pid) {
        if (pid == hybrid_partition_index) continue;
        base_ids.push_back(pid);
    }

        

    // Track best attempt
    size_t best_failed = std::numeric_limits<size_t>::max();
    index_t best_pivot = blob_node->multi_partitions.size();
    std::vector<index_t> best_circle;

    auto try_with_pivot = [&](index_t pivot_partition_index) -> void {

        // ids excluding pivot & hybrid
        std::vector<index_t> ids;
        ids.reserve(blob_node->multi_partitions.size());
        for (index_t pid : base_ids) {
            if (pid == pivot_partition_index) continue;
            ids.push_back(pid);
        }

        size_t failed_counts = 0;

        std::sort(ids.begin(), ids.end(),
                    [&](index_t a, index_t b) {
                        return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
                                                            gene_trees, iter_limit,
                                                            failed_counts);
                    });

        // auto comp = [&](index_t a, index_t b) {
        //                 return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
        //                                                     gene_trees, iter_limit,
        //                                                     failed_counts);
        //             };


        // riffleSort(ids, 0.0, comp);

        std::vector<index_t> circle_ordering;
        circle_ordering.reserve(blob_node->multi_partitions.size());
        circle_ordering.push_back(hybrid_partition_index);
        circle_ordering.push_back(pivot_partition_index);
        circle_ordering.insert(circle_ordering.end(), ids.begin(), ids.end());

        // Keep best (min failed). Tie-break: keep the first found.
        if (best_pivot == blob_node->multi_partitions.size() || failed_counts < best_failed) {
            best_failed = failed_counts;
            best_pivot = pivot_partition_index;
            best_circle = std::move(circle_ordering);
        }
    };

    // 1) Try the first pivot

    for (auto pivot_partition_index : base_ids) {
        std::cout << "Trying pivot partition id : " << pivot_partition_index << std::endl;
        try_with_pivot(pivot_partition_index);
        std::cout << "Finished trying pivot partition id : " << pivot_partition_index << " with failed counts : " << best_failed << std::endl;
        if (best_failed == 0) break; // early stop if perfect

    }
    

        // Apply best result
    if (best_pivot != blob_node->multi_partitions.size()) {
        blob_node->circle_ordering = best_circle;

        std::cout << "blob id " << blob_node->index
                    << " chosen pivot=" << best_pivot
                    << " failed_count=" << best_failed
                    << " circle ordering:\n";

        for (index_t pid : blob_node->circle_ordering) {
            std::cout << "  bucket " << pid << ": [ ";
            if (pid >= 0 && pid < (index_t)blob_node->multi_partitions.size()) {
                for (index_t taxon : blob_node->multi_partitions[pid]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
            } else {
                std::cout << "(invalid bucket id) ";
            }
            std::cout << "]\n";
        }
    } else {
        std::cout << "Warning: blob id " << blob_node->index
                    << " no valid pivot found for circle sorting.\n";
    }
}

// circle sorting with enumerating pivots using expected qCFS
void SpeciesTree::circle_sorting_enmuerate_pivots(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit, Node * blob_node) {

    if (blob_node->multi_partitions.size() == 4) {
        std::cout << "blob id : " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization " << std::endl;
        return; 
    } 

    const index_t hybrid_partition_index = blob_node->hybrid_index;

    // Build base list of candidate partitions (excluding hybrid and pivot later)
    std::vector<index_t> base_ids;
    base_ids.reserve(blob_node->multi_partitions.size());
    for (index_t pid = 0; pid < blob_node->multi_partitions.size(); ++pid) {
        if (pid == hybrid_partition_index) continue;
        base_ids.push_back(pid);
    }

        

    // Track best attempt
    size_t best_failed = std::numeric_limits<size_t>::max();
    index_t best_pivot = blob_node->multi_partitions.size();
    std::vector<index_t> best_circle;

    auto try_with_pivot = [&](index_t pivot_partition_index) -> void {

        // ids excluding pivot & hybrid
        std::vector<index_t> ids;
        ids.reserve(blob_node->multi_partitions.size());
        for (index_t pid : base_ids) {
            if (pid == pivot_partition_index) continue;
            ids.push_back(pid);
        }

        size_t failed_counts = 0;

        std::sort(ids.begin(), ids.end(),
                    [&](index_t a, index_t b) {
                        return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
                                                            qCFs_table, iter_limit,
                                                            failed_counts);
                    });

        std::vector<index_t> circle_ordering;
        circle_ordering.reserve(blob_node->multi_partitions.size());
        circle_ordering.push_back(hybrid_partition_index);
        circle_ordering.push_back(pivot_partition_index);
        circle_ordering.insert(circle_ordering.end(), ids.begin(), ids.end());

        // Keep best (min failed). Tie-break: keep the first found.
        if (best_pivot == blob_node->multi_partitions.size() || failed_counts <= best_failed) {
            best_failed = failed_counts;
            best_pivot = pivot_partition_index;
            best_circle = std::move(circle_ordering);
        }
    };

    // 1) Try the first pivot

    for (auto pivot_partition_index : base_ids) {
        std::cout << "Trying pivot partition id : " << pivot_partition_index << std::endl;
        try_with_pivot(pivot_partition_index);
        std::cout << "Finished trying pivot partition id : " << pivot_partition_index << " with failed counts : " << best_failed << std::endl;
        if (best_failed == 0) break; // early stop if perfect

    }
    

        // Apply best result
    if (best_pivot != blob_node->multi_partitions.size()) {
        blob_node->circle_ordering = best_circle;

        std::cout << "blob id " << blob_node->index
                    << " chosen pivot=" << best_pivot
                    << " failed_count=" << best_failed
                    << " circle ordering:\n";

        for (index_t pid : blob_node->circle_ordering) {
            std::cout << "  bucket " << pid << ": [ ";
            if (pid >= 0 && pid < (index_t)blob_node->multi_partitions.size()) {
                for (index_t taxon : blob_node->multi_partitions[pid]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
            } else {
                std::cout << "(invalid bucket id) ";
            }
            std::cout << "]\n";
        }
    } else {
        std::cout << "Warning: blob id " << blob_node->index
                    << " no valid pivot found for circle sorting.\n";
    }
}

void SpeciesTree::circle_sorting(std::vector<Tree *> &gene_trees, unsigned long int iter_limit, Node * blob_node) {

    if (blob_node->multi_partitions.size() == 4) {
        std::cout << "blob id : " << blob_node->index << " has the degree of less than 5, unable to identify the hybridization " << std::endl;
        // std::cout << "Using default circle ordering for 4 partitions. " << std::endl;
        // blob_node->circle_ordering.push_back(0);
        // blob_node->circle_ordering.push_back(1);
        // blob_node->circle_ordering.push_back(3);
        // blob_node->circle_ordering.push_back(2);
            
        return; 
    } 

    const index_t hybrid_partition_index = blob_node->hybrid_index;

    // Build base list of candidate partitions (excluding hybrid and pivot later)
    std::vector<index_t> base_ids;
    base_ids.reserve(blob_node->multi_partitions.size());
    for (index_t pid = 0; pid < blob_node->multi_partitions.size(); ++pid) {
        if (pid == hybrid_partition_index) continue;
        base_ids.push_back(pid);
    }

        

    // Track best attempt
    size_t best_failed = std::numeric_limits<size_t>::max();
    index_t best_pivot = -1;
    std::vector<index_t> best_circle;

    auto try_with_pivot = [&](index_t pivot_partition_index) -> void {

        // ids excluding pivot & hybrid
        std::vector<index_t> ids;
        ids.reserve(blob_node->multi_partitions.size());
        for (index_t pid : base_ids) {
            if (pid == pivot_partition_index) continue;
            ids.push_back(pid);
        }

        size_t failed_counts = 0;

        std::sort(ids.begin(), ids.end(),
                    [&](index_t a, index_t b) {
                        return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
                                                            gene_trees, iter_limit,
                                                            failed_counts);
                    });

        // auto comp = [&](index_t a, index_t b) {
        //                 return is_bucket_i_less_than_bucket_j(a, b, pivot_partition_index, blob_node, 
        //                                                     gene_trees, iter_limit,
        //                                                     failed_counts);
        //             };


        // riffleSort(ids, 0.0, comp);

        std::vector<index_t> circle_ordering;
        circle_ordering.reserve(blob_node->multi_partitions.size());
        circle_ordering.push_back(hybrid_partition_index);
        circle_ordering.push_back(pivot_partition_index);
        circle_ordering.insert(circle_ordering.end(), ids.begin(), ids.end());

        // Keep best (min failed). Tie-break: keep the first found.
        if (failed_counts < best_failed) {
            best_failed = failed_counts;
            best_pivot = pivot_partition_index;
            best_circle = std::move(circle_ordering);
        }
    };

    // 1) Try the first pivot
    try_with_pivot(blob_node->pivots[0]);

    // 2) If not perfect, try other pivots and pick min failed
    if (best_failed > 0) {
        for (size_t pi = 1; pi < 2; ++pi) {
            try_with_pivot(blob_node->pivots[pi]);
            if (best_failed == 0) break; // early stop if perfect
        }
    }

        // Apply best result
    if (best_pivot >= 0) {
        blob_node->circle_ordering = best_circle;

        std::cout << "blob id " << blob_node->index
                    << " chosen pivot=" << best_pivot
                    << " failed_count=" << best_failed
                    << " circle ordering:\n";

        for (index_t pid : blob_node->circle_ordering) {
            std::cout << "  bucket " << pid << ": [ ";
            if (pid >= 0 && pid < (index_t)blob_node->multi_partitions.size()) {
                for (index_t taxon : blob_node->multi_partitions[pid]) {
                    std::cout << dict->index2label(taxon) << " ";
                }
            } else {
                std::cout << "(invalid bucket id) ";
            }
            std::cout << "]\n";
        }
    } else {
        std::cout << "Warning: blob id " << blob_node->index
                    << " no valid pivot found for circle sorting.\n";
    }
}
  

// void SpeciesTree::circle_sorting(std::vector<Tree *> gene_trees, size_t iter_limit, std::vector<Node *> &hybrid_blob_nodes) {

//     for (Node *blob_node : hybrid_blob_nodes) {
//         if (blob_node->multi_partitions.size() > 4) {
//             std::vector<index_t> multi_partitions_index_no_pivot_hybrid;
//             index_t hybrid_partition_index = blob_node->hybrid_index;
//             index_t pivot_partition_index = blob_node->pivots[0];
//             std::vector<index_t> multi_partitions_index_no_pivot;
//             std::vector<index_t> circle_ordering;
//             circle_ordering.push_back(hybrid_partition_index);
//             circle_ordering.push_back(pivot_partition_index);
//             for (index_t i = 0; i < blob_node->multi_partitions.size(); i++) {
                
//                 if (i != hybrid_partition_index && i != pivot_partition_index) {
//                     multi_partitions_index_no_pivot.push_back(i);
//                 }
//             }

//             size_t failed_counts = 0;
//             std::sort(multi_partitions_index_no_pivot.begin(), multi_partitions_index_no_pivot.end(), 
//                 [&](index_t a, index_t b) {
//                 return is_bucket_i_less_than_bucket_j(a, b, blob_node, gene_trees, iter_limit, failed_counts);
//             });

//             circle_ordering.insert(circle_ordering.end(), multi_partitions_index_no_pivot.begin(), multi_partitions_index_no_pivot.end());

//             blob_node->circle_ordering = circle_ordering;

//             if (failed_counts > 0) {
//                 size_t failed_counts_another = 0;

//             }
            
//             std::cout << "blob id " << blob_node->index << " circle ordering:\n";
//             for (index_t pid : blob_node->circle_ordering) {
//                 std::cout << "  bucket " << pid << ": [ ";
//                 if (pid >= 0 && pid < (index_t)blob_node->multi_partitions.size()) {
//                     for (index_t taxon : blob_node->multi_partitions[pid]) {
//                         std::cout << dict->index2label(taxon) << " ";
//                     }
//                 } else {
//                     std::cout << "(invalid bucket id) ";
//                 }
//                 std::cout << "]\n";
//             }
//         }
        
//     }

// }


// compute network with expected qcfs table
SpeciesTree::SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, unsigned long int iter_limit_blob) {
        
        std::cout << "Contracting branches with alpha = " << alpha << " and beta = " << beta << std::endl;
        
        this->dict = dict;
        
        std::vector<Node *> internal;
        
        std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
        
        input->get_bipartitions(&internal, &bips);
        
        std::cout << bips.size() << " branches to test" << std::endl;
        
        std::unordered_set<Node *> false_positive_beta;
        
        std::unordered_set<Node *> false_positive_alpha;

        // std::unordered_set<Node *> false_positive;

        for (index_t i = 0; i < internal.size(); i ++) {
                
            // if (internal[i]->min_pvalue < alpha) {
            //     false_positive_alpha.insert(internal[i]);
            //     // false_positive.insert(internal[i]);
            // } else if (internal[i]->max_pvalue > beta) {
            //     false_positive_beta.insert(internal[i]);
            //     // false_positive.insert(internal[i]);
            // }
                
        }

        std::cout << "number of branches in the input tree: " << internal.size() << std::endl;

        auto [new_root, root_min] = hybrid_info_tree(input->root, dict, false_positive_alpha, false_positive_beta);
        root = new_root;
        root->minimizers = root_min;
        index2node[root->index] = root;
        std::cout << "Number of branches to contract under alpha criterion: " << false_positive_alpha.size() << std::endl;
        std::cout << "Number of branches to contract under beta criterion: " << false_positive_beta.size() << std::endl;
        

        std::cout << "The tree of blobs : " << input->display_tree_basic(root) << std::endl;
        std::unordered_set<Node *> full_leaf_nodes;
        get_leaf_set(root, &full_leaf_nodes);
        std::unordered_set<Node *> full_leaf_indices;
        for (Node *leaf : full_leaf_nodes) {
            full_leaf_indices.insert(leaf);
        }


        add_r_libpaths_and_load(RINS);

        std::vector<Node *> hybrid_blob_nodes;

        compute_taxon2parition_mapping(qCFs_table, root, dict, hybrid_blob_nodes, full_leaf_indices, iter_limit_blob, alpha);

        // for debuging hybrid_infor_tree
        
        // root = build_refinement(input->root, false_positive);
        // std::cout << "TOB using build_refinement: " << input->display_tree_basic(root) << std::endl;

        

        std::sort(hybrid_blob_nodes.begin(), hybrid_blob_nodes.end(), [](const Node * a, Node * b){
            return a->multi_partitions.size() > b->multi_partitions.size();
        });
        std::vector<std::unordered_set<index_t>> banned_buckets;

        

        std::cout << "Printing output number of nontrival blobs:" << hybrid_blob_nodes.size() << std::endl;
        
        for (Node * blob_node : hybrid_blob_nodes) {
            qCFs_average_cache.clear();
            hybrid_voting(qCFs_table, dict, blob_node, iter_limit_blob, banned_buckets);
            std::cout << "Done with hybrid voting for blob id : " << blob_node->index << std::endl;
            
            
            if (blob_node->hybrid_index == blob_node->multi_partitions.size()) {
                std::cout << "blob id : " << blob_node->index << " has no hybridization detected, skip circle sorting " << std::endl;
                continue;
            }

            // pivot_scan(gene_trees, dict, blob_node, iter_limit_blob);
            // circle_sorting(gene_trees, iter_limit_blob, blob_node);

            circle_sorting_enmuerate_pivots(qCFs_table, iter_limit_blob, blob_node);
            std::cout << "Done with circle sorting for blob id : " << blob_node->index << std::endl;

        }
        
        
    }

// given tob-annotatiation tree and compute the m-partion maps and its inverse for all  nontrival blobs and the hybrization parititons
SpeciesTree::SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, std::vector<Tree *> &gene_trees, unsigned long int iter_limit_blob) {
        
        std::cout << "Contracting branches with alpha = " << alpha << " and beta = " << beta << std::endl;
        
        this->dict = dict;
        
        std::vector<Node *> internal;
        
        std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
        
        input->get_bipartitions(&internal, &bips);
        
        std::cout << bips.size() << " branches to test" << std::endl;
        
        std::unordered_set<Node *> false_positive_beta;
        
        std::unordered_set<Node *> false_positive_alpha;

        // std::unordered_set<Node *> false_positive;

        for (index_t i = 0; i < internal.size(); i ++) {
                
            if (internal[i]->min_pvalue < alpha) {
                false_positive_alpha.insert(internal[i]);
                // false_positive.insert(internal[i]);
            } 
            // else if (internal[i]->max_pvalue > beta) {
                // false_positive_beta.insert(internal[i]);
                // false_positive.insert(internal[i]);
            // }
                
        }

        // std::unordered_set<Node*> internal_set(internal.begin(), internal.end());

        // std::function<void(Node*)> show_missing = [&](Node* u){
        //     if (!u) return;
        //     if (!u->children.empty() && !internal_set.count(u)) {
        //     std::cerr << "[MISSING] ptr=" << u
        //             << " deg=" << u->children.size()
        //             << " min_p=" << u->min_pvalue
        //             << " max_p=" << u->max_pvalue
        //             << "\n";
        //         }
        //     for (auto* c : u->children) show_missing(c);
        //     };
        // show_missing(input->root);

        std::cout << "number of branches in the input tree: " << internal.size() << std::endl;

        auto [new_root, root_min] = hybrid_info_tree(input->root, dict, false_positive_alpha, false_positive_beta);
        root = new_root;
        root->minimizers = root_min;
        index2node[root->index] = root;
        std::cout << "Number of branches to contract under alpha criterion: " << false_positive_alpha.size() << std::endl;
        std::cout << "Number of branches to contract under beta criterion: " << false_positive_beta.size() << std::endl;
        

        std::cout << "The tree of blobs : " << input->display_tree_basic(root) << std::endl;
        std::unordered_set<Node *> full_leaf_nodes;
        get_leaf_set(root, &full_leaf_nodes);
        std::unordered_set<Node *> full_leaf_indices;
        for (Node *leaf : full_leaf_nodes) {
            full_leaf_indices.insert(leaf);
        }


        add_r_libpaths_and_load(RINS);
        for (Tree *t : gene_trees) t->LCA_preprocessing();

        std::vector<Node *> hybrid_blob_nodes;

        compute_taxon2parition_mapping(gene_trees, root, dict, hybrid_blob_nodes, full_leaf_indices, iter_limit_blob, alpha);

        // for debuging hybrid_infor_tree
        
        // root = build_refinement(input->root, false_positive);
        // std::cout << "TOB using build_refinement: " << input->display_tree_basic(root) << std::endl;

        

        std::sort(hybrid_blob_nodes.begin(), hybrid_blob_nodes.end(), [](const Node * a, Node * b){
            return a->multi_partitions.size() > b->multi_partitions.size();
        });
        std::vector<std::unordered_set<index_t>> banned_buckets;

        

        std::cout << "Printing output number of nontrival blobs:" << hybrid_blob_nodes.size() << std::endl;
        
        for (Node * blob_node : hybrid_blob_nodes) {
            qCFs_average_cache.clear();
            hybrid_voting(gene_trees, dict, blob_node, iter_limit_blob, banned_buckets);
            std::cout << "Done with hybrid voting for blob id : " << blob_node->index << std::endl;
            
            
            if (blob_node->hybrid_index == blob_node->multi_partitions.size()) {
                std::cout << "blob id : " << blob_node->index << " has no hybridization detected, skip circle sorting " << std::endl;
                continue;
            }

            // pivot_scan(gene_trees, dict, blob_node, iter_limit_blob);
            // circle_sorting(gene_trees, iter_limit_blob, blob_node);

            circle_sorting_enmuerate_pivots(gene_trees, iter_limit_blob, blob_node);
            std::cout << "Done with circle sorting for blob id : " << blob_node->index << std::endl;

        }
        
        
    }


// load p-value and just testing it get tob
SpeciesTree::SpeciesTree(Tree *input, Dict *dict, weight_t alpha, weight_t beta, bool enable_split_test) {
    std::cout << "Contracting branches with alpha = " << alpha << " and beta = " << beta << std::endl;
    if (enable_split_test) {
        std::cout << "Also consider branches with nonzero split match and mismatch count as false positive for refinement" << std::endl;
    }
    this->dict = dict;
    
    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    
    input->get_bipartitions(&internal, &bips);
    
    std::cout << bips.size() << " branches to test" << std::endl;
    std::unordered_set<Node *> false_positive;
    
    for (index_t i = 0; i < internal.size(); i ++) {
        if (internal[i]->min_pvalue < alpha || internal[i]->max_pvalue > beta) {
            false_positive.insert(internal[i]);
        } else if (enable_split_test && internal[i]->split_match_count > 0 && internal[i]->split_mismatch_count > 0) {
            std::cout << "Branch id " << i << " has nonzero split match and mismatch count, consider it as false positive for refinement. split match count: " << internal[i]->split_match_count << " split mismatch count: " << internal[i]->split_mismatch_count << std::endl;
            false_positive.insert(internal[i]);
        }
    }
    
    root = build_refinement(input->root, false_positive);
}
///////////above is for net-cs //////////////////////////////////

// Map the requested pair (a,b) to its coordinate in the qCF array. After the
// four taxa are sorted, qCF coordinates are 01|23, 02|13, and 03|12. This
// also fixes the complementary pair, so it identifies the full quartet split.
// This function only maps coordinates; it does not infer a topology.
static int topology_for_pair(const index_t *sorted, index_t a, index_t b) {
    const auto is_pair = [a, b](index_t u, index_t v) {
        return (a == u && b == v) || (a == v && b == u);
    };
    if (is_pair(sorted[0], sorted[1]) || is_pair(sorted[2], sorted[3])) return 0;
    if (is_pair(sorted[0], sorted[2]) || is_pair(sorted[1], sorted[3])) return 1;
    if (is_pair(sorted[0], sorted[3]) || is_pair(sorted[1], sorted[2])) return 2;
    return -1;
}

// Parse an oracle spec such as "t1", "t1+cf" or "t1+cf|maj". Terms joined by
// '+' are conjunctions within one track; tracks joined by '|' are alternatives,
// each carrying its own row threshold and corroboration requirement.
OracleSpec OracleSpec::parse(const std::string &text, weight_t query_alpha,
                             weight_t cf_max, weight_t margin,
                             std::string *error) {
    OracleSpec spec;
    spec.query_alpha = query_alpha;
    spec.cf_max = cf_max;
    spec.margin = margin;

    std::size_t pos = 0;
    while (pos <= text.size()) {
        const std::size_t bar = text.find('|', pos);
        const std::string track_text =
            text.substr(pos, bar == std::string::npos ? std::string::npos : bar - pos);
        OracleTerm term = {false, false, false, false, false};
        std::size_t tpos = 0;
        while (tpos <= track_text.size()) {
            const std::size_t plus = track_text.find('+', tpos);
            std::string name = track_text.substr(
                tpos, plus == std::string::npos ? std::string::npos : plus - tpos);
            // Trim, so "t1 + cf" parses like "t1+cf".
            while (!name.empty() && std::isspace((unsigned char) name.front())) name.erase(name.begin());
            while (!name.empty() && std::isspace((unsigned char) name.back())) name.pop_back();
            if (name == "t1") term.use_t1 = true;
            else if (name == "t3") term.use_t3 = true;
            else if (name == "cf") term.use_cf = true;
            else if (name == "maj") term.use_majority = true;
            else if (name == "sym") term.use_sym = true;
            else if (!name.empty()) {
                if (error) *error = "unknown oracle test: " + name;
                return spec;
            }
            if (plus == std::string::npos) break;
            tpos = plus + 1;
        }
        if (!term.use_t1 && !term.use_t3 && !term.use_majority && !term.use_sym) {
            // "cf" alone is a guard, not a test; it can only narrow something.
            if (error) *error = "each oracle track needs t1, t3, sym or maj";
            return spec;
        }
        spec.tracks.push_back(term);
        if (bar == std::string::npos) break;
        pos = bar + 1;
    }
    if (spec.tracks.empty() && error) *error = "empty oracle spec";
    return spec;
}

std::string OracleSpec::to_string() const {
    std::string out;
    for (std::size_t i = 0; i < tracks.size(); ++i) {
        if (i) out += "|";
        std::string term;
        if (tracks[i].use_t1) term += "t1";
        if (tracks[i].use_t3) term += term.empty() ? "t3" : "+t3";
        if (tracks[i].use_sym) term += term.empty() ? "sym" : "+sym";
        if (tracks[i].use_cf) term += term.empty() ? "cf" : "+cf";
        if (tracks[i].use_majority) term += term.empty() ? "maj" : "+maj";
        out += term;
    }
    return out;
}

// Statistical stand-in for the oracle predicate Query(x,y,rho,r) == xy|rho r.
// Returns true when the 4-set agrees with the target quartet under `track` of
// `oracle`, false when it contradicts it.
//
// The "t1" track is the original surrogate, unchanged: T1 tests this specific
// quartet tree, so its concordant qCF count is passed first, and false means T1
// rejected the target at query_alpha. A rejection can reflect a 4-blob, a
// different tree topology, or sampling error; it is not by itself a unique
// diagnosis of a 4-blob, and this finite-sample surrogate can also fail to
// detect one.
//
// The optional "cf" guard and "maj" test are described on OracleSpec. The qCF
// counts and T1 p-values are cached per 4-set, so a track that does not need T1
// never pays for the R call -- which also makes "t1+cf" cheaper than "t1",
// since a highly concordant quartet short-circuits before the p-value.
bool SpeciesTree::query_pairs_together(std::vector<Tree *> &input, index_t x,
                                       index_t y, index_t rho, index_t r,
                                       const OracleSpec &oracle,
                                       std::size_t track) {
    index_t sorted[4] = {x, y, rho, r};
    std::sort(sorted, sorted + 4);
    for (int i = 1; i < 4; ++i) {
        if (sorted[i - 1] == sorted[i]) return false;
    }

    const int target = topology_for_pair(sorted, x, y);
    if (target < 0) return false;

    const quartet_t quartet = join(sorted);
    auto qcf_it = row_sweep_qcfs_cache.find(quartet);
    if (qcf_it == row_sweep_qcfs_cache.end()) {
        weight_t counts[3] = {0, 0, 0};
        get_qCFs(input, sorted, counts);
        qcf_it = row_sweep_qcfs_cache.emplace(
            quartet, std::array<weight_t, 3>{counts[0], counts[1], counts[2]}
        ).first;
    }
    const std::array<weight_t, 3> &qcfs = qcf_it->second;
    // The oracle model always answers; its empirical surrogate cannot run T1
    // when no input gene tree resolves this four-taxon set.
    const weight_t total = qcfs[0] + qcfs[1] + qcfs[2];
    if (total == 0) return false;//todo, is this valid?

    const OracleTerm &term = oracle.tracks[track < oracle.tracks.size() ? track : 0];

    if (term.use_majority) {
        const weight_t best = std::max(qcfs[0], std::max(qcfs[1], qcfs[2]));
        if (!((best - qcfs[target]) / total > oracle.margin)) return true;
        if (!term.use_t1 && !term.use_t3 && !term.use_sym) return false;
    }

    // A strongly concordant quartet cannot be evidence against the split, no
    // matter how significant its minor asymmetry is.
    if (term.use_cf && qcfs[target] / total >= oracle.cf_max) return true;

    if (term.use_sym) {
        // Exchangeability of the two minors is implied by "e is a cut edge"
        // on its own -- the near-side lineages that fail to coalesce cross e
        // as an interchangeable pair -- so asymmetry between them is evidence
        // against the split, and it is the reticulation signature itself.
        // Chi-square with 1 df, whose survival function is erfc(sqrt(x/2));
        // closed form, so this track never calls into R.
        const weight_t m1 = qcfs[(target + 1) % 3];
        const weight_t m2 = qcfs[(target + 2) % 3];
        const weight_t minor_total = m1 + m2;
        bool asymmetric = false;
        if (minor_total > 0) {
            const double d = (double) (m1 - m2);
            const double chi = d * d / (double) minor_total;
            asymmetric = std::erfc(std::sqrt(chi / 2.0)) < oracle.query_alpha;
        }
        if (!asymmetric) return true;                  // symmetric: agrees
        if (!term.use_t1 && !term.use_t3) return false;
    }

    if (term.use_t3) {
        // T3 sorts the counts internally, so its p-value depends only on the
        // 4-set, not on which topology is the target: one cache slot each.
        auto t3_it = row_sweep_t3_pvalues_cache.find(quartet);
        if (t3_it == row_sweep_t3_pvalues_cache.end()) {
            weight_t counts[3] = {qcfs[0], qcfs[1], qcfs[2]};
            t3_it = row_sweep_t3_pvalues_cache.emplace(
                quartet, (weight_t) pvalue(counts)).first;
        }
        if (t3_it->second >= oracle.query_alpha) return true;   // a tree does fit
        if (!term.use_t1 && !term.use_sym) return false;
    }

    if (!term.use_t1) return true;

    auto p_it = row_sweep_t1_pvalues_cache.find(quartet);
    if (p_it == row_sweep_t1_pvalues_cache.end()) {
        const weight_t missing = std::numeric_limits<weight_t>::quiet_NaN();
        p_it = row_sweep_t1_pvalues_cache.emplace(
            quartet, std::array<weight_t, 3>{missing, missing, missing}
        ).first;
    }
    weight_t &p = p_it->second[target];
    if (std::isnan(p)) {
        // T1 treats entry 0 as the specified species-tree topology; the two
        // discordant topology counts may appear in either remaining order.
        weight_t t1_counts[3] = {
            qcfs[target], qcfs[(target + 1) % 3], qcfs[(target + 2) % 3]
        };
        p = pvalue_t1(t1_counts);
    }

    // MSCquartets marks p < alpha as a rejection of the specified T1 tree.
    return p >= oracle.query_alpha;
}

// The anchors the sweep uses: `count` taxa spread evenly through R, always
// including R[0] so that count == 1 is the original single-anchor sweep.
// Duplicates are dropped, so the result may be shorter than `count`.
static std::vector<std::size_t> row_sweep_anchors(std::size_t m,
                                                  unsigned long int count) {
    std::vector<std::size_t> out;
    if (count == 0 || count >= m) {
        for (std::size_t i = 0; i < m; ++i) out.push_back(i);
        return out;
    }
    for (unsigned long int i = 0; i < count; ++i) {
        std::size_t a = (std::size_t) ((i * m) / count);
        if (std::find(out.begin(), out.end(), a) == out.end()) out.push_back(a);
    }
    return out;
}

// RowSweepTest(A, B) -> true = ACCEPT, false = REJECT
//
// Takes taxon indices directly instead of Node* pointers, so it can be driven
// from a plain-text bipartition file without needing a base tree.
//
// With the defaults (oracle "t1", one track, tau = (1+delta)/4, heavy = 1,
// anchors = 1) this is exactly the algorithm as originally published: fix
// rho = R[0], count for each row (x, r) how many y in S\{x} contradict, and
// reject as soon as one row exceeds tau*(s-1).
//
// The generalisations, each inert at its default:
//   - several oracle tracks, each with its own tau and heavy, rejecting if any
//     track fires; a noisy indicator can then be admitted under a stricter rule
//     than a reliable one;
//   - heavy[t] > 1 demands bad rows at that many distinct partners r, which
//     rejects the single anomalous (rho, r) pair that manufactures most false
//     rejections, at the cost of detections on edges whose blob signal reaches
//     only one partner;
//   - anchors > 1 repeats the sweep from several taxa of R.
bool SpeciesTree::row_sweep_test_idx(std::vector<Tree *> &input,
                                     const std::vector<index_t> &A,
                                     const std::vector<index_t> &B,
                                     const RowSweepParams &params) {
    const std::vector<index_t> &S = (A.size() >= B.size()) ? A : B;  // larger side
    const std::vector<index_t> &R = (A.size() >= B.size()) ? B : A;  // smaller side
    const std::size_t s = S.size();
    const std::size_t m = R.size();
    if (s < 2 || m < 2) return true;   // not testable -> ACCEPT (trivial split)

    const std::size_t n_tracks = params.oracle.tracks.size();

    // Random/Pooled give every partner column its own anchor, so that a single
    // anomalous taxon pair no longer contaminates every row of the edge. The
    // published rule (Fixed) conditions every row on the one pair (R[0], r).
    if (params.row_mode != RowSweepRowMode::Fixed) {
        // Deterministic per edge: the same bipartition always draws the same
        // anchors, so a run is reproducible from the seed alone.
        std::uint64_t h = params.seed;
        for (index_t v : S) h = h * 1000003u + (std::uint64_t) v + 1u;
        for (index_t v : R) h = h * 1000003u + (std::uint64_t) v + 1u;
        std::mt19937_64 rng(h);

        // One column per partner r in R, each with an anchor drawn from R\{r}.
        std::vector<index_t> rho_of(m);
        for (std::size_t j = 0; j < m; ++j) {
            std::size_t a = (std::size_t) (rng() % (m - 1));
            if (a >= j) a++;                       // skip r itself
            rho_of[j] = R[a];
        }

        for (std::size_t t = 0; t < n_tracks; ++t) {
            const double tau = (double) params.tau[t];
            const unsigned long int heavy =
                params.heavy[t] == 0 ? 1 : params.heavy[t];

            std::vector<std::vector<std::size_t>> c(s, std::vector<std::size_t>(m, 0));
            for (std::size_t xi = 0; xi < s; xi++) {
                for (std::size_t yi = xi + 1; yi < s; yi++) {
                    for (std::size_t j = 0; j < m; j++) {
                        if (!query_pairs_together(input, S[xi], S[yi], rho_of[j],
                                                  R[j], params.oracle, t)) {
                            c[xi][j]++;
                            c[yi][j]++;
                        }
                    }
                }
            }

            if (params.row_mode == RowSweepRowMode::Pooled) {
                // One row per x over every column: a single bad pair now
                // contributes at most 1/m of the statistic instead of all of it.
                for (std::size_t xi = 0; xi < s; xi++) {
                    std::size_t total = 0;
                    for (std::size_t j = 0; j < m; j++) total += c[xi][j];
                    if ((double) total / (double) ((s - 1) * m) > tau) return false;
                }
            } else {
                std::size_t corroborating = 0;
                for (std::size_t j = 0; j < m; j++) {
                    for (std::size_t xi = 0; xi < s; xi++) {
                        if ((double) c[xi][j] / (double) (s - 1) > tau) {
                            corroborating++;
                            break;
                        }
                    }
                }
                const std::size_t needed = std::min((std::size_t) heavy, m);
                if (corroborating >= needed && corroborating > 0) return false;
            }
        }
        return true;
    }

    const std::vector<std::size_t> anchors = row_sweep_anchors(m, params.anchors);

    for (std::size_t ai = 0; ai < anchors.size(); ++ai) {
        const index_t rho = R[anchors[ai]];

        // The partners of this anchor, in R's order minus the anchor itself.
        std::vector<index_t> rest;
        rest.reserve(m - 1);
        for (std::size_t j = 0; j < m; ++j)
            if (j != anchors[ai]) rest.push_back(R[j]);
        const std::size_t r_free = rest.size();

        for (std::size_t t = 0; t < n_tracks; ++t) {
            const double tau = (double) params.tau[t];
            const unsigned long int heavy =
                params.heavy[t] == 0 ? 1 : params.heavy[t];

            std::vector<std::vector<std::size_t>> c(s, std::vector<std::size_t>(r_free, 0));
            for (std::size_t xi = 0; xi < s; xi++) {
                for (std::size_t yi = xi + 1; yi < s; yi++) {
                    for (std::size_t ri = 0; ri < r_free; ri++) {
                        // Count the empirical surrogate for q != xy|rho r. An
                        // exact 4-blob answer belongs here; with observations it
                        // is counted when the oracle track detects it. No-data
                        // queries are also treated as false.
                        if (!query_pairs_together(input, S[xi], S[yi], rho,
                                                  rest[ri], params.oracle, t)) {
                            c[xi][ri]++;
                            c[yi][ri]++;
                        }
                    }
                }
            }

            // A partner corroborates when at least one of its rows is bad; the
            // requirement is capped at the partners available, so an edge whose
            // smaller side holds two taxa behaves as it always did.
            std::size_t corroborating = 0;
            for (std::size_t ri = 0; ri < r_free; ri++) {
                for (std::size_t xi = 0; xi < s; xi++) {
                    if ((double) c[xi][ri] / (double) (s - 1) > tau) {
                        corroborating++;
                        break;
                    }
                }
            }
            const std::size_t needed = std::min((std::size_t) heavy, r_free);
            if (corroborating >= needed && corroborating > 0) return false;  // REJECT
        }
    }
    return true;                                                             // ACCEPT
}

// The continuous statistic behind row_sweep_test_idx, for the Fixed row mode.
//
// The rule rejects when `corroborating >= heavy`, where a partner corroborates
// when at least one of its rows exceeds tau. Writing v_r for a partner's
// largest row fraction and v_(k) for the k-th largest of those,
//
//     #{r : v_r > tau} >= k   <=>   v_(k) > tau,
//
// so the heavy-th largest partner maximum *is* the quantity the rule compares
// against tau, and taking the max over anchors reproduces the OR across them.
// Returning it lets a second stage re-examine an edge at a different threshold
// without recomputing anything: every 4-set it needs is already cached.
std::vector<weight_t> SpeciesTree::row_sweep_scores_idx(
        std::vector<Tree *> &input, const std::vector<index_t> &A,
        const std::vector<index_t> &B, const RowSweepParams &params) {
    const std::vector<index_t> &S = (A.size() >= B.size()) ? A : B;
    const std::vector<index_t> &R = (A.size() >= B.size()) ? B : A;
    const std::size_t s = S.size();
    const std::size_t m = R.size();
    const std::size_t n_tracks = params.oracle.tracks.size();
    std::vector<weight_t> score(n_tracks, (weight_t) 0);
    if (s < 2 || m < 2 || params.row_mode != RowSweepRowMode::Fixed)
        return std::vector<weight_t>();

    const std::vector<std::size_t> anchors = row_sweep_anchors(m, params.anchors);
    for (std::size_t ai = 0; ai < anchors.size(); ++ai) {
        const index_t rho = R[anchors[ai]];
        std::vector<index_t> rest;
        rest.reserve(m - 1);
        for (std::size_t j = 0; j < m; ++j)
            if (j != anchors[ai]) rest.push_back(R[j]);
        const std::size_t r_free = rest.size();

        for (std::size_t t = 0; t < n_tracks; ++t) {
            const unsigned long int heavy =
                params.heavy[t] == 0 ? 1 : params.heavy[t];
            std::vector<std::vector<std::size_t>> c(s, std::vector<std::size_t>(r_free, 0));
            for (std::size_t xi = 0; xi < s; xi++)
                for (std::size_t yi = xi + 1; yi < s; yi++)
                    for (std::size_t ri = 0; ri < r_free; ri++)
                        if (!query_pairs_together(input, S[xi], S[yi], rho,
                                                  rest[ri], params.oracle, t)) {
                            c[xi][ri]++;
                            c[yi][ri]++;
                        }

            std::vector<double> per_partner(r_free, 0.0);
            for (std::size_t ri = 0; ri < r_free; ri++) {
                double best = 0.0;
                for (std::size_t xi = 0; xi < s; xi++)
                    best = std::max(best, (double) c[xi][ri] / (double) (s - 1));
                per_partner[ri] = best;
            }
            std::sort(per_partner.begin(), per_partner.end(), std::greater<double>());
            const std::size_t k = std::min((std::size_t) heavy, r_free);
            if (k >= 1 && (weight_t) per_partner[k - 1] > score[t])
                score[t] = (weight_t) per_partner[k - 1];
        }
    }
    return score;
}

// Annotate every edge in a binary refinement with the row-sweep decision and
// return a copy with rejected edges contracted. ACCEPT is stored as qtt_p=1
// and REJECT as qtt_p=0 on the supplied refinement. Row sweep has no quartet
// star test, so qst_p is fixed at zero.
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict,
                         SpeciesTree* display, const RowSweepParams &row_sweep) {
    std::cout << "Constructing tree of blobs using row-sweep search"
              << " (oracle " << row_sweep.oracle.to_string() << ")" << std::endl;

    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();

    this->dict = display->dict;
    display->refine();

    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    display->get_bipartitions(&internal, &bips);
    std::cout << bips.size() << " branches to test" << std::endl;

    std::unordered_set<Node *> false_positive;
    for (std::size_t i = 0; i < internal.size(); ++i) {
        Node *edge = internal[i];
        edge->blob_id = i;
        std::cout << "Testing branch id " << i << ", ";

        std::vector<index_t> A;
        std::vector<index_t> B;
        A.reserve(bips[i].first.size());
        B.reserve(bips[i].second.size());
        for (Node *leaf : bips[i].first) A.push_back(leaf->index);
        for (Node *leaf : bips[i].second) B.push_back(leaf->index);

        bool accept = false;
        if (edge->isfake) {
            std::cout << "fake ***" << std::endl;
        } else {
            accept = row_sweep_test_idx(input, A, B, row_sweep);
            std::cout << "row-sweep: " << (accept ? "ACCEPT" : "REJECT")
                      << "; qtt_p: " << (accept ? 1 : 0) << std::endl;
        }
        if (!accept) false_positive.insert(edge);
        edge->min_pvalue = accept ? 1.0 : 0.0;
        edge->max_pvalue = 0.0;
        edge->min_f[0] = edge->min_f[1] = edge->min_f[2] = 0.0;
        edge->split_match_count = 0;
        edge->split_mismatch_count = 0;

        // Keep the generic p-value Newick representation well-formed. These
        // four taxa identify the tested edge but are not a row-sweep minimizer.
        std::vector<index_t> witness;
        witness.reserve(A.size() + B.size());
        if (A.size() >= 2 && B.size() >= 2) {
            witness = {A[0], A[1], B[0], B[1]};
        } else {
            witness.insert(witness.end(), A.begin(), A.end());
            witness.insert(witness.end(), B.begin(), B.end());
        }
        if (witness.empty()) witness.push_back(0);
        for (std::size_t j = 0; j < 4; ++j)
            edge->minimizer[j] = witness[std::min(j, witness.size() - 1)];

    }

    // Optional: record the continuous statistic for every testable edge, so a
    // whole tau curve can be swept offline from a single run instead of one run
    // per threshold. Costs no queries -- every 4-set is already cached -- and
    // runs after all decisions are made, so it cannot affect them.
    if (!row_sweep.score_out.empty()) {
        std::ofstream fout(row_sweep.score_out);
        if (fout.fail()) {
            std::cout << "\nERROR: Unable to write to " << row_sweep.score_out
                      << std::endl;
            exit(1);
        }
        fout << "edge_id\tisfake\trejected\tscores\tside\n";
        for (std::size_t i = 0; i < internal.size(); ++i) {
            Node *edge = internal[i];
            std::vector<index_t> A, B;
            for (Node *leaf : bips[i].first) A.push_back(leaf->index);
            for (Node *leaf : bips[i].second) B.push_back(leaf->index);
            std::vector<weight_t> score;
            if (!edge->isfake) score = row_sweep_scores_idx(input, A, B, row_sweep);
            fout << i << "\t" << (edge->isfake ? 1 : 0) << "\t"
                 << (false_positive.count(edge) ? 1 : 0) << "\t";
            for (std::size_t t = 0; t < score.size(); ++t)
                fout << (t ? "," : "") << score[t];
            if (score.empty()) fout << "NA";
            fout << "\t";
            const std::vector<index_t> &side = (A.size() <= B.size()) ? A : B;
            for (std::size_t j = 0; j < side.size(); ++j)
                fout << (j ? "," : "") << dict->index2label(side[j]);
            fout << "\n";
        }
        fout.close();
        std::cout << "Wrote per-edge row-sweep scores to " << row_sweep.score_out
                  << std::endl;
    }

    // Second stage: re-examine, at a looser threshold, only the edges adjacent
    // to one the first stage rejected. Blob-internal edges of T' form connected
    // subtrees (a degree-k blob contributes k-3 of them), so that neighbourhood
    // is enriched ~4x for them relative to true edges, which is what pays for
    // the looser bar. Every 4-set is already cached, so this costs no queries.
    bool two_stage = false;
    for (std::size_t t = 0; t < row_sweep.tau2.size() && t < row_sweep.tau.size(); ++t)
        if (row_sweep.tau2[t] < row_sweep.tau[t]) two_stage = true;
    if (two_stage && row_sweep.row_mode == RowSweepRowMode::Fixed) {
        std::unordered_set<Node *> stage1 = false_positive;
        std::size_t promoted = 0, considered = 0;
        for (std::size_t i = 0; i < internal.size(); ++i) {
            Node *edge = internal[i];
            if (edge->isfake || stage1.count(edge)) continue;
            // Neighbours share an endpoint: the parent, the two children, and
            // the sibling. Only internal edges of the refinement are tested.
            bool near = false;
            if (edge->parent != NULL && stage1.count(edge->parent)) near = true;
            for (Node *c : edge->children)
                if (stage1.count(c)) near = true;
            if (edge->parent != NULL)
                for (Node *c : edge->parent->children)
                    if (c != edge && stage1.count(c)) near = true;
            if (!near) continue;
            ++considered;

            std::vector<index_t> A, B;
            for (Node *leaf : bips[i].first) A.push_back(leaf->index);
            for (Node *leaf : bips[i].second) B.push_back(leaf->index);
            std::vector<weight_t> score =
                row_sweep_scores_idx(input, A, B, row_sweep);
            for (std::size_t t = 0; t < score.size(); ++t) {
                const weight_t tau2 = t < row_sweep.tau2.size()
                    ? row_sweep.tau2[t] : row_sweep.tau[t];
                if (score[t] > tau2) {
                    if (!false_positive.count(edge)) {
                        false_positive.insert(edge);
                        edge->min_pvalue = 0.0;
                        ++promoted;
                        std::cout << "Stage 2 REJECT branch id " << i
                                  << " (score " << score[t] << " > tau2 "
                                  << tau2 << ")" << std::endl;
                    }
                    break;
                }
            }
        }
        std::cout << "Stage 2 examined " << considered << " branches next to a"
                  << " stage-1 rejection and rejected " << promoted
                  << " of them" << std::endl;
    }

    if (display->root->children.size() == 2)
        false_positive.insert(display->root->children[1]);

    root = build_refinement(display->root, false_positive);
}

// Corner sets of the internal edge e = (edge, parent(edge)) of the binary
// refinement T'. On success corners[0], corners[1] are A_1(e), A_2(e) -- the two
// components at the lower endpoint -- and corners[2], corners[3] are B_1(e),
// B_2(e), the two components at the upper endpoint. The rooted encoding stores a
// degree-2 root; when the upper endpoint is that root, the true other endpoint of
// e is the root of the sibling subtree, so the upper corners come from the
// sibling's two children. Returns false when e is not an internal edge of the
// unrooted tree underlying T' (one of the four corners would be empty), i.e. when
// e induces a trivial split of X.
bool SpeciesTree::corner_sets_for_edge(Tree *refinement, Node *edge,
                                       std::vector<index_t> corners[4]) {
    for (int i = 0; i < 4; ++i) corners[i].clear();
    if (edge == NULL || edge->parent == NULL || edge->children.size() != 2)
        return false;

    Node *corner_root[4] = {edge->children[0], edge->children[1], NULL, NULL};
    bool complement_last = false;

    Node *upper = edge->parent;
    if (upper->parent == NULL && upper->children.size() == 2) {
        Node *sibling = (upper->children[0] == edge) ? upper->children[1]
                                                     : upper->children[0];
        if (sibling->children.size() != 2) return false;  // pendant, trivial split
        corner_root[2] = sibling->children[0];
        corner_root[3] = sibling->children[1];
    } else if (upper->parent == NULL) {
        if (upper->children.size() != 3) return false;    // not binary
        int slot = 2;
        for (Node *child : upper->children)
            if (child != edge) corner_root[slot++] = child;
    } else {
        if (upper->children.size() != 2) return false;    // not binary
        corner_root[2] = (upper->children[0] == edge) ? upper->children[1]
                                                      : upper->children[0];
        complement_last = true;   // B_2 is everything above the upper endpoint
    }

    std::unordered_set<index_t> covered;
    const int explicit_corners = complement_last ? 3 : 4;
    for (int i = 0; i < explicit_corners; ++i) {
        std::unordered_set<Node *> leaf_set;
        refinement->get_leaf_set(corner_root[i], &leaf_set);
        for (Node *leaf : leaf_set) {
            corners[i].push_back(leaf->index);
            covered.insert(leaf->index);
        }
        if (corners[i].empty()) return false;
    }
    if (complement_last) {
        std::vector<Node *> leaves;
        refinement->get_leaves(refinement->root, &leaves);
        for (Node *leaf : leaves)
            if (covered.find(leaf->index) == covered.end())
                corners[3].push_back(leaf->index);
        if (corners[3].empty()) return false;
    }

    // get_leaf_set returns an unordered set; sorting makes the row populations
    // -- and hence the Phase 1 samples drawn from the seed -- reproducible.
    for (int i = 0; i < 4; ++i) std::sort(corners[i].begin(), corners[i].end());
    return true;
}

// Draw a uniformly random k-element subset of {0,...,N-1} without materializing
// the population, by running k steps of a Fisher-Yates shuffle over a sparse
// representation of the identity permutation. Positions below the current step
// are never read again, so only the entries actually touched are stored.
static void corner_row_sample_indices(std::uint64_t population,
                                      std::size_t k,
                                      std::mt19937_64 &rng,
                                      std::vector<std::uint64_t> &out) {
    out.clear();
    if (population == 0) return;
    if ((std::uint64_t) k > population) k = (std::size_t) population;
    out.reserve(k);

    std::unordered_map<std::uint64_t, std::uint64_t> permuted;
    const auto value_at = [&permuted](std::uint64_t position) {
        auto it = permuted.find(position);
        return it == permuted.end() ? position : it->second;
    };
    for (std::size_t step = 0; step < k; ++step) {
        std::uniform_int_distribution<std::uint64_t> pick(step, population - 1);
        const std::uint64_t position = pick(rng);
        const std::uint64_t drawn = value_at(position);
        permuted[position] = value_at(step);
        out.push_back(drawn);
    }
}

// The cached qCF triple of a 4-set, computing it if this is the first request.
// Shares row_sweep_qcfs_cache with query_pairs_together, so a 4-set already
// queried as a contradiction test costs nothing here.
bool SpeciesTree::quartet_counts(std::vector<Tree *> &input, index_t a,
                                 index_t b, index_t c, index_t d,
                                 weight_t out[3]) {
    index_t sorted[4] = {a, b, c, d};
    std::sort(sorted, sorted + 4);
    for (int i = 1; i < 4; ++i) if (sorted[i - 1] == sorted[i]) return false;
    const quartet_t quartet = join(sorted);
    auto it = row_sweep_qcfs_cache.find(quartet);
    if (it == row_sweep_qcfs_cache.end()) {
        weight_t counts[3] = {0, 0, 0};
        get_qCFs(input, sorted, counts);
        it = row_sweep_qcfs_cache.emplace(
            quartet, std::array<weight_t, 3>{counts[0], counts[1], counts[2]}).first;
    }
    for (int i = 0; i < 3; ++i) out[i] = it->second[i];
    return out[0] + out[1] + out[2] > 0;
}

// Pooled cross-corner resolution score. Extracted verbatim from the branch-cut
// constructor so that corner row and the row sweep can use it too: it tests
// whether T' picked the locally best of the three corner pairings, which is the
// ASTRAL-error half of N, and is independent of the cluster certificate that
// finds blobs. Every 4-set goes through quartet_counts, so it shares the qCF
// cache with the oracle and costs no R call.
bool SpeciesTree::corner_resolution_score(std::vector<Tree *> &input,
                                          const std::vector<index_t> *corners,
                                          unsigned long int samples,
                                          std::mt19937_64 &rng,
                                          double *margin, double *z,
                                          std::size_t *pooled,
                                          std::vector<double> *per_set,
                                          std::vector<std::array<index_t, 4>> *per_quad) {
    if (per_set) per_set->clear();
    if (per_quad) per_quad->clear();
    if (corners[0].empty() || corners[1].empty()
        || corners[2].empty() || corners[3].empty()) return false;

    double w[3] = {0.0, 0.0, 0.0};
    const std::uint64_t pop =
        (std::uint64_t) corners[0].size() * corners[1].size()
        * corners[2].size() * corners[3].size();
    std::size_t want = (std::size_t) pop;
    if (samples > 0) {
        // A multiple of the corner-size SUM, not of one side: the population is
        // a fourfold product, so a per-side budget would grow cubically.
        const std::uint64_t budget = (std::uint64_t) samples
            * (std::uint64_t) (corners[0].size() + corners[1].size()
                               + corners[2].size() + corners[3].size());
        if (budget < pop) want = (std::size_t) budget;
    }
    std::vector<std::uint64_t> picks;
    corner_row_sample_indices(pop, want, rng, picks);
    std::size_t used = 0;
    for (std::uint64_t code : picks) {
        const index_t d0 = corners[3][code % corners[3].size()];
        code /= corners[3].size();
        const index_t c0 = corners[2][code % corners[2].size()];
        code /= corners[2].size();
        const index_t b0 = corners[1][code % corners[1].size()];
        code /= corners[1].size();
        const index_t a0 = corners[0][code % corners[0].size()];
        weight_t counts[3];
        if (!quartet_counts(input, a0, b0, c0, d0, counts)) continue;
        index_t sorted4[4] = {a0, b0, c0, d0};
        std::sort(sorted4, sorted4 + 4);
        const int t0 = topology_for_pair(sorted4, a0, b0);
        const int t1 = topology_for_pair(sorted4, a0, c0);
        const int t2 = topology_for_pair(sorted4, a0, d0);
        if (t0 < 0 || t1 < 0 || t2 < 0) continue;
        w[0] += counts[t0];
        w[1] += counts[t1];
        w[2] += counts[t2];
        ++used;
        // Per-4-set margin. The POOLED margin is a location statistic and can
        // only see an edge whose resolution is wrong throughout -- an ASTRAL
        // error. A blob-internal edge is a MIXTURE: Lem. branch-restriction
        // makes only the 4-sets that meet four distinct branches at the blob
        // vertex anomalous, and the rest look exactly like a retained edge. A
        // mixture moves the DISPERSION of the per-4-set margins, not their
        // mean, which is why pooling hides it.
        if (per_set) {
            const double s0 = (double) counts[t0], s1 = (double) counts[t1],
                         s2 = (double) counts[t2];
            const double tot = s0 + s1 + s2;
            if (tot > 0.0) {
                per_set->push_back((s0 - std::max(s1, s2)) / tot);
                if (per_quad)
                    per_quad->push_back(std::array<index_t, 4>{a0, b0, c0, d0});
            }
        }
    }
    const double wbest = std::max(w[1], w[2]);
    const double wsum = w[0] + w[1] + w[2];
    if (wsum <= 0) return false;
    if (z) *z = (w[0] - wbest) / std::sqrt(std::max(w[0] + wbest, 1.0));
    if (margin) *margin = (w[0] - wbest) / wsum;
    if (pooled) *pooled = used;
    return true;
}

// Every edge-induced cluster of T', as leaf-index bitmasks. Removing an edge of
// T' splits its leaves in two, and both sides are clusters; the branch-cut
// certificate is one particular cluster (the reticulation branch), so the
// algorithm has to be able to name all of them.
void SpeciesTree::branch_cut_clusters(Node *root, std::size_t nwords,
                                      std::vector<std::vector<std::uint64_t>> &out,
                                      std::vector<std::uint64_t> &universe) {
    std::vector<std::vector<std::uint64_t>> below;
    // Post-order: a node's leaf set is the union of its children's.
    std::vector<Node *> stack{root}, order;
    while (!stack.empty()) {
        Node *v = stack.back();
        stack.pop_back();
        order.push_back(v);
        for (Node *c : v->children) stack.push_back(c);
    }
    std::unordered_map<Node *, std::size_t> slot;
    below.assign(order.size(), std::vector<std::uint64_t>(nwords, 0));
    for (std::size_t i = 0; i < order.size(); ++i) slot[order[i]] = i;
    for (std::size_t i = order.size(); i-- > 0;) {
        Node *v = order[i];
        if (v->children.empty()) {
            below[i][v->index >> 6] |= (std::uint64_t) 1 << (v->index & 63);
        } else {
            for (Node *c : v->children) {
                const std::vector<std::uint64_t> &cb = below[slot[c]];
                for (std::size_t w = 0; w < nwords; ++w) below[i][w] |= cb[w];
            }
        }
    }
    universe = below[slot[root]];
    out.clear();
    for (std::size_t i = 0; i < order.size(); ++i) {
        if (order[i] == root) continue;
        out.push_back(below[i]);
        std::vector<std::uint64_t> comp(nwords);
        for (std::size_t w = 0; w < nwords; ++w) comp[w] = universe[w] & ~below[i][w];
        out.push_back(comp);
    }
}

static inline bool bits_subset(const std::vector<std::uint64_t> &sub,
                               const std::vector<std::uint64_t> &sup) {
    for (std::size_t w = 0; w < sub.size(); ++w)
        if (sub[w] & ~sup[w]) return false;
    return true;
}

static inline bool bits_equal(const std::vector<std::uint64_t> &a,
                              const std::vector<std::uint64_t> &b) {
    for (std::size_t w = 0; w < a.size(); ++w) if (a[w] != b[w]) return false;
    return true;
}

static inline bool bits_any(const std::vector<std::uint64_t> &a) {
    for (std::uint64_t w : a) if (w) return true;
    return false;
}

static inline bool bits_test(const std::vector<std::uint64_t> &a, index_t i) {
    return (a[i >> 6] >> (i & 63)) & 1;
}

// Localised margin contrast.
//
// The contradiction certificate and the pooled resolution margin sit at the two
// extremes of the same statistic. The certificate is cluster-aware but throws
// away everything except the bit "did this 4-set reject at query_alpha", so an
// edge whose four-cycle quarnets are never individually significant is SILENT --
// and 80% of the blob edges branch cut misses are silent, with a cluster rate of
// exactly zero. The pooled margin keeps the continuous signal but averages over
// all 4-sets, and Lem. branch-restriction says only the 4-sets that meet four
// distinct branches at the blob vertex are affected, so pooling dilutes a
// mixture into its mean. The dispersion of the per-4-set margins sees the
// mixture but not WHERE it is, so it pays the full variance of an unlabelled
// two-component problem.
//
// This statistic keeps both: for every cluster U of T' that lies properly inside
// one corner, split the sampled cross-corner 4-sets by whether that corner's
// taxon falls in U, and compare the two groups' margins. If U is the
// reticulation branch, Cor. branch-corner-alignment puts it inside exactly one
// corner and the "in" group is precisely the four-cycle family, so its margins
// drop while the "out" group stays at the retained-edge level. On a retained
// edge no cluster separates anything, whatever the branch length -- a short
// branch lowers ALL the margins together, which is the failure mode a fixed bar
// on the pooled margin cannot survive.
//
// Reports the maximum over clusters of the mean difference (`out` minus `in`)
// and of the corresponding Welch t.
static void localised_margin_contrast(
        const std::vector<std::vector<std::uint64_t>> &clusters,
        const std::vector<std::uint64_t> corner_mask[4],
        const std::vector<std::array<index_t, 4>> &quads,
        const std::vector<double> &mus,
        std::size_t min_group,
        double *best_d, double *best_t) {
    *best_d = 0.0;
    *best_t = 0.0;
    if (quads.size() != mus.size() || mus.size() < 2 * min_group) return;
    for (const std::vector<std::uint64_t> &U : clusters) {
        if (!bits_any(U)) continue;
        int corner = -1;
        for (int c = 0; c < 4; ++c) {
            if (bits_subset(U, corner_mask[c]) && !bits_equal(U, corner_mask[c])) {
                corner = c;
                break;
            }
        }
        if (corner < 0) continue;
        std::size_t n_in = 0, n_out = 0;
        double s_in = 0.0, s_out = 0.0, q_in = 0.0, q_out = 0.0;
        for (std::size_t j = 0; j < mus.size(); ++j) {
            const double v = mus[j];
            if (bits_test(U, quads[j][corner])) {
                ++n_in; s_in += v; q_in += v * v;
            } else {
                ++n_out; s_out += v; q_out += v * v;
            }
        }
        if (n_in < min_group || n_out < min_group) continue;
        const double m_in = s_in / (double) n_in, m_out = s_out / (double) n_out;
        const double v_in = std::max(0.0, q_in / (double) n_in - m_in * m_in)
                            * (double) n_in / (double) (n_in - 1);
        const double v_out = std::max(0.0, q_out / (double) n_out - m_out * m_out)
                             * (double) n_out / (double) (n_out - 1);
        const double d = m_out - m_in;
        if (d > *best_d) *best_d = d;
        // Floor the standard error at one part in a thousand of the margin
        // scale: two groups that happen to be internally constant would
        // otherwise give an unbounded t on a handful of 4-sets.
        const double se = std::sqrt(v_in / (double) n_in + v_out / (double) n_out
                                    + 1e-6);
        const double t = d / se;
        if (t > *best_t) *best_t = t;
    }
}

// Per-cluster version of the contrast above: for every cluster of T' that lies
// properly inside one corner, the gap between the mean per-4-set margin off the
// cluster and on it, indexed to match `clusters`. A cluster that does not fit
// inside a corner, or that leaves too few 4-sets on either side, gets -infinity
// so that no rule built on this can fire there.
static void per_cluster_margin_gap(
        const std::vector<std::vector<std::uint64_t>> &clusters,
        const std::vector<std::uint64_t> corner_mask[4],
        const std::vector<std::array<index_t, 4>> &quads,
        const std::vector<double> &mus,
        std::size_t min_group,
        std::vector<double> &out) {
    out.assign(clusters.size(), -std::numeric_limits<double>::infinity());
    if (quads.size() != mus.size() || mus.size() < 2 * min_group) return;
    for (std::size_t ui = 0; ui < clusters.size(); ++ui) {
        const std::vector<std::uint64_t> &U = clusters[ui];
        if (!bits_any(U)) continue;
        int corner = -1;
        for (int c = 0; c < 4; ++c) {
            if (bits_subset(U, corner_mask[c]) && !bits_equal(U, corner_mask[c])) {
                corner = c;
                break;
            }
        }
        if (corner < 0) continue;
        std::size_t n_in = 0, n_out = 0;
        double s_in = 0.0, s_out = 0.0;
        for (std::size_t j = 0; j < mus.size(); ++j) {
            if (bits_test(U, quads[j][corner])) { ++n_in; s_in += mus[j]; }
            else { ++n_out; s_out += mus[j]; }
        }
        if (n_in < min_group || n_out < min_group) continue;
        out[ui] = s_out / (double) n_out - s_in / (double) n_in;
    }
}



// P(X >= c) for X ~ Hypergeometric(N draws in total, K of them marked, n of them
// crossing), i.e. Fisher's exact upper tail for the 2x2 table
//
//               contradicts   agrees
//   crosses U        c        n - c
//   does not       K - c    N-K-n+c
//
// This is the null of BranchCutParams::contrast_alpha: "a contradiction is no
// more likely on a coordinate crossing U than on any other coordinate of this
// side". Conditioning on both margins removes the edge's own contradiction rate,
// which is the nuisance parameter that makes the absolute bar tau fragile on
// short branches.
//
// Evaluated by starting from log pmf(c) and walking up with the exact ratio
//   pmf(x+1)/pmf(x) = (K-x)(n-x) / ((x+1)(N-K-n+x+1)),
// so the cost is one lgamma quadruple plus O(tail length) multiplications.
static double hypergeom_upper_tail(std::size_t N, std::size_t K,
                                   std::size_t n, std::size_t c) {
    if (c == 0) return 1.0;
    if (K > N || n > N) return 1.0;
    const std::size_t xmax = std::min(n, K);
    if (c > xmax) return 0.0;
    // Lower end of the support: x >= n - (N - K).
    const std::size_t xmin = (n + K > N) ? (n + K - N) : 0;
    if (c <= xmin) return 1.0;
    auto lchoose = [](double a, double b) {
        return std::lgamma(a + 1.0) - std::lgamma(b + 1.0) - std::lgamma(a - b + 1.0);
    };
    const double lognorm = lchoose((double) N, (double) n);
    double term = std::exp(lchoose((double) K, (double) c)
                           + lchoose((double) (N - K), (double) (n - c))
                           - lognorm);
    double total = term;
    for (std::size_t x = c; x < xmax; ++x) {
        const double num = (double) (K - x) * (double) (n - x);
        const double den = (double) (x + 1) * (double) (N - K - n + x + 1);
        if (den <= 0.0) break;
        term *= num / den;
        total += term;
        if (term < 1e-300) break;
    }
    if (total > 1.0) total = 1.0;
    if (!(total >= 0.0)) total = 0.0;
    return total;
}

// The bar a cluster has to clear, once the localised margin gap is allowed to
// move it. `cluster_margin` splits the clusters into corroborated and not:
//   corroborated (gap > cluster_margin)         -> the lowered bar `tau_low`
//   contradicted (gap finite, <= cluster_margin)-> the raised bar `tau_high`
//   uninformative (too few 4-sets either side)  -> the published bar `tau`
// The third case matters: a singleton reticulation branch cannot be split
// against enough 4-sets to have a gap at all, and it must not be penalised for
// that. Either of `tau_low` / `tau_high` left at 0 falls back to `tau`.
static inline double cluster_bar(double tau, const BranchCutParams &bc,
                                 const std::vector<double> &cgap,
                                 std::size_t ui, std::size_t track) {
    // Track 0 only. The corroborating quantity is a depressed cross-corner
    // margin, which is what a FOUR-CYCLE quarnet looks like, so it corroborates
    // the certificate track (t1/sym) and nothing else. `maj` fires when the
    // target quartet loses its plurality -- most often a short branch -- so a
    // cluster it flags is not a reticulation branch and must not inherit the
    // lowered bar. Same argument as the seed restriction in Alg.
    // seed-and-propagate.
    if (track != 0) return tau;
    if (bc.cluster_margin == 0.0 || ui >= cgap.size()) return tau;
    const double g = cgap[ui];
    if (!std::isfinite(g)) return tau;
    if (g > (double) bc.cluster_margin)
        return bc.tau_low != 0.0 ? (double) bc.tau_low : tau;
    return bc.tau_high != 0.0 ? (double) bc.tau_high : tau;
}

// BranchCutContraction (Algorithm "branch-cut-contraction"): contract the
// refinement-only edges of T' using fully bad cut certificates.
//
// For each internal edge e with corners A_1, A_2 | B_1, B_2, sample coordinates
// from Xi_A(e) = C(A,2) x B_1 x B_2 and Xi_B(e) = C(B,2) x A_1 x A_2, then scan
// every edge-induced cluster U of T' inside A (resp. B). A coordinate *crosses*
// U when exactly one member of its near-side pair lies in U. If e is
// refinement-only and U is the reticulation branch at the blob vertex, every
// crossing coordinate is a perfect contradiction; if e is retained, no
// coordinate contradicts except through oracle error.
//
// See BranchCutParams for the two departures from the algorithm as published
// (a tunable tau below 1/2, and a minimum-support guard), and why the recorded
// evidence requires them.
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict,
                         SpeciesTree *display,
                         const BranchCutParams &branch_cut) {
    std::cout << "Constructing tree of blobs using branch-cut contraction"
              << " (oracle " << branch_cut.oracle.to_string()
              << ", tau ";
    for (std::size_t t = 0; t < branch_cut.tau.size(); ++t)
        std::cout << (t ? "," : "") << branch_cut.tau[t];
    std::cout << ", min support " << branch_cut.min_support << ")" << std::endl;

    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();

    this->dict = display->dict;
    display->refine();

    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    display->get_bipartitions(&internal, &bips);
    std::cout << internal.size() << " branches to test" << std::endl;

    const std::size_t ntaxa = (std::size_t) display->dict->size();
    const std::size_t nwords = (ntaxa + 63) / 64;
    std::vector<std::vector<std::uint64_t>> clusters;
    std::vector<std::uint64_t> universe;
    display->branch_cut_clusters(display->root, nwords, clusters, universe);

    std::mt19937_64 rng(branch_cut.seed);
    std::unordered_set<Node *> false_positive;

    // Seed-and-propagate (theory.tex Alg. seed-and-propagate). The screening
    // rule below scores an edge by the MAXIMUM over ~2n clusters, so its null
    // is inflated by that maximisation -- measured on n15+n25, the 99th
    // percentile of the statistic on retained edges is 0.080 for the max but
    // only 0.0030 for a single fixed cluster, a 27x tighter tail. Lem.
    // propagation says why: for a retained edge Phi(f, L) is null for ANY L
    // properly inside one side, so once L is known no maximisation is needed.
    //
    // So: phase 1 screens at `tau` and remembers the cluster that fired, phase 2
    // re-tests every surviving edge against those seeds at the much lower
    // `propagate_tau`. Since the seed rate never exceeds the max, phase 2 can
    // only fire where propagate_tau < rate <= tau, making it a strict addition
    // to the screening rule rather than a retuning of it.
    const bool propagate = branch_cut.propagate_tau != 0.0;
    struct EdgeEvidence {
        std::vector<std::pair<index_t, index_t>> pairs;
        std::vector<char> bad;
        std::vector<std::uint64_t> nearmask;
    };
    std::vector<std::vector<EdgeEvidence>> evidence(propagate ? internal.size() : 0);
    std::vector<int> seed_of(propagate ? internal.size() : 0, -1);
    std::vector<char> screened(propagate ? internal.size() : 0, 0);

    // The contrast test and the score dump both need every eligible cluster, so
    // neither can use the "stop at the first firing cluster" short circuit.
    const bool contrast = branch_cut.contrast_alpha > 0.0;
    // Number of split-sample corroboration groups; 1 = the published rule.
    const std::size_t G = (branch_cut.corroborate > 1)
        ? std::min<std::size_t>((std::size_t) branch_cut.corroborate, 64) : 1;
    // How many of the Geff groups have to vote to reject. frac = 1 is the
    // conjunction; smaller values give the two-level majority test.
    const double corr_frac = (branch_cut.corroborate_frac > 0.0
                              && branch_cut.corroborate_frac <= 1.0)
        ? (double) branch_cut.corroborate_frac : 1.0;
    // The bar one group must clear; defaults to the pooled tau of its track.
    const double corr_bar_opt = (double) branch_cut.corroborate_bar;
    auto inner_bar = [&](double tau_t) {
        return corr_bar_opt > 0.0 ? corr_bar_opt : tau_t;
    };
    auto need_groups = [&](std::size_t Ge) -> std::size_t {
        std::size_t k = (std::size_t) std::ceil(corr_frac * (double) Ge - 1e-9);
        if (k < 1) k = 1;
        if (k > Ge) k = Ge;
        return k;
    };
    if (G > 1)
        std::cout << "branch cut: split-sample corroboration, " << G
                  << " groups (a cluster must fire in every group)" << std::endl;
    const bool dumping = !branch_cut.score_out.empty()
                         || !branch_cut.quad_out.empty();
    const bool quad_dumping = !branch_cut.quad_out.empty();
    std::vector<std::vector<std::array<index_t, 4>>> dump_quads(
        quad_dumping ? internal.size() : 0);
    std::vector<std::vector<double>> dump_mus(
        quad_dumping ? internal.size() : 0);
    const bool collect = contrast || dumping || propagate;
    // Per-edge record for the offline sweep: one block per (track, side).
    struct SideRecord {
        int track, side;
        std::size_t M_all, C_all;
        std::vector<std::uint32_t> ui, M, C;
    };
    std::vector<std::vector<SideRecord>> dump(dumping ? internal.size() : 0);
    std::vector<std::array<double, 4>> dump_res(
        dumping ? internal.size() : 0, {0.0, 0.0, 0.0, 0.0});
    std::vector<std::array<double, 7>> dump_disp(
        dumping ? internal.size() : 0, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    std::vector<std::array<double, 2>> dump_mlc(
        dumping ? internal.size() : 0, {0.0, 0.0});
    if (contrast)
        std::cout << "branch-cut: hypergeometric contrast on, "
                  << "Bonferroni level " << branch_cut.contrast_alpha << std::endl;

    // Per-taxon anchor quality.  Ass. protocol still holds: these 4-sets are
    // drawn before any edge is tested and without reference to any edge.
    //
    // Under the pure MSC a quartet's frequencies depend only on its internal
    // branch, so anchor IDENTITY would not matter at all.  It matters here
    // because the gene trees are ESTIMATED: a taxon that is placed unreliably
    // depresses the resolution margin of every 4-set containing it.  The score
    // is that taxon's mean per-4-set margin, which needs no true topology.
    std::vector<double> anchor_q;
    if (branch_cut.anchor_quality > 0 && ntaxa >= 4) {
        std::vector<double> sum(ntaxa, 0.0);
        std::vector<std::uint32_t> cnt(ntaxa, 0);
        std::mt19937_64 qrng(branch_cut.seed ^ 0x51F3A9C7D2B84E61ULL);
        std::uniform_int_distribution<std::size_t> pick(0, ntaxa - 1);
        std::size_t used = 0;
        for (std::size_t s = 0; s < (std::size_t) branch_cut.anchor_quality; ++s) {
            index_t q[4];
            for (int k = 0; k < 4; ) {
                const index_t cand = (index_t) pick(qrng);
                bool dup = false;
                for (int j = 0; j < k; ++j) if (q[j] == cand) { dup = true; break; }
                if (!dup) q[k++] = cand;
            }
            weight_t counts[3];
            if (!quartet_counts(input, q[0], q[1], q[2], q[3], counts)) continue;
            const double c0 = (double) counts[0], c1 = (double) counts[1],
                         c2 = (double) counts[2];
            const double tot = c0 + c1 + c2;
            if (tot <= 0.0) continue;
            const double hi = std::max(c0, std::max(c1, c2));
            const double lo = std::min(c0, std::min(c1, c2));
            const double margin = (hi - (tot - hi - lo)) / tot;
            for (int k = 0; k < 4; ++k) { sum[q[k]] += margin; ++cnt[q[k]]; }
            ++used;
        }
        anchor_q.assign(ntaxa, 1.0);
        double mean = 0.0; std::size_t seen = 0;
        for (std::size_t t = 0; t < ntaxa; ++t)
            if (cnt[t] > 0) { anchor_q[t] = sum[t] / cnt[t]; mean += anchor_q[t]; ++seen; }
        if (seen > 0) mean /= (double) seen;
        // Normalise to mean 1, then floor at a small positive weight: a poor
        // anchor is down-weighted but never excluded, so every parallel class
        // stays reachable and Lem. cycle-cover's availability count is intact.
        double wmin = 1e9, wmax = -1e9;
        const double gam = branch_cut.anchor_power;
        for (std::size_t t = 0; t < ntaxa; ++t) {
            double v = (mean > 0.0) ? anchor_q[t] / mean : 1.0;
            if (gam != 1.0) v = std::pow(std::max(v, 1e-6), gam);
            if (v < 0.05) v = 0.05;
            anchor_q[t] = v;
            wmin = std::min(wmin, v);
            wmax = std::max(wmax, v);
        }
        std::cout << "branch cut: anchor quality from " << used
                  << " global 4-sets, gamma " << gam << ", weight range ["
                  << wmin << ", " << wmax << "]" << std::endl;
        if (!branch_cut.anchor_quality_out.empty()) {
            std::ofstream qf(branch_cut.anchor_quality_out);
            qf << "taxon\tname\tmean_margin\tn_sets\tweight\n";
            for (std::size_t t = 0; t < ntaxa; ++t)
                qf << t << '\t' << display->dict->index2label((index_t) t) << '\t'
                   << (cnt[t] ? sum[t] / cnt[t] : 0.0) << '\t' << cnt[t]
                   << '\t' << anchor_q[t] << '\n';
        }
    }

    for (std::size_t i = 0; i < internal.size(); ++i) {
        Node *edge = internal[i];
        edge->blob_id = i;
        std::cout << "Testing branch id " << i << ", ";

        if (edge->isfake) {
            std::cout << "fake ***" << std::endl;
            false_positive.insert(edge);
            edge->min_pvalue = 0.0;
            edge->max_pvalue = 0.0;
            edge->min_f[0] = edge->min_f[1] = edge->min_f[2] = 0.0;
            edge->split_match_count = edge->split_mismatch_count = 0;
            for (std::size_t j = 0; j < 4; ++j) edge->minimizer[j] = 0;
            continue;
        }

        std::vector<index_t> corners[4];
        bool reject = false;
        weight_t worst = 0.0;
        double best_p = 1.0;
        const std::size_t n_tracks = branch_cut.oracle.tracks.size();
        // CLUSTER-LEVEL CORROBORATION.
        //
        // Threshold work on the cluster statistic alone is finished: every
        // restriction of it trades about one rescued true edge for three lost
        // blob edges, because the maximisation over ~2n clusters isolates four
        // noise contradictions on a true edge as cleanly as it isolates
        // twenty-two real ones on a blob edge. What the statistic cannot supply
        // is a SECOND, independent read on the same cluster.
        //
        // There is one available. The cycle-cover coordinates and the
        // cross-corner 4-sets of the resolution test are different query
        // families, and Lem. reticulation-cut predicts the same cluster for
        // both: if U is the reticulation branch then the coordinates crossing U
        // contradict AND the 4-sets whose corner taxon lies in U have depressed
        // margins. A noise cluster on a true edge has no reason to do the
        // second. So the bar on the contradiction rate can be lowered to
        // `tau_low` for clusters whose localised margin gap clears
        // `cluster_margin`, which is a strict addition to the published rule.
        std::vector<double> cgap;
        const bool want_cmargin = branch_cut.cluster_margin != 0.0
                                  && (branch_cut.tau_low != 0.0
                                      || branch_cut.tau_high != 0.0);
        bool cgap_ready = false;
        // Computed on first demand: only a cluster whose count already sits in
        // the band between the lowest and highest bar can have its verdict
        // changed by the gap, and on most edges no cluster gets close, so the
        // extra cross-corner sample is never drawn.
        auto ensure_cgap = [&]() {
            if (cgap_ready) return;
            cgap_ready = true;
            std::mt19937_64 cm_rng;
            std::uint64_t z = (std::uint64_t) branch_cut.seed
                + 0x9E3779B97F4A7C15ULL * (std::uint64_t) (i + 1)
                + 0xD6E8FEB86659FD93ULL;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            cm_rng.seed(z ^ (z >> 31));
            std::mt19937_64 &cdraw = branch_cut.fixed_streams ? cm_rng : rng;
            double dummy_m = 0.0, dummy_z = 0.0;
            std::size_t dummy_p = 0;
            std::vector<double> cm_mus;
            std::vector<std::array<index_t, 4>> cm_quads;
            if (!corner_resolution_score(input, corners, branch_cut.samples,
                                         cdraw, &dummy_m, &dummy_z, &dummy_p,
                                         &cm_mus, &cm_quads))
                return;
            std::vector<std::uint64_t> corner_mask[4];
            for (int c = 0; c < 4; ++c) {
                corner_mask[c].assign(nwords, 0);
                for (index_t t : corners[c])
                    corner_mask[c][t >> 6] |= (std::uint64_t) 1 << (t & 63);
            }
            per_cluster_margin_gap(clusters, corner_mask, cm_quads, cm_mus,
                                   (std::size_t) branch_cut.mlc_min_group, cgap);
        };
        // The lowest bar any cluster can be held to, so the lazy computation
        // above is only triggered by a cluster that could actually move.
        double min_bar = 1.0;
        if (want_cmargin) {
            min_bar = branch_cut.tau.empty() ? 1.0 : (double) branch_cut.tau[0];
            if (branch_cut.tau_low != 0.0)
                min_bar = std::min(min_bar, (double) branch_cut.tau_low);
        }
        if (display->corner_sets_for_edge(display, edge, corners)) {
          for (std::size_t track = 0; track < n_tracks && (!reject || propagate || dumping); ++track) {
            const double tau_t = (double) branch_cut.tau[
                track < branch_cut.tau.size() ? track : 0];
            const int mode_t = (track < branch_cut.mode.size())
                ? branch_cut.mode[track] : 0;
            // Two symmetric sides: pairs from `near`, anchors one per far corner.
            for (int side = 0; side < 2 && (!reject || propagate || dumping); ++side) {
                // One stream per (edge, track, side) when asked for, so that a
                // flag which changes a decision cannot reshuffle the samples of
                // every later edge and turn a paired comparison into an unpaired
                // one. SplitMix-style mixing of the three indices with the seed.
                std::mt19937_64 local_rng;
                if (branch_cut.fixed_streams) {
                    std::uint64_t z = (std::uint64_t) branch_cut.seed
                        + 0x9E3779B97F4A7C15ULL * (std::uint64_t) (i + 1)
                        + 0xBF58476D1CE4E5B9ULL * (std::uint64_t)
                              (branch_cut.shared_coords ? 1 : track + 1)
                        + 0x94D049BB133111EBULL * (std::uint64_t) (side + 1);
                    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
                    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
                    local_rng.seed(z ^ (z >> 31));
                }
                std::mt19937_64 &draw = branch_cut.fixed_streams ? local_rng : rng;
                std::vector<index_t> near;
                const std::vector<index_t> &f1 = corners[side == 0 ? 2 : 0];
                const std::vector<index_t> &f2 = corners[side == 0 ? 3 : 1];
                near.insert(near.end(), corners[side == 0 ? 0 : 2].begin(),
                            corners[side == 0 ? 0 : 2].end());
                near.insert(near.end(), corners[side == 0 ? 1 : 3].begin(),
                            corners[side == 0 ? 1 : 3].end());
                const std::size_t m = near.size();
                if (m < 2 || f1.empty() || f2.empty()) continue;

                std::vector<std::uint64_t> nearmask(nwords, 0);
                for (index_t t : near)
                    nearmask[t >> 6] |= (std::uint64_t) 1 << (t & 63);

                const std::uint64_t npairs = (std::uint64_t) m * (m - 1) / 2;
                const std::uint64_t population =
                    npairs * (std::uint64_t) f1.size() * (std::uint64_t) f2.size();
                std::size_t want = (std::size_t) population;
                if (branch_cut.samples > 0) {
                    const std::uint64_t budget =
                        (std::uint64_t) branch_cut.samples * (std::uint64_t) m;
                    if (budget < population) want = (std::size_t) budget;
                }
                std::vector<std::pair<index_t, index_t>> pairs;
                std::vector<char> bad;
                // Split-sample corroboration: which group of cycles supplied
                // each coordinate.  See BranchCutParams::corroborate.
                std::vector<std::uint8_t> grp;
                std::size_t Geff = 1;

                if (branch_cut.cycles > 0) {
                    // Lem. cycle-cover: h edge-disjoint spanning cycles, one per
                    // anchor pair. A spanning cycle crosses every nonempty
                    // proper subset at least twice, so M_U >= 2h for EVERY
                    // cluster -- the depth the theorem's proof requires, which
                    // random sampling does not give.
                    const std::uint64_t t_A =
                        (std::uint64_t) f1.size() * (std::uint64_t) f2.size();
                    std::size_t h = (std::size_t) branch_cut.cycles;
                    // Availability (Lem. cycle-cover): h <= t*floor((r-1)/2),
                    // and for r == 2 a spanning 2-cycle needs 2h <= t. Taking
                    // one cycle per anchor pair needs only h <= t_A.
                    // The one-cycle-per-anchor-pair construction caps h at t_A,
                    // and t_A = |B1||B2| is TINY for an edge whose far corners
                    // are small -- a far cherry gives t_A = 1, so `h` collapses
                    // to 1 however large --branchcut-cycles is set. That cap,
                    // not any statistical saturation, is why h = 64 and h = 200
                    // measured the same as h = 32: for the edges that needed
                    // more depth the extra cycles were never drawn.
                    //
                    // Lem. cycle-cover itself allows h <= t_A * floor((m-1)/2),
                    // because K_m carries floor((m-1)/2) edge-disjoint spanning
                    // cycles (Walecki) inside EACH parallel class. `cycle_reuse`
                    // lifts the cap to that bound and takes the extra cycles as
                    // fresh random vertex orders rather than a Walecki
                    // decomposition; those are edge-disjoint only in
                    // expectation, so a coordinate that repeats one already
                    // drawn on this side is dropped instead of being counted
                    // twice (a repeat is the SAME 4-set, so the oracle returns
                    // the same answer and counting it twice would inflate M_U
                    // with no new information). Every cycle still crosses every
                    // cluster at least twice before de-duplication.
                    // Lem. cycle-cover's own availability bound; the
                    // one-cycle-per-anchor construction stops at t_A.
                    const std::uint64_t h_avail = (m == 2)
                        ? (t_A / 2)
                        : std::min<std::uint64_t>(
                              t_A * std::max<std::uint64_t>(1, (m - 1) / 2),
                              (std::uint64_t) 1 << 32);
                    // Only lift the cap where there are enough distinct anchor
                    // pairs for the extra depth to be independent evidence.
                    const bool may_reuse = branch_cut.cycle_reuse
                        && t_A >= (std::uint64_t) branch_cut.reuse_min_anchors;
                    std::size_t h_max = (std::size_t) (may_reuse
                        ? h_avail
                        : ((m == 2) ? (t_A / 2)
                                    : std::min<std::uint64_t>(t_A, t_A * ((m - 1) / 2))));
                    // Targeted lift: give a THIN edge just enough depth to leave
                    // the first-order regime of Prop. count-floor, without
                    // paying for the full `cycles`. Never lowers h_max.
                    if (branch_cut.min_depth > 0) {
                        const std::size_t want =
                            std::min<std::uint64_t>(branch_cut.min_depth, h_avail);
                        if (want > h_max) h_max = want;
                    }
                    if (h > h_max) h = h_max;
                    pairs.reserve(h * m);
                    bad.reserve(h * m);
                    // Only as many groups as there are cycles to split: an edge
                    // whose far corners are a cherry has t_A = 1, hence h = 1,
                    // and asking two groups of it would make the certificate
                    // uncorroborable and the edge unrejectable.
                    if (G > 1) { grp.reserve(h * m); Geff = std::min(G, h); }
                    // Walecki's decomposition is exactly edge-disjoint, so it
                    // never needs the dedup -- and never loses a crossing to it.
                    const bool use_walecki = branch_cut.walecki && m >= 3;
                    std::unordered_set<std::uint64_t> seen;
                    const bool dedup = !use_walecki
                                       && (may_reuse || branch_cut.min_depth > 0)
                                       && (std::uint64_t) h > t_A;
                    if (dedup) seen.reserve(h * m * 2);
                    std::vector<std::size_t> order(m);
                    for (std::size_t i = 0; i < m; ++i) order[i] = i;
                    // Under Walecki the cycles of a class are deterministic
                    // rotations of a zig-zag on one vertex labelling, so the
                    // randomness sits in the labelling rather than in each
                    // cycle. Edge-disjointness is only ever needed WITHIN a
                    // class -- two cycles in different classes carry different
                    // anchor pairs, hence different 4-sets, whatever near-side
                    // edges they share -- so each class gets its OWN labelling.
                    // Sharing one labelling across all classes costs accuracy
                    // (measured at n50: err 4.05 -> 4.45, J -0.024) because it
                    // makes every class probe the same near-side pairs.
                    // Labellings are drawn before any query (Ass. protocol).
                    const std::size_t tA_sz = f1.size() * f2.size();
                    const std::size_t ncls = std::min(h, tA_sz);
                    // The anchor pairs actually used. The enumeration
                    // f1i = k % |B1|, f2i = (k / |B1|) % |B2| walks B1 first, so
                    // whenever |B1| >= ncls EVERY class shares b2 = B2[0]: the
                    // edge leans on one anchor taxon and its errors are
                    // perfectly correlated across classes, which is the
                    // independence Ass. independence assumes between them.
                    // Sampling the pairs uniformly costs nothing and removes it.
                    std::vector<std::size_t> classes;
                    if (branch_cut.anchor_spread && tA_sz > 0 && ncls > 0) {
                        classes.reserve(ncls);
                        if (ncls >= tA_sz) {
                            for (std::size_t k = 0; k < tA_sz; ++k) classes.push_back(k);
                        } else {
                            std::unordered_set<std::size_t> got;
                            // Quality weighting draws b1 and b2 independently
                            // with weight q(b), so a pair carries q(b1)q(b2).
                            if (!anchor_q.empty()) {
                                std::vector<double> w1(f1.size()), w2(f2.size());
                                for (std::size_t k = 0; k < f1.size(); ++k)
                                    w1[k] = anchor_q[f1[k]];
                                for (std::size_t k = 0; k < f2.size(); ++k)
                                    w2[k] = anchor_q[f2[k]];
                                std::discrete_distribution<std::size_t>
                                    d1(w1.begin(), w1.end()), d2(w2.begin(), w2.end());
                                // Bounded: a heavily skewed weight vector can
                                // keep redrawing the same few pairs, so fall
                                // through to the uniform top-up below.
                                const std::size_t guard = 64 * ncls + 64;
                                for (std::size_t g = 0;
                                     g < guard && classes.size() < ncls; ++g) {
                                    const std::size_t v = d1(draw)
                                                        + f1.size() * d2(draw);
                                    if (got.insert(v).second) classes.push_back(v);
                                }
                            }
                            std::uniform_int_distribution<std::size_t> d(0, tA_sz - 1);
                            while (classes.size() < ncls) {
                                const std::size_t v = d(draw);
                                if (got.insert(v).second) classes.push_back(v);
                            }
                        }
                    }
                    std::vector<std::vector<std::size_t>> perms;
                    if (use_walecki) {
                        perms.assign(ncls, order);
                        for (std::size_t p = 0; p < ncls; ++p)
                            for (std::size_t i = m; i > 1; --i) {
                                std::uniform_int_distribution<std::size_t> d(0, i - 1);
                                std::swap(perms[p][i - 1], perms[p][d(draw)]);
                            }
                    }
                    for (std::size_t c = 0; c < h; ++c) {
                        // Distinct parallel class => the cycles are edge-disjoint.
                        // Walecki walks the classes first so that cycle j > 0 of
                        // a class is reached only once every class holds one.
                        // `cls` is the class SLOT (which labelling); `anc` is the
                        // anchor pair that slot carries. They coincide unless
                        // the pairs are being sampled rather than enumerated.
                        const std::size_t cls = use_walecki ? (c % tA_sz)
                                                            : (classes.empty() ? c : (c % ncls));
                        const std::size_t anc = classes.empty() ? cls : classes[cls % ncls];
                        const std::size_t f1i = (std::size_t) (anc % f1.size());
                        const std::size_t f2i = (std::size_t) ((anc / f1.size()) % f2.size());
                        if (use_walecki) {
                            // Cycle j of K_m: hub perm[0], the other m-1 vertices
                            // on a circle, zig-zag from a rotating start
                            //   a, a+1, a-1, a+2, a-2, ...   (mod M = m-1).
                            // Rotating `a` gives floor((m-1)/2) pairwise
                            // edge-disjoint spanning cycles -- Lem. cycle-cover's
                            // availability bound, met exactly.
                            const std::size_t M = m - 1;
                            const std::size_t j = c / tA_sz;
                            const std::vector<std::size_t> &perm = perms[cls];
                            std::size_t k = 0;
                            order[k++] = perm[0];
                            order[k++] = perm[1 + (j % M)];
                            for (std::size_t t = 1; k < m; ++t) {
                                order[k++] = perm[1 + ((j + t) % M)];
                                if (k < m)
                                    order[k++] = perm[1 + ((j + M - t) % M)];
                            }
                        } else {
                        // Vertex order drawn before any query (Ass. protocol).
                        for (std::size_t i = m; i > 1; --i) {
                            std::uniform_int_distribution<std::size_t> d(0, i - 1);
                            std::swap(order[i - 1], order[d(draw)]);
                        }
                        }
                        for (std::size_t i = 0; i < m; ++i) {
                            const index_t x = near[order[i]];
                            const index_t y = near[order[(i + 1) % m]];
                            if (m == 2 && i == 1) break;   // avoid the duplicate edge
                            if (dedup) {
                                const std::uint64_t lo = std::min(x, y), hi = std::max(x, y);
                                const std::uint64_t key = (((std::uint64_t) f1i) << 48)
                                    ^ (((std::uint64_t) f2i) << 32) ^ (lo << 16) ^ hi;
                                if (!seen.insert(key).second) continue;
                            }
                            const bool agree = query_pairs_together(
                                input, x, y, f1[f1i], f2[f2i], branch_cut.oracle, track);
                            pairs.push_back(std::make_pair(x, y));
                            bad.push_back(agree ? 0 : 1);
                            if (Geff > 1)
                                grp.push_back((std::uint8_t) (c % Geff));
                        }
                    }
                } else {
                std::vector<std::uint64_t> picks;
                corner_row_sample_indices(population, want, draw, picks);

                // Decode each sampled coordinate, query it, and record the pair.
                pairs.reserve(picks.size());
                bad.reserve(picks.size());
                if (G > 1) { grp.reserve(picks.size()); Geff = G; }
                for (std::uint64_t code : picks) {
                    const std::uint64_t f2i = code % f2.size();
                    code /= f2.size();
                    const std::uint64_t f1i = code % f1.size();
                    code /= f1.size();
                    // Unrank code in [0, npairs) to an unordered pair of `near`.
                    std::uint64_t xi = 0, rem = code;
                    while (rem >= m - 1 - xi) { rem -= (m - 1 - xi); ++xi; }
                    const std::uint64_t yi = xi + 1 + rem;
                    const index_t x = near[xi], y = near[yi];
                    const bool agree = query_pairs_together(
                        input, x, y, f1[f1i], f2[f2i], branch_cut.oracle, track);
                    pairs.push_back(std::make_pair(x, y));
                    bad.push_back(agree ? 0 : 1);
                    if (Geff > 1)
                        grp.push_back((std::uint8_t) ((pairs.size() - 1) % Geff));
                }
                }

                if (propagate) {
                    EdgeEvidence ev;
                    ev.pairs = pairs;
                    ev.bad = bad;
                    ev.nearmask = nearmask;
                    evidence[i].push_back(ev);
                }

                // Side totals: the conditioning margins of the contrast test,
                // and the statistic of a `global` track.
                std::size_t M_all = pairs.size(), C_all = 0;
                for (char b : bad) C_all += (std::size_t) (b ? 1 : 0);

                if (mode_t == 1) {
                    // GLOBAL track. `maj` asks whether the target quartet keeps
                    // its plurality, which is a property of the edge and not of
                    // any cluster, so localising it only adds the maximisation
                    // penalty. One rate test on the whole side instead.
                    if (M_all >= (std::size_t) branch_cut.min_support) {
                        const weight_t rate = (weight_t) C_all / (weight_t) M_all;
                        if (rate > worst) worst = rate;
                        const std::size_t Ct = (C_all > (std::size_t) branch_cut.trim)
                            ? C_all - (std::size_t) branch_cut.trim : 0;
                        bool fires_g = (double) Ct > tau_t * (double) M_all;
                        // A global track has no cluster to corroborate, so the
                        // side total itself has to hold up in every group.
                        if (fires_g && Geff > 1) {
                            std::vector<std::uint32_t> Mg(Geff, 0), Cg(Geff, 0);
                            for (std::size_t j = 0; j < bad.size(); ++j) {
                                ++Mg[grp[j]]; Cg[grp[j]] += (std::uint32_t) bad[j];
                            }
                            const double ib = inner_bar(tau_t);
                            std::size_t voted = 0;
                            for (std::size_t g = 0; g < Geff; ++g)
                                if (Mg[g] > 0
                                    && (double) Cg[g] > ib * (double) Mg[g]) ++voted;
                            fires_g = voted >= need_groups(Geff);
                        }
                        if (fires_g) reject = true;
                    }
                    if (dumping) {
                        SideRecord rec;
                        rec.track = (int) track; rec.side = side;
                        rec.M_all = M_all; rec.C_all = C_all;
                        dump[i].push_back(rec);
                    }
                    continue;
                }

                // Cluster scan. When the contrast test or the dump is on the
                // first firing cluster cannot settle the edge, because the
                // Bonferroni factor needs the number of eligible clusters and
                // the dump needs all of them; so gather first, decide after.
                std::vector<std::uint32_t> cl_ui, cl_M, cl_C;
                std::vector<char> cl_cor;   // corroborated in every group?
                for (std::size_t ui = 0; ui < clusters.size(); ++ui) {
                    const std::vector<std::uint64_t> &U = clusters[ui];
                    if (!bits_any(U) || !bits_subset(U, nearmask)) continue;
                    if (bits_equal(U, nearmask)) continue;   // U must be proper
                    std::size_t M = 0, C = 0;
                    std::uint32_t Mg[64], Cg[64];
                    if (Geff > 1)
                        for (std::size_t g = 0; g < Geff; ++g) { Mg[g] = 0; Cg[g] = 0; }
                    for (std::size_t j = 0; j < pairs.size(); ++j) {
                        const bool in_x = bits_test(U, pairs[j].first);
                        const bool in_y = bits_test(U, pairs[j].second);
                        if (in_x == in_y) continue;          // does not cross U
                        ++M;
                        C += bad[j];
                        if (Geff > 1) {
                            const std::size_t g = grp[j];
                            ++Mg[g]; Cg[g] += (std::uint32_t) bad[j];
                        }
                    }
                    if (M < (std::size_t) branch_cut.min_support) continue;
                    if (collect) {
                        cl_ui.push_back((std::uint32_t) ui);
                        cl_M.push_back((std::uint32_t) M);
                        cl_C.push_back((std::uint32_t) C);
                        if (Geff > 1) {
                            const double ib = inner_bar(tau_t);
                            std::size_t voted = 0;
                            for (std::size_t g = 0; g < Geff; ++g)
                                if (Mg[g] > 0
                                    && (double) Cg[g] > ib * (double) Mg[g]) ++voted;
                            cl_cor.push_back(voted >= need_groups(Geff) ? 1 : 0);
                        }
                        continue;
                    }
                    // Trimmed statistic: drop up to `trim` contradicting
                    // coordinates before comparing with the threshold.
                    const std::size_t Ct = (C > (std::size_t) branch_cut.trim)
                        ? C - (std::size_t) branch_cut.trim : 0;
                    const weight_t rate = (weight_t) C / (weight_t) M;
                    if (rate > worst) {
                        worst = rate;
                        // The cluster that best explains this edge is the seed
                        // candidate phase 2 propagates. Recorded only when
                        // propagating, so the early exit below is unaffected.
                        //
                        // Only track 0 may supply a seed. Lem. seed certifies a
                        // cluster because the perfect answer there is a
                        // FOUR-CYCLE, which is what the certificate track
                        // (t1/sym) tests. `maj` is an "unsupported" detector --
                        // it fires on any quartet that loses its plurality, most
                        // often a short branch -- so a cluster it flags is not a
                        // reticulation branch and propagating it costs
                        // sensitivity with no certificate behind it.
                        if (propagate && track == 0) seed_of[i] = (int) ui;
                    }
                    if (want_cmargin && track == 0
                        && (double) Ct > min_bar * (double) M) ensure_cgap();
                    const double bar_u =
                        cluster_bar(tau_t, branch_cut, cgap, ui, track);
                    bool fires = (double) Ct > bar_u * (double) M
                                 && C >= (std::size_t) branch_cut.cmin;
                    // Corroboration: the SAME cluster must clear the bar in
                    // every group of cycles, not just in the pooled count.
                    if (fires && Geff > 1) {
                        const double ib = inner_bar(bar_u);
                        std::size_t voted = 0;
                        for (std::size_t g = 0; g < Geff; ++g)
                            if (Mg[g] > 0
                                && (double) Cg[g] > ib * (double) Mg[g]) ++voted;
                        fires = voted >= need_groups(Geff);
                    }
                    if (fires) {
                        reject = true;
                        // Without propagation the first firing cluster settles
                        // the edge; with it, keep scanning so the seed is the
                        // argmax rather than whichever cluster happened first.
                        if (!propagate) break;
                    }
                }

                if (collect) {
                    const std::size_t n_elig = cl_ui.size();
                    for (std::size_t q = 0; q < n_elig; ++q) {
                        const std::size_t M = cl_M[q], C = cl_C[q];
                        const weight_t rate = (weight_t) C / (weight_t) M;
                        if (rate > worst) {
                            worst = rate;
                            if (propagate && track == 0) seed_of[i] = (int) cl_ui[q];
                        }
                        const std::size_t Ct = (C > (std::size_t) branch_cut.trim)
                            ? C - (std::size_t) branch_cut.trim : 0;
                        if (want_cmargin && track == 0
                            && (double) Ct > min_bar * (double) M) ensure_cgap();
                        if (!((double) Ct >
                              cluster_bar(tau_t, branch_cut, cgap, cl_ui[q], track)
                              * (double) M)) continue;
                        if (C < (std::size_t) branch_cut.cmin) continue;
                        if (Geff > 1 && q < cl_cor.size() && !cl_cor[q]) continue;
                        if (contrast) {
                            // Fisher's exact upper tail against the side's own
                            // contradiction rate, Bonferroni-adjusted over the
                            // clusters this edge maximised over.
                            double p = hypergeom_upper_tail(M_all, C_all, M, C)
                                       * (double) (n_elig ? n_elig : 1);
                            if (p > 1.0) p = 1.0;
                            if (p < best_p) best_p = p;
                            if (!(p < (double) branch_cut.contrast_alpha)) continue;
                        }
                        reject = true;
                    }
                    if (dumping) {
                        SideRecord rec;
                        rec.track = (int) track; rec.side = side;
                        rec.M_all = M_all; rec.C_all = C_all;
                        rec.ui = cl_ui; rec.M = cl_M; rec.C = cl_C;
                        dump[i].push_back(rec);
                    }
                }
            }
          }
        }

        // Edge-level corner-resolution test: did T' pick the locally best of the
        // three corner pairings? Aimed at ASTRAL errors, which the cluster
        // certificate cannot see because they are not blobs at all.
        double res_z = std::numeric_limits<double>::infinity();
        double res_margin = std::numeric_limits<double>::infinity();
        const bool want_res = (branch_cut.resolution_z != 0.0
                               || branch_cut.resolution_margin != 0.0);
        // min, q10, q25, median, mean, sd, count of the per-4-set margins.
        std::array<double, 7> disp = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double mlc_d = 0.0, mlc_t = 0.0;
        const bool want_mlc = dumping || branch_cut.mlc_d != 0.0
                              || branch_cut.mlc_t != 0.0;
        const bool want_disp = dumping || branch_cut.margin_q10 != 0.0
                               || branch_cut.margin_sd != 0.0
                               || branch_cut.margin_min != 0.0 || want_mlc;
        // The dump wants the margin for every edge, including the ones the
        // cluster rule already rejected -- but only when the streams are fixed,
        // since under the shared stream an extra draw would shift every later
        // edge and the dump would no longer describe the run that produced it.
        const bool force_res = dumping && branch_cut.fixed_streams;
        if ((!reject && (want_res || want_disp)) || force_res) {
            // `rng` is shared with the cluster sampling above, so this call must
            // stay in the same place and draw the same number of values as the
            // block it replaced, or every --branchcut-blob result shifts.
            std::mt19937_64 res_local;
            if (branch_cut.fixed_streams) {
                std::uint64_t z = (std::uint64_t) branch_cut.seed
                    + 0x9E3779B97F4A7C15ULL * (std::uint64_t) (i + 1)
                    + 0xD6E8FEB86659FD93ULL;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
                res_local.seed(z ^ (z >> 31));
            }
            std::mt19937_64 &rdraw = branch_cut.fixed_streams ? res_local : rng;
            std::size_t res_pooled = 0;
            std::vector<double> mus;
            std::vector<std::array<index_t, 4>> quads;
            const bool ok = corner_resolution_score(input, corners,
                                                    branch_cut.samples, rdraw,
                                                    &res_margin, &res_z, &res_pooled,
                                                    want_disp ? &mus : nullptr,
                                                    want_mlc ? &quads : nullptr);
            if (want_mlc && !quads.empty()) {
                std::vector<std::uint64_t> corner_mask[4];
                for (int c = 0; c < 4; ++c) {
                    corner_mask[c].assign(nwords, 0);
                    for (index_t t : corners[c])
                        corner_mask[c][t >> 6] |= (std::uint64_t) 1 << (t & 63);
                }
                localised_margin_contrast(clusters, corner_mask, quads, mus,
                                          (std::size_t) branch_cut.mlc_min_group,
                                          &mlc_d, &mlc_t);
            }
            if (want_disp && mus.size() >= 2) {
                std::vector<double> srt(mus);
                std::sort(srt.begin(), srt.end());
                const std::size_t k = srt.size();
                double mean = 0.0;
                for (double v : srt) mean += v;
                mean /= (double) k;
                double var = 0.0;
                for (double v : srt) var += (v - mean) * (v - mean);
                var /= (double) (k - 1);
                disp[0] = srt[0];                                   // min
                disp[1] = srt[(std::size_t) (0.10 * (k - 1))];      // q10
                disp[2] = srt[(std::size_t) (0.25 * (k - 1))];      // q25
                disp[3] = srt[k / 2];                               // median
                disp[4] = mean;
                disp[5] = std::sqrt(var);                           // sd
                disp[6] = (double) k;
            }
            if (dumping) {
                dump_res[i] = {ok ? 1.0 : 0.0, res_margin, res_z, (double) res_pooled};
                dump_disp[i] = disp;
                dump_mlc[i] = {mlc_d, mlc_t};
                if (quad_dumping) { dump_quads[i] = quads; dump_mus[i] = mus; }
            }
            if (!reject && want_res && ok
                && res_pooled >= (std::size_t) branch_cut.resolution_min_pooled) {
                if (branch_cut.resolution_margin != 0.0
                    && res_margin < (double) branch_cut.resolution_margin)
                    reject = true;
                if (branch_cut.resolution_z != 0.0
                    && res_z < (double) branch_cut.resolution_z)
                    reject = true;
            }
            // Dispersion rules. Applied only to an edge the cluster certificate
            // and the pooled margin both accepted, so they are a strict
            // addition aimed at the blob edges those two cannot see.
            if (!reject && disp[6] >= (double) branch_cut.margin_min_sets) {
                if (branch_cut.margin_q10 != 0.0
                    && disp[1] < (double) branch_cut.margin_q10) reject = true;
                if (branch_cut.margin_sd != 0.0
                    && disp[5] > (double) branch_cut.margin_sd) reject = true;
                if (branch_cut.margin_min != 0.0
                    && disp[0] < (double) branch_cut.margin_min) reject = true;
                if (branch_cut.mlc_d != 0.0
                    && mlc_d > (double) branch_cut.mlc_d) reject = true;
                if (branch_cut.mlc_t != 0.0
                    && mlc_t > (double) branch_cut.mlc_t) reject = true;
            }
        }

        std::cout << "branch-cut: " << (reject ? "REJECT" : "ACCEPT")
                  << " (best cluster rate " << worst;
        if (contrast) std::cout << ", best contrast p " << best_p;
        std::cout << ", resolution margin " << res_margin << ")" << std::endl;
        if (propagate) screened[i] = reject ? 1 : 0;
        if (reject) false_positive.insert(edge);
        edge->min_pvalue = reject ? 0.0 : 1.0;
        edge->max_pvalue = 0.0;
        edge->min_f[0] = edge->min_f[1] = edge->min_f[2] = 0.0;
        edge->split_match_count = edge->split_mismatch_count = 0;

        std::vector<index_t> A, B;
        for (Node *leaf : bips[i].first) A.push_back(leaf->index);
        for (Node *leaf : bips[i].second) B.push_back(leaf->index);
        std::vector<index_t> witness;
        if (A.size() >= 2 && B.size() >= 2) witness = {A[0], A[1], B[0], B[1]};
        else {
            witness.insert(witness.end(), A.begin(), A.end());
            witness.insert(witness.end(), B.begin(), B.end());
        }
        if (witness.empty()) witness.push_back(0);
        for (std::size_t j = 0; j < 4; ++j)
            edge->minimizer[j] = witness[std::min(j, witness.size() - 1)];
    }

    // Phase 2: propagation. The clusters that fired during screening are the
    // candidate reticulation branches; a blob contributes a connected clump of
    // collapsed edges that all share one, so a seed found at any edge of the
    // clump certifies the rest. Testing a FIXED cluster carries none of the
    // maximisation penalty of phase 1, which is what licenses the lower bar.
    if (propagate) {
        std::vector<int> seeds;
        for (std::size_t i = 0; i < internal.size(); ++i)
            if (screened[i] && seed_of[i] >= 0) seeds.push_back(seed_of[i]);
        std::sort(seeds.begin(), seeds.end());
        seeds.erase(std::unique(seeds.begin(), seeds.end()), seeds.end());

        const double tau2 = (double) branch_cut.propagate_tau;
        std::size_t promoted = 0;
        for (std::size_t i = 0; i < internal.size(); ++i) {
            Node *edge = internal[i];
            if (edge->isfake || screened[i]) continue;
            bool reject = false;
            for (int ui : seeds) {
                const std::vector<std::uint64_t> &U = clusters[ui];
                for (const EdgeEvidence &ev : evidence[i]) {
                    if (!bits_any(U) || !bits_subset(U, ev.nearmask)) continue;
                    if (bits_equal(U, ev.nearmask)) continue;
                    std::size_t M = 0, C = 0;
                    for (std::size_t j = 0; j < ev.pairs.size(); ++j) {
                        const bool in_x = bits_test(U, ev.pairs[j].first);
                        const bool in_y = bits_test(U, ev.pairs[j].second);
                        if (in_x == in_y) continue;
                        ++M;
                        C += ev.bad[j];
                    }
                    if (M < (std::size_t) branch_cut.min_support) continue;
                    const std::size_t Ct = (C > (std::size_t) branch_cut.trim)
                        ? C - (std::size_t) branch_cut.trim : 0;
                    if ((double) Ct > tau2 * (double) M) { reject = true; break; }
                }
                if (reject) break;
            }
            if (reject) {
                ++promoted;
                false_positive.insert(edge);
                edge->min_pvalue = 0.0;
            }
        }
        std::cout << "propagation: " << seeds.size() << " seed cluster(s), "
                  << promoted << " further edge(s) contracted at tau2 = "
                  << tau2 << std::endl;
    }

    // Per-edge sufficient statistics, so a whole multi-dimensional threshold
    // sweep (tau, min support, trim, contrast level, resolution margin, per
    // track and per side) costs one run per replicate rather than one run per
    // point. Written after every decision, so it cannot affect any of them.
    if (dumping) {
        std::ofstream fout(branch_cut.score_out);
        if (fout.fail()) {
            std::cout << "\nERROR: Unable to write to " << branch_cut.score_out
                      << std::endl;
            exit(1);
        }
        fout << "edge_id\tisfake\trejected\tres_ok\tres_margin\tres_z"
                "\tres_pooled\tm_min\tm_q10\tm_q25\tm_med\tm_mean\tm_sd"
                "\tm_n\tmlc_d\tmlc_t\tevidence\tside\n";
        for (std::size_t i = 0; i < internal.size(); ++i) {
            Node *edge = internal[i];
            fout << i << "\t" << (edge->isfake ? 1 : 0) << "\t"
                 << (false_positive.count(edge) ? 1 : 0) << "\t";
            if (dump_res[i][0] != 0.0)
                fout << 1 << "\t" << dump_res[i][1] << "\t" << dump_res[i][2]
                     << "\t" << (std::size_t) dump_res[i][3];
            else
                fout << 0 << "\tNA\tNA\t0";
            for (std::size_t q = 0; q < 6; ++q) fout << "\t" << dump_disp[i][q];
            fout << "\t" << (std::size_t) dump_disp[i][6];
            fout << "\t" << dump_mlc[i][0] << "\t" << dump_mlc[i][1];
            fout << "\t";
            // evidence := block '|' block, block := track:side:M_all:C_all:list
            //             list  := ui,M,C ';' ui,M,C ...
            if (dump[i].empty()) fout << "NA";
            for (std::size_t b = 0; b < dump[i].size(); ++b) {
                const SideRecord &r = dump[i][b];
                if (b) fout << "|";
                fout << r.track << ":" << r.side << ":" << r.M_all << ":"
                     << r.C_all << ":";
                for (std::size_t q = 0; q < r.ui.size(); ++q)
                    fout << (q ? ";" : "") << r.ui[q] << "," << r.M[q] << ","
                         << r.C[q];
            }
            fout << "\t";
            std::vector<index_t> A, B;
            for (Node *leaf : bips[i].first) A.push_back(leaf->index);
            for (Node *leaf : bips[i].second) B.push_back(leaf->index);
            const std::vector<index_t> &sd = (A.size() <= B.size()) ? A : B;
            for (std::size_t j = 0; j < sd.size(); ++j)
                fout << (j ? "," : "") << dict->index2label(sd[j]);
            fout << "\n";
        }
        fout.close();
        std::cout << "Wrote per-edge branch-cut evidence to "
                  << branch_cut.score_out << std::endl;
    }

    if (quad_dumping) {
        std::ofstream fout(branch_cut.quad_out);
        if (fout.fail()) {
            std::cout << "\nERROR: Unable to write to " << branch_cut.quad_out
                      << std::endl;
            exit(1);
        }
        fout << "edge_id\ta1\ta2\tb1\tb2\tmu\n";
        for (std::size_t i = 0; i < internal.size(); ++i)
            for (std::size_t j = 0; j < dump_mus[i].size(); ++j)
                fout << i << "\t"
                     << dict->index2label(dump_quads[i][j][0]) << "\t"
                     << dict->index2label(dump_quads[i][j][1]) << "\t"
                     << dict->index2label(dump_quads[i][j][2]) << "\t"
                     << dict->index2label(dump_quads[i][j][3]) << "\t"
                     << dump_mus[i][j] << "\n";
        fout.close();
        std::cout << "Wrote per-4-set margins to " << branch_cut.quad_out
                  << std::endl;
    }

    if (display->root->children.size() == 2)
        false_positive.insert(display->root->children[1]);

    root = build_refinement(display->root, false_positive);
}

// CornerRowContraction (Algorithm "refinement-edge-sweep"): edge-by-edge
// reconstruction of the tree of blobs T from a binary refinement T'. Every
// internal edge of T' is tested with rows built only from that edge's four corner
// sets, and every flagged edge is contracted.
//
// Phase 1 pre-samples every row of every edge from a private seed before a single
// oracle query is issued, as the theorem's query-independence assumption requires;
// Phase 2 then classifies each edge, terminating as soon as one row flags it.
//
// The oracle O is replaced by the same configurable surrogate the row sweep
// uses: query_pairs_together(x, y, c_1, c_2, oracle) is false when that oracle
// spec judges the 4-set to contradict xy | c_1 c_2, which is what the algorithm
// counts as a contradiction. With the default spec "t1" that is exactly "T1
// rejects at query_alpha", i.e. the original behaviour. Corner row is the method
// most likely to profit from a spec aimed squarely at 4-blob signal, since a
// blob is all it is looking for. Query results are cached inside that surrogate, so a 4-set shared
// by several rows or edges is only evaluated once.
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict,
                         SpeciesTree *display,
                         const CornerRowParams &corner_row) {
    std::cout << "Constructing tree of blobs using corner-row contraction"
              << std::endl;

    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();

    this->dict = display->dict;
    display->refine();

    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    display->get_bipartitions(&internal, &bips);
    std::cout << internal.size() << " branches to test" << std::endl;

    std::vector<Node *> all_leaves;
    display->get_leaves(display->root, &all_leaves);
    std::vector<index_t> taxa;
    taxa.reserve(all_leaves.size());
    for (Node *leaf : all_leaves) taxa.push_back(leaf->index);
    std::sort(taxa.begin(), taxa.end());   // fix the order of the rows x in X
    const std::size_t n = taxa.size();

    // Lemma "row-population size" guarantees |Omega_x(e)| >= n - 3, so every row
    // admits a k-element subset exactly when 1 <= k <= n - 3.
    const std::size_t k_max = (n >= 4) ? (n - 3) : 0;
    std::size_t k = (std::size_t) corner_row.k;
    if (k == 0) k = k_max;
    if (k > k_max) {
        std::cout << "WARNING: row sample size k=" << k << " exceeds n-3="
                  << k_max << "; using k=" << k_max << std::endl;
        k = k_max;
    }
    // Heavy sampling spends the slack in rows whose population exceeds k: row x
    // draws min(heavy * k, |Omega_x(e)|) tuples instead of k. Note that
    // max_{j <= heavy} min(j*k, |Omega_x(e)|) collapses to min(heavy*k,
    // |Omega_x(e)|), since min(j*k, N) is nondecreasing in j. Both tails of
    // Theorem "refinement-edge-guarantee" are e^{-(row size) * D(tau||.)}, and
    // its lemmas hold for any sample size at most the row population, so a
    // larger row can only sharpen the guarantee.
    const unsigned long int heavy = (corner_row.heavy == 0) ? 1 : corner_row.heavy;
    const std::size_t n_tracks = corner_row.oracle.tracks.size();
    if (corner_row.tau.size() != n_tracks) {
        std::cout << "\nERROR: corner-row needs one tau per oracle track ("
                  << n_tracks << " expected, " << corner_row.tau.size()
                  << " given)" << std::endl;
        exit(1);
    }
    const double tau = (double) corner_row.tau[0];
    if (!(tau > 0.0 && tau < 1.0)) {
        std::cout << "\nERROR: corner-row tau must lie strictly between 0 and 1: "
                  << tau << std::endl;
        exit(1);
    }
    std::cout << "corner-row parameters: k=" << k << ", heavy=" << heavy
              << (corner_row.cross ? ", cross" : "")
              << ", tau=" << tau
              << ", query alpha=" << corner_row.query_alpha
              << ", seed=" << corner_row.seed << std::endl;

    // Separate stream for the resolution test, so that turning it on cannot
    // perturb the phase 1 row sample and every pre-existing run reproduces.
    std::mt19937_64 resolution_rng(corner_row.seed ^ 0x9E3779B97F4A7C15ULL);

    // Corner sets of every internal edge, computed up front.
    std::vector<std::array<std::vector<index_t>, 4>> edge_corners(internal.size());
    std::vector<bool> testable(internal.size(), false);
    for (std::size_t i = 0; i < internal.size(); ++i) {
        internal[i]->blob_id = i;
        if (internal[i]->isfake) continue;   // arbitrary resolution of a polytomy
        std::vector<index_t> corners[4];
        if (!corner_sets_for_edge(display, internal[i], corners)) continue;
        for (int c = 0; c < 4; ++c) edge_corners[i][c].swap(corners[c]);
        testable[i] = true;
    }

    // Phase 1: pre-sample every row of every edge, before any oracle query. A row
    // of edge e indexed by x in A draws from Omega_x^A(e) = (A \ {x}) x B_1 x B_2
    // and by x in B from Omega_x^B(e) = (B \ {x}) x A_1 x A_2.
    struct CornerRowSample {
        index_t x;
        std::vector<std::array<index_t, 3>> tuples;
    };
    std::vector<std::vector<CornerRowSample>> rows(internal.size());

    std::mt19937_64 rng(corner_row.seed);
    std::vector<std::uint64_t> drawn;
    std::size_t sampled_rows = 0;
    std::size_t sampled_tuples = 0;
    std::size_t exhaustive_rows = 0;
    for (std::size_t i = 0; i < internal.size(); ++i) {
        if (!testable[i]) continue;
        const std::vector<index_t> &A1 = edge_corners[i][0];
        const std::vector<index_t> &A2 = edge_corners[i][1];
        const std::vector<index_t> &B1 = edge_corners[i][2];
        const std::vector<index_t> &B2 = edge_corners[i][3];

        std::vector<index_t> A(A1);
        A.insert(A.end(), A2.begin(), A2.end());
        std::vector<index_t> B(B1);
        B.insert(B.end(), B2.begin(), B2.end());
        std::sort(A.begin(), A.end());
        std::sort(B.begin(), B.end());
        const std::unordered_set<index_t> A_set(A.begin(), A.end());

        const std::unordered_set<index_t> A1_set(A1.begin(), A1.end());
        const std::unordered_set<index_t> B1_set(B1.begin(), B1.end());

        rows[i].reserve(n);
        for (index_t x : taxa) {
            const bool x_in_A = A_set.find(x) != A_set.end();
            const std::vector<index_t> &same = x_in_A ? A : B;
            const std::vector<index_t> &first = x_in_A ? B1 : A1;
            const std::vector<index_t> &second = x_in_A ? B2 : A2;

            // Omega_x(e) = (A \ {x}) x B_1 x B_2 draws the partner y from the
            // whole near side, so a y in x's own corner gives a 4-set spanning
            // only three corners. Such a 4-set has a tree quarnet whatever the
            // edge is, so it can never contradict, and it dilutes the row's
            // contradiction fraction by the size of x's corner. --corner-cross
            // draws y from the opposite near corner instead, so every query
            // spans all four corners. The row population shrinks accordingly.
            std::vector<index_t> partners;
            if (corner_row.cross) {
                const std::vector<index_t> &opposite =
                    x_in_A ? (A1_set.count(x) ? A2 : A1)
                           : (B1_set.count(x) ? B2 : B1);
                partners = opposite;
            } else {
                partners.reserve(same.size() - 1);
                for (index_t y : same)
                    if (y != x) partners.push_back(y);
            }

            const std::uint64_t population = (std::uint64_t) partners.size() *
                first.size() * second.size();
            if (population == 0) continue;   // no queryable tuple for this row
            const std::uint64_t row_size =
                std::min<std::uint64_t>((std::uint64_t) heavy * k, population);
            if (row_size == population) exhaustive_rows++;
            corner_row_sample_indices(population, (std::size_t) row_size, rng,
                                      drawn);

            CornerRowSample row;
            row.x = x;
            row.tuples.reserve(drawn.size());
            for (std::uint64_t code : drawn) {
                const std::uint64_t yi = code % partners.size();
                code /= partners.size();
                const std::uint64_t ci = code % first.size();
                code /= first.size();
                row.tuples.push_back({partners[yi], first[ci], second[code]});
            }
            sampled_tuples += row.tuples.size();
            rows[i].push_back(row);
            sampled_rows++;
        }
    }
    std::cout << "Phase 1: sampled " << sampled_rows << " rows ("
              << sampled_tuples << " tuples, " << exhaustive_rows
              << " exhaustive) over " << internal.size() << " branches"
              << std::endl;

    // Phase 2: classify each edge, with early termination.
    std::unordered_set<Node *> contracted;
    for (std::size_t i = 0; i < internal.size(); ++i) {
        Node *edge = internal[i];
        std::cout << "Testing branch id " << i << ", ";

        index_t witness[4] = {0, 0, 0, 0};
        bool have_witness = false;
        bool flagged = false;

        if (edge->isfake) {
            // Introduced by refine() to break a polytomy of the input, so it is
            // not an edge of T'; contract it rather than query it.
            flagged = true;
            std::cout << "fake ***" << std::endl;
        } else if (!testable[i]) {
            // A trivial split is a split of every tree on X, so it is retained.
            std::cout << "corner-row: ACCEPT (trivial split, not tested)"
                      << std::endl;
        } else {
            // Each oracle track is swept over the same rows with its own tau,
            // and the edge is flagged if any of them finds a bad row.
            for (std::size_t t = 0; t < n_tracks && !flagged; ++t) {
                const double tau_t = (double) corner_row.tau[t];
                // How many distinct rows x must be bad before this track flags
                // the edge. At 1 the first bad row ends the sweep, which is the
                // published behaviour; above 1 a single anomalous taxon can no
                // longer contract an edge on its own -- the corner-row analogue
                // of the row sweep's --rowsweep-heavy, which demands bad rows at
                // several partners.
                const unsigned long int need =
                    (t < corner_row.corroborate.size()) ? corner_row.corroborate[t] : 1;
                unsigned long int bad_rows = 0;
                for (const CornerRowSample &row : rows[i]) {
                    // Rows may differ in size under heavy sampling, so the
                    // flagging rule c > tau * k uses each row's own sample size.
                    const double threshold = tau_t * (double) row.tuples.size();
                    std::size_t contradictions = 0;
                    for (const std::array<index_t, 3> &tuple : row.tuples) {
                        if (!query_pairs_together(input, row.x, tuple[0], tuple[1],
                                                  tuple[2], corner_row.oracle, t)) {
                            contradictions++;
                            if ((double) contradictions > threshold) {
                                // The count only grows, so the remaining queries
                                // of this row cannot change its verdict.
                                if (!have_witness) {
                                    witness[0] = row.x;
                                    witness[1] = tuple[0];
                                    witness[2] = tuple[1];
                                    witness[3] = tuple[2];
                                    have_witness = true;
                                }
                                bad_rows++;
                                break;
                            }
                        }
                    }
                    if (bad_rows >= need) { flagged = true; break; }
                }
            }
            // Corner-resolution test, orthogonal to the contradiction test
            // above: it asks whether T' picked the locally best of the three
            // corner pairings, which is what an ASTRAL error gets wrong and a
            // blob does not. Uses its own generator so the phase 1 sample, and
            // hence every pre-existing corner-row result, is untouched.
            if (!flagged && corner_row.resolution_margin != 0.0) {
                double res_margin = std::numeric_limits<double>::infinity();
                if (corner_resolution_score(input, edge_corners[i].data(),
                                            corner_row.resolution_samples,
                                            resolution_rng, &res_margin,
                                            nullptr, nullptr)
                    && res_margin < (double) corner_row.resolution_margin) {
                    flagged = true;
                }
            }
            std::cout << "corner-row: " << (flagged ? "REJECT" : "ACCEPT")
                      << "; qtt_p: " << (flagged ? 0 : 1) << std::endl;
        }

        if (flagged) contracted.insert(edge);
        edge->min_pvalue = flagged ? 0.0 : 1.0;
        edge->max_pvalue = 0.0;
        edge->min_f[0] = edge->min_f[1] = edge->min_f[2] = 0.0;
        edge->split_match_count = 0;
        edge->split_mismatch_count = 0;

        // Keep the generic p-value Newick representation well-formed. For a
        // contracted edge these are the four taxa of the query that pushed a row
        // past the threshold; otherwise they merely identify the edge.
        if (!have_witness) {
            const std::vector<Node *> &A = bips[i].first;
            const std::vector<Node *> &B = bips[i].second;
            std::vector<index_t> fallback;
            fallback.reserve(A.size() + B.size());
            if (A.size() >= 2 && B.size() >= 2) {
                fallback = {A[0]->index, A[1]->index, B[0]->index, B[1]->index};
            } else {
                for (Node *leaf : A) fallback.push_back(leaf->index);
                for (Node *leaf : B) fallback.push_back(leaf->index);
            }
            if (fallback.empty()) fallback.push_back(0);
            for (std::size_t j = 0; j < 4; ++j)
                witness[j] = fallback[std::min(j, fallback.size() - 1)];
        }
        for (std::size_t j = 0; j < 4; ++j) edge->minimizer[j] = witness[j];
    }

    std::cout << "Contracting " << contracted.size() << " of "
              << internal.size() << " branches" << std::endl;

    // Suppress the degree-2 root of the rooted encoding, as the row sweep does.
    if (display->root->children.size() == 2)
        contracted.insert(display->root->children[1]);

    root = build_refinement(display->root, contracted);
}

// Instrumentation for offline study of the row-sweep decision rule.
//
// Issues exactly the queries row_sweep_test_idx would issue on every internal
// edge of this refinement, but records the evidence instead of thresholding it:
// the qCF counts and the MSCquartets T1 p-value of every distinct 4-set. No
// decision is made and no existing code path is affected, so a dump costs the
// same as one row-sweep run and reproduces its inputs exactly.
//
// Two files are written:
//   <prefix>.meta.json     taxon labels, refined newick, and per-edge S/R/rho
//   <prefix>.quartets.bin  little-endian records, one per distinct 4-set:
//                            int32 t[4]   taxon indices, ascending
//                            float64 q[3] qCF counts for t0t1|t2t3, t0t2|t1t3,
//                                         t0t3|t1t2
//                            float64 p[3] T1 p-value of the same three
//                                         topologies; NaN when not queried
// With all_targets set, the three p-values are always filled in, which costs
// extra R calls but lets an offline rule score topologies the sweep never asked
// about.
void SpeciesTree::dump_row_sweep_evidence(std::vector<Tree *> &input,
                                          const std::string &out_prefix,
                                          const OracleSpec &oracle,
                                          bool all_targets,
                                          bool all_anchors) {
    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();
    refine();

    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    get_bipartitions(&internal, &bips);
    std::cout << bips.size() << " branches to dump" << std::endl;

    std::ofstream meta(out_prefix + ".meta.json");
    if (meta.fail()) {
        std::cout << "\nERROR: Unable to write to " << out_prefix << ".meta.json" << std::endl;
        exit(1);
    }
    meta << "{\n  \"query_alpha\": " << oracle.query_alpha
         << ",\n  \"all_targets\": " << (all_targets ? "true" : "false")
         << ",\n  \"all_anchors\": " << (all_anchors ? "true" : "false")
         << ",\n  \"ntaxa\": " << dict->size() << ",\n  \"taxa\": [";
    for (index_t idx = 0; idx < dict->size(); ++idx)
        meta << (idx ? ", " : "") << "\"" << dict->index2label(idx) << "\"";
    meta << "],\n  \"refinement_newick\": \"" << to_string_basic() << "\",\n  \"edges\": [\n";

    for (std::size_t i = 0; i < internal.size(); ++i) {
        Node *edge = internal[i];
        std::vector<index_t> A, B;
        A.reserve(bips[i].first.size());
        B.reserve(bips[i].second.size());
        for (Node *leaf : bips[i].first) A.push_back(leaf->index);
        for (Node *leaf : bips[i].second) B.push_back(leaf->index);

        // Mirror row_sweep_test_idx's choice of sides, anchor, and loop order.
        const std::vector<index_t> &S = (A.size() >= B.size()) ? A : B;
        const std::vector<index_t> &R = (A.size() >= B.size()) ? B : A;
        const bool testable = (S.size() >= 2 && R.size() >= 2);

        meta << (i ? ",\n" : "") << "    {\"id\": " << i
             << ", \"isfake\": " << (edge->isfake ? "true" : "false")
             << ", \"testable\": " << (testable ? "true" : "false")
             << ", \"S\": [";
        for (std::size_t j = 0; j < S.size(); ++j) meta << (j ? "," : "") << S[j];
        meta << "], \"R\": [";
        for (std::size_t j = 0; j < R.size(); ++j) meta << (j ? "," : "") << R[j];
        meta << "]}";

        if (!testable || edge->isfake) continue;

        const index_t rho = R[0];
        const std::size_t r_free = R.size() - 1;
        for (std::size_t xi = 0; xi < S.size(); ++xi) {
            for (std::size_t yi = xi + 1; yi < S.size(); ++yi) {
                if (all_anchors) {
                    // Every pair from R, not just the ones containing the
                    // arbitrary anchor R[0], so an offline rule can re-anchor
                    // the sweep or average over anchors. Costs a factor |R|/2.
                    for (std::size_t ai = 0; ai < R.size(); ++ai)
                        for (std::size_t bi = ai + 1; bi < R.size(); ++bi)
                            query_pairs_together(input, S[xi], S[yi], R[ai],
                                                 R[bi], oracle);
                } else {
                    for (std::size_t ri = 0; ri < r_free; ++ri) {
                        // Populates the qCF and T1 caches; the answer is discarded.
                        query_pairs_together(input, S[xi], S[yi], rho, R[ri + 1],
                                             oracle);
                    }
                }
            }
        }
        std::cout << "Dumped branch id " << i << " (|S|=" << S.size()
                  << ", |R|=" << R.size() << "); distinct 4-sets so far: "
                  << row_sweep_qcfs_cache.size() << std::endl;
    }
    meta << "\n  ]\n}\n";
    meta.close();

    std::ofstream bin(out_prefix + ".quartets.bin", std::ios::binary);
    if (bin.fail()) {
        std::cout << "\nERROR: Unable to write to " << out_prefix << ".quartets.bin" << std::endl;
        exit(1);
    }
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    std::size_t written = 0;
    for (const auto &entry : row_sweep_qcfs_cache) {
        index_t packed[4];
        split(packed, entry.first);
        std::sort(packed, packed + 4);

        auto p_it = row_sweep_t1_pvalues_cache.find(entry.first);
        double p[3] = {nan_value, nan_value, nan_value};
        if (p_it != row_sweep_t1_pvalues_cache.end()) {
            for (int j = 0; j < 3; ++j) p[j] = (double) p_it->second[j];
        }
        if (all_targets && entry.second[0] + entry.second[1] + entry.second[2] > 0) {
            for (int j = 0; j < 3; ++j) {
                if (!std::isnan(p[j])) continue;
                weight_t counts[3] = {entry.second[j], entry.second[(j + 1) % 3],
                                      entry.second[(j + 2) % 3]};
                p[j] = (double) pvalue_t1(counts);
            }
        }

        int32_t taxa[4];
        for (int j = 0; j < 4; ++j) taxa[j] = (int32_t) packed[j];
        double q[3] = {(double) entry.second[0], (double) entry.second[1],
                       (double) entry.second[2]};
        bin.write(reinterpret_cast<const char *>(taxa), sizeof(taxa));
        bin.write(reinterpret_cast<const char *>(q), sizeof(q));
        bin.write(reinterpret_cast<const char *>(p), sizeof(p));
        ++written;
    }
    bin.close();
    std::cout << "Wrote " << written << " quartet records to "
              << out_prefix << ".quartets.bin" << std::endl;
}

// Parse a comma-separated list of taxon names into indices, using the dict-derived
// name2index map built by the caller. Exits with an error on an unrecognized name,
// rather than silently minting a new index (which Dict::label2index would do).
static std::vector<index_t> row_sweep_parse_side(const std::string &field,
        const std::unordered_map<std::string, index_t> &name2index) {
    std::vector<index_t> out;
    std::stringstream ss(field);
    std::string name;
    while (std::getline(ss, name, ',')) {
        if (name.empty()) continue;
        auto it = name2index.find(name);
        if (it == name2index.end()) {
            std::cerr << "ERROR: unknown taxon in --rowsweep file: " << name << std::endl;
            exit(1);
        }
        out.push_back(it->second);
    }
    return out;
}

// Reads a TSV of "id\tA\tB" bipartitions (A, B = comma-separated taxon names),
// runs row_sweep_test_idx on each, and writes "id\tprediction" (1 = ACCEPT/split, 0 = REJECT).
void SpeciesTree::run_split_experiment(std::vector<Tree *> &input,
        const std::unordered_map<std::string, index_t> &name2index,
        const std::string &bipartition_file,
        const std::string &output_file,
        const RowSweepParams &params) {
    add_r_libpaths_and_load(RINS);                 // load MSCquartets before any pvalue() call
    for (Tree *t : input) t->LCA_preprocessing();  // required by get_quartet()
    std::ifstream fin(bipartition_file);
    if (fin.fail()) {
        std::cout << "\nERROR: Unable to open " << bipartition_file << std::endl;
        exit(1);
    }
    std::ofstream fout(output_file);
    if (fout.fail()) {
        std::cout << "\nERROR: Unable to write to " << output_file << std::endl;
        exit(1);
    }
    fout << "id\tprediction\n";
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string id, aField, bField;
        std::getline(ss, id, '\t');
        std::getline(ss, aField, '\t');
        std::getline(ss, bField, '\t');
        if (id == "id") continue;                       // header
        std::vector<index_t> A = row_sweep_parse_side(aField, name2index);
        std::vector<index_t> B = row_sweep_parse_side(bField, name2index);
        bool accept = row_sweep_test_idx(input, A, B, params);
        fout << id << "\t" << (accept ? 1 : 0) << "\n";
    }
    fin.close();
    fout.close();
    std::cout << "Wrote row-sweep predictions to " << output_file << std::endl;
}

// 2f2a search algorithm O(n^3)
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict,
                         SpeciesTree* display, QCFWriter *qcf_writer) {//display means refinement of tree of blobs
    std::cout << "Constructing tree of blobs using 2-fix-2-alter search" << std::endl;

    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();
    this->dict = display->dict;
    display->refine();

    //std::string mode = "n";
    //display->annotate(input, mode);

    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    display->get_bipartitions(&internal, &bips);
    std::cout << bips.size() << " branches to test" << std::endl;
    std::unordered_set<Node *> false_positive;

    for (std::size_t i = 0; i < internal.size(); i ++) {
        std::cout << "Testing branch id " << i << ", ";
        internal[i]->blob_id = i;

        if (internal[i]->isfake) {
            false_positive.insert(internal[i]);
            std::cout << "fake ***" << std::endl;
            continue;
        }

        // Search for min p-value for quartet tree tests
        weight_t min;
        index_t minimizer[4];
        size_t split_match_count = 0;
        size_t split_mismatch_count = 0;

        std:: cout << "testing branch with split : [";
        for (Node *n : bips[i].first) {
            std::cout << dict->index2label(n->index) << ",";
        }
        std::cout << "]" << std::endl << "and [";
        for (Node *n : bips[i].second) {
            std::cout << dict->index2label(n->index) << ",";
        }
        std::cout << "]" << std::endl;


        min = search_2f2a(input, bips[i].first, bips[i].second,
                          minimizer, split_match_count, split_mismatch_count,
                          static_cast<index_t>(i), qcf_writer);
        internal[i]->min_pvalue = min;

        // Get qCFs that yielded the min p-value
        weight_t min_f[3];
        get_qCFs(input, minimizer, min_f);
        internal[i]->min_f[0] = min_f[0];
        internal[i]->min_f[1] = min_f[1];
        internal[i]->min_f[2] = min_f[2];
        
        internal[i]->minimizer[0] = minimizer[0];
        internal[i]->minimizer[1] = minimizer[1];
        internal[i]->minimizer[2] = minimizer[2];
        internal[i]->minimizer[3] = minimizer[3];
        internal[i]->split_match_count = split_match_count;
        internal[i]->split_mismatch_count = split_mismatch_count;

        // Appy quartet star test
        weight_t max = -1.0;
        if ((min_f[0] + min_f[1] + min_f[2]) > 0)
            max = pvalue_star(min_f);
        //if ((internal[i]->f[0] + internal[i]->f[1] + internal[i]->f[2]) > 0)
        //    max = pvalue_star(internal[i]->f);
        //if (iter_limit != 0)
        //    max = search_star(input, bips[i].first, bips[i].second, iter_limit);
        //else
        //    max = search_star(input, bips[i].first, bips[i].second);
        internal[i]->max_pvalue = max;

        // Write to standard out
        std::cout << "QTT: " << min << "; ";
        std::cout << "QST: " << max << "; ";
        std::cout << "qCF: [" << min_f[0] << "/" << min_f[1] << "/" << min_f[2] << "]; ";
        std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        std::cout << "split match/mismatch count: [" << split_match_count << "/" << split_mismatch_count << "]";
        std::cout << std::endl;

        //std::cout << "QTT: " << min << " ";
        //std::cout << "qCF: [" << internal[i]->min_f[0] << "/" << internal[i]->min_f[1] << "/" << internal[i]->min_f[2] << "] ";
        //std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        //std::cout << "QST: " << max << " ";
        //std::cout << "[" << internal[i]->f[0] << "/" << internal[i]->f[1] << "/" << internal[i]->f[2] << "]";
        //std::cout << std::endl;

        // Write to a table...
    }
    if (display->root->children.size() == 2) 
        false_positive.insert(display->root->children[1]);

    root = build_refinement(display->root, false_positive);
}

// 3 fix 1 alter search algorithm O(n^2)
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, 
                         SpeciesTree* display, 
                         unsigned long int iter_limit_blob, 
                         bool three_fix_one_alter,
                         bool two_fix_two_alter, 
                         bool is_quard,
                         const std::string &output_qcfs_table_file) {
    std::unique_ptr<QCFWriter> qcf_writer;
    if (!output_qcfs_table_file.empty()) {
        qcf_writer = std::make_unique<QCFWriter>(output_qcfs_table_file, dict);
    }
    
    if (!three_fix_one_alter && !is_quard && !two_fix_two_alter) {
       SpeciesTree(input, dict, display, iter_limit_blob, qcf_writer.get());
       return;
    }

    if (two_fix_two_alter && three_fix_one_alter) {
        std::cerr << "Error: both two_fix_two_alter and three_fix_one_alter cannot be true at the same time." << std::endl;
        return;
    }

    if (two_fix_two_alter && is_quard) {
        std::cerr << "Error: both two_fix_two_alter and is_quard cannot be true at the same time." << std::endl;
        return;
    }

    if (two_fix_two_alter) {
        SpeciesTree(input, dict, display, qcf_writer.get());
        return;
    }
    
    if (is_quard){
        std::cout << "Constructing tree of blobs using quard search" << std::endl;
    } else if (three_fix_one_alter) {
        std::cout << "Constructing tree of blobs using 3-fix-1-alter search" << std::endl;
    }

    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();

    display->refine();

    this->dict = display->dict;

    //std::string mode = "n";
    //display->annotate(input, mode);
    
    std::vector<Node *> internal;
    
    std::vector<std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>>> quads;
    display->get_quardpartitions(&internal, &quads, dict);
    std::cout << quads.size() << " branches to test" << std::endl;

    std::unordered_set<Node *> false_positive;
    for (index_t i = 0; i < internal.size(); i ++) {
        std::cout << "Testing branch id " << i << ", ";
        internal[i]->blob_id = i;

        if (internal[i]->isfake) {
            false_positive.insert(internal[i]);
            std::cout << "fake ***" << std::endl;
            continue;
        }


        // Search for min p-value for quartet tree tests
        weight_t min;
        index_t minimizer[4];
        size_t match_count = 0;
        size_t mismatch_count = 0;
        if (three_fix_one_alter) {
            
            min = search_3f1a(input, &quads[i], minimizer, i, qcf_writer.get());
        }   
        else if (is_quard)
            min = search_quard(input, &quads[i], minimizer, i, qcf_writer.get());
        internal[i]->min_pvalue = min;

        // Get qCFs that yielded the min p-value
        weight_t min_f[3];
        get_qCFs(input, minimizer, min_f);
        internal[i]->min_f[0] = min_f[0];
        internal[i]->min_f[1] = min_f[1];
        internal[i]->min_f[2] = min_f[2];

        internal[i]->minimizer[0] = minimizer[0];
        internal[i]->minimizer[1] = minimizer[1];
        internal[i]->minimizer[2] = minimizer[2];
        internal[i]->minimizer[3] = minimizer[3];
        if (two_fix_two_alter) {
            internal[i]->split_match_count = match_count;
            internal[i]->split_mismatch_count = mismatch_count;
        }
        // Apply quartet star test
        weight_t max = -1.0;
        if ((min_f[0] + min_f[1] + min_f[2]) > 0)
            max = pvalue_star(min_f);
        internal[i]->max_pvalue = max;
        //if ((internal[i]->f[0] + internal[i]->f[1] + internal[i]->f[2]) > 0)
        //    max = pvalue_star(internal[i]->f);

        // Write to standard out
        std::cout << "QTT: " << min << "; ";
        std::cout << "QST: " << max << "; ";
        std::cout << "qCF: [" << min_f[0] << "/" << min_f[1] << "/" << min_f[2] << "]; ";
        std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "]";
        std::cout << std::endl;
        //std::cout << "QTT: " << min << " ";
        //std::cout << "qCF: [" << internal[i]->min_f[0] << "/" << internal[i]->min_f[1] << "/" << internal[i]->min_f[2] << "] ";
        //std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        //std::cout << "QST: " << max << " ";
        //std::cout << "[" << internal[i]->f[0] << "/" << internal[i]->f[1] << "/" << internal[i]->f[2] << "]";
        //std::cout << std::endl;
    }

    if (display->root->children.size() == 2) {
        false_positive.insert(display->root->children[1]);
    }

    this->dict = display->dict;
    root = build_refinement(display->root, false_positive);
}


// O(n*cn*klogn), O(kn^3logn) if c is O(n)
SpeciesTree::SpeciesTree(std::vector<Tree *> &input, Dict *dict, 
                         SpeciesTree* display,
                         unsigned long int iter_limit_blob,
                         QCFWriter *qcf_writer) {
    std::cout << "Constructing tree of blobs" << std::endl;
    add_r_libpaths_and_load(RINS);
    for (Tree *t : input) t->LCA_preprocessing();
    this->dict = display->dict;
    display->refine();

    //std::string mode = "n";
    //display->annotate(input, mode);

    std::vector<Node *> internal;
    std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> bips;
    display->get_bipartitions(&internal, &bips);
    std::cout << bips.size() << " branches to test" << std::endl;
    std::unordered_set<Node *> false_positive;
    size_t iter_limit = iter_limit_blob;

    for (std::size_t i = 0; i < internal.size(); i ++) {
        std::cout << "Testing branch id " << i << ", ";
        internal[i]->blob_id = i;

        if (internal[i]->isfake) {
            false_positive.insert(internal[i]);
            std::cout << "fake ***" << std::endl;
            continue;
        }

        // Search for min p-value for quartet tree tests
        weight_t min;
        index_t minimizer[4];
        if (iter_limit != 0)
            min = search(input, bips[i].first, bips[i].second, iter_limit,
                         minimizer, static_cast<index_t>(i), qcf_writer);
        else
            min = search(input, bips[i].first, bips[i].second, minimizer,
                         static_cast<index_t>(i), qcf_writer);
        internal[i]->min_pvalue = min;

        // Get qCFs that yielded the min p-value
        weight_t min_f[3];
        get_qCFs(input, minimizer, min_f);
        internal[i]->min_f[0] = min_f[0];
        internal[i]->min_f[1] = min_f[1];
        internal[i]->min_f[2] = min_f[2];
        
        internal[i]->minimizer[0] = minimizer[0];
        internal[i]->minimizer[1] = minimizer[1];
        internal[i]->minimizer[2] = minimizer[2];
        internal[i]->minimizer[3] = minimizer[3];
    

        // Appy quartet star test
        weight_t max = -1.0;
        if ((min_f[0] + min_f[1] + min_f[2]) > 0)
            max = pvalue_star(min_f);
        //if ((internal[i]->f[0] + internal[i]->f[1] + internal[i]->f[2]) > 0)
        //    max = pvalue_star(internal[i]->f);
        //if (iter_limit != 0)
        //    max = search_star(input, bips[i].first, bips[i].second, iter_limit);
        //else
        //    max = search_star(input, bips[i].first, bips[i].second);
        internal[i]->max_pvalue = max;

        // Write to standard out
        std::cout << "QTT: " << min << "; ";
        std::cout << "QST: " << max << "; ";
        std::cout << "qCF: [" << min_f[0] << "/" << min_f[1] << "/" << min_f[2] << "]; ";
        std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        std::cout << std::endl;

        //std::cout << "QTT: " << min << " ";
        //std::cout << "qCF: [" << internal[i]->min_f[0] << "/" << internal[i]->min_f[1] << "/" << internal[i]->min_f[2] << "] ";
        //std::cout << "minimizer: [" << dict->index2label(minimizer[0]) << "/" << dict->index2label(minimizer[1]) << "/" << dict->index2label(minimizer[2]) << "/" << dict->index2label(minimizer[3]) << "] ";
        //std::cout << "QST: " << max << " ";
        //std::cout << "[" << internal[i]->f[0] << "/" << internal[i]->f[1] << "/" << internal[i]->f[2] << "]";
        //std::cout << std::endl;

        // Write to a table...
    }
    if (display->root->children.size() == 2) 
        false_positive.insert(display->root->children[1]);

    root = build_refinement(display->root, false_positive);
}


std::string SpeciesTree::to_string_pvalue() {
    return display_tree_pvalue(root) + ";";
}

std::string SpeciesTree::display_tree_pvalue(Node *root) {
    if (root->children.size() == 0) 
        return dict->index2label(root->index);
    std::string s = "(";
    for (Node * node : root->children) 
        s += display_tree_pvalue(node) + ",";
    s[s.size() - 1] = ')';
    if (root->parent != NULL && (root->parent->parent != NULL || (root->parent->parent == NULL && root == root->parent->children[0]))) {
        std::ostringstream ss;
        // serialize the minimizer and qCFs
        ss << std::scientific << std::setprecision(12) 
                              << "'[" 
                              << "blob_id=" << std::to_string(root->blob_id)
                              << ";qtt_p=" << (double) root->min_pvalue
                              << ";qst_p=" << (double) root->max_pvalue 
                              << ";qcf_1=" << std::to_string((int) root->min_f[0])
                              << ";qcf_2=" << std::to_string((int) root->min_f[1])
                              << ";qcf_3=" << std::to_string((int) root->min_f[2])
                              << ";minimizer=" << dict->index2label(root->minimizer[0]) << "/" << dict->index2label(root->minimizer[1]) << "/" << dict->index2label(root->minimizer[2]) << "/" << dict->index2label(root->minimizer[3])
                             << ";split_match_count=" << root->split_match_count
                             << ";split_mismatch_count=" << root->split_mismatch_count
                              << "]'";
        return s + ss.str();
    }
    else {
        return s;
    }
}

Node *SpeciesTree::build_refinement(Node *root, std::unordered_set<Node *> false_positive) {
    if (root->children.size() == 0) {
        Node *new_root = new Node(root->index, false);
        index2node[new_root->index] = new_root;
        return new_root;
    }
    else {
        Node *new_root = new Node(pseudonym(), false_positive.find(root) != false_positive.end());
        for (Node *child : root->children) {
            Node *new_child = build_refinement(child, false_positive);
            if (! new_child->isfake) {
                new_root->children.push_back(new_child);
            }
            else {
                for (Node *grand_child : new_child->children) 
                    new_root->children.push_back(grand_child);
                new_child->children.clear();
                delete new_child;
            }
        }
        for (Node *new_child : new_root->children) 
            new_child->parent = new_root;
        return new_root;
    }
}

void Tree::get_bipartition(Node *root, std::vector<Node *> *A, std::vector<Node *> *B) {
    std::vector<Node *> leaves;
    get_leaves(this->root, &leaves);
    std::unordered_set<Node *> leaf_set;
    get_leaf_set(root, &leaf_set);
    for (Node *leaf : leaves) {
        if (leaf_set.find(leaf) == leaf_set.end()) 
            A->push_back(leaf);
        else 
            B->push_back(leaf);
    }
}

void Tree::get_bipartitions(Node *root, std::vector<Node *> *internal, std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> *bips) {
    if (root->children.size() == 0) return ;
    if (root->parent != NULL && (root->parent->parent != NULL || root->parent->children.size() > 2 || root == root->parent->children[0])) {
        std::vector<Node *> A, B;
        get_bipartition(root, &A, &B);
        bips->push_back(std::make_pair(A, B));
        internal->push_back(root); // std::cout << "?" << root->support << std::endl;
    }
    for (Node *child : root->children) {
        get_bipartitions(child, internal, bips);
    }
}

void Tree::get_bipartitions(std::vector<Node *> *internal, std::vector<std::pair<std::vector<Node *>, std::vector<Node *>>> *bips) {
    get_bipartitions(root, internal, bips);
}

std::string Tree::display_bipartition(std::vector<Node *> &A, std::vector<Node *> &B) {
    std::string s = "";
    for (Node *a : A) s += std::to_string(a->index) + " ";
    s += "| ";
    for (Node *b : B) s += std::to_string(b->index) + " ";
    return s;
}

void Tree::get_quardpartition(Node *root, std::vector<Node *> *A, std::vector<Node *> *B, std::vector<Node *> *C, std::vector<Node *> *D, Dict *dict) {
    //std::cout << "Getting quardpartition for node " << root->index << std::endl;

    std::vector<Node *> leaves;
    get_leaves(this->root, &leaves);

    if (root->children.size() < 2) {
            std::cerr << "Error: less than two children" << std::endl;
            exit(1);
        }

    if (root->children.size() > 2) {
        std::cerr << "Error: more than two children" << std::endl;
        exit(1);
    }
    
    Node *c_1, *c_2, *s_1;

    c_1 = root->children[0];
    c_2 = root->children[1];

    if (root != this->root && root->parent != this->root) {
        s_1 = root->parent->children[0] == root ? root->parent->children[1] : root->parent->children[0];

    } else if (root != this->root && root->parent == this->root){
        s_1 = (this->root->children[0] == root) ? this->root->children[1]->children[0] : this->root->children[0]->children[0];
    }

    std::unordered_set<Node *> leaf_set1;
    std::unordered_set<Node *> leaf_set2;
    std::unordered_set<Node *> leaf_set3;
    get_leaf_set(c_1, &leaf_set1);
    get_leaf_set(c_2, &leaf_set2);
    get_leaf_set(s_1, &leaf_set3);
    for (Node *leaf : leaves) {
        if (leaf_set1.find(leaf) != leaf_set1.end()) 
            A->push_back(leaf);
        else if (leaf_set2.find(leaf) != leaf_set2.end()) 
            B->push_back(leaf);
        else if (leaf_set3.find(leaf) != leaf_set3.end()) 
            C->push_back(leaf);
        else 
            D->push_back(leaf);
    }
    // std::string quard_str = display_quardpartitions(*A, *B, *C, *D, dict);
    // // checking sum of sizes 
    // if (A->size() + B->size() + C->size() + D->size() != leaves.size()) {
    //     std::cerr << "Error in quardpartition: sizes do not add up!" << std::endl;
    //     exit(1);
    // }
    // std::cout << "Quardpartition: " << quard_str << std::endl;
    
}

void Tree::get_quardpartitions(Node *root, std::vector<Node *> *internal, std::vector<std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>>> *quards, std::unordered_set<uint64_t>* seen_edges, Dict *dict) {
    
    if (root->children.size() == 0) return ;

    
    if (root->parent != NULL && (root->parent->parent != NULL || root->parent->children.size() > 2 || root == root->parent->children[0])) {
        std::vector<Node *> A, B, C, D;
        get_quardpartition(root, &A, &B, &C, &D, dict);
        auto quard = std::make_tuple(A, B, C, D); 
        internal->push_back(root); // std::cout << "?" << root->support << std::endl;
        quards->push_back(quard);
    }
    
    for (Node *child : root->children) {
        get_quardpartitions(child, internal, quards, seen_edges, dict);
    }
}

void Tree::get_quardpartitions(std::vector<Node *> *internal, 
                               std::vector<std::tuple<std::vector<Node *>,
                               std::vector<Node *>, std::vector<Node *>,std::vector<Node *>>> *quards, Dict *dict) {
    std::unordered_set<uint64_t> seen_edges;
    get_quardpartitions(root, internal, quards, &seen_edges, dict);
}


std::string Tree::display_quardpartitions(std::vector<Node *> &A, std::vector<Node *> &B, std::vector<Node *> &C, std::vector<Node *> &D, Dict *dict) {
    std::string s = "";
    for (Node *a : A) s += dict->index2label(a->index) + ",";
    s = s.substr(0, s.size() - 1);
    s += "| ";
    for (Node *b : B) s += dict->index2label(b->index) + ",";
    s = s.substr(0, s.size() - 1);
    s += "| ";
    for (Node *c : C) s += dict->index2label(c->index) + ",";
    s = s.substr(0, s.size() - 1);
    s += "| ";
    for (Node *d : D) s += dict->index2label(d->index) + ",";
    s = s.substr(0, s.size() - 1);
    return s;
}


weight_t SpeciesTree::search(std::vector<Tree *> &input, 
                             std::vector<Node *> &A, 
                             std::vector<Node *> &B, 
                             index_t* minimizer,
                             index_t branch_id,
                             QCFWriter *qcf_writer) {
    index_t i[4];
    weight_t min = -1;
    size_t count = 0;
    for (i[0] = 0; i[0] < A.size(); i[0] ++) {
        for (i[1] = i[0] + 1; i[1] < A.size(); i[1] ++) {
            for (i[2] = 0; i[2] < B.size(); i[2] ++) {
                for (i[3] = i[2] + 1; i[3] < B.size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = A[i[0]]->index;
                    temp[1] = A[i[1]]->index;
                    temp[2] = B[i[2]]->index;
                    temp[3] = B[i[3]]->index;
                    count ++;

                    auto [score, qcfs] = get_pvalue_and_qCFs(
                        input, temp, branch_id, qcf_writer
                    );
                    if (min < 0 || score < min) {
                        min = score;
                        minimizer[0] = temp[0]; 
                        minimizer[1] = temp[1]; 
                        minimizer[2] = temp[2]; 
                        minimizer[3] = temp[3];
                    }
                }
            }
        }
    }
    //std::cout << "naive iter: " << count << std::endl;
    return min;
}


weight_t SpeciesTree::search_quard(std::vector<Tree *> &input, std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>> *quad, index_t* minimizer, index_t branch_id, QCFWriter *qcf_writer) {
    index_t i[4] = {0, 0, 0, 0};
    weight_t min = -1;
    auto& t = *quad;

    for (i[0] = 0; i[0] < std::get<0>(t).size(); i[0] ++) {
        for (i[1] = 0; i[1] < std::get<1>(t).size(); i[1] ++) {
            for (i[2] = 0; i[2] < std::get<2>(t).size(); i[2] ++) {
                for (i[3] = 0; i[3] < std::get<3>(t).size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = std::get<0>(t)[i[0]]->index;
                    temp[1] = std::get<1>(t)[i[1]]->index;
                    temp[2] = std::get<2>(t)[i[2]]->index;
                    temp[3] = std::get<3>(t)[i[3]]->index;
                    auto [score, qcfs] = get_pvalue_and_qCFs(
                        input, temp, branch_id, qcf_writer
                    );
                    if (min < 0 || score < min) {
                        min = score;
                        minimizer[0] = temp[0]; minimizer[1] = temp[1]; minimizer[2] = temp[2]; minimizer[3] = temp[3];
                    }
                }
            }
        }
    }

    return min;
}


weight_t SpeciesTree::search_quard(std::vector<Tree *> &input, std::vector<std::vector<index_t>> &quad, index_t* minimizer) {
    index_t i[4] = {0, 0, 0, 0};
    weight_t min = -1;

    for (i[0] = 0; i[0] < quad[0].size(); i[0] ++) {
        for (i[1] = 0; i[1] < quad[1].size(); i[1] ++) {
            for (i[2] = 0; i[2] < quad[2].size(); i[2] ++) {
                for (i[3] = 0; i[3] < quad[3].size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = quad[0][i[0]];
                    temp[1] = quad[1][i[1]];
                    temp[2] = quad[2][i[2]];
                    temp[3] = quad[3][i[3]];
                    // std::cout << "debugging line 755 " << std::endl;
                    weight_t score = get_pvalue(input, temp);
                    // std::cout << "done with get_pvalue " << std::endl;
                    if (min < 0 || score < min) {
                        min = score;
                        minimizer[0] = temp[0]; minimizer[1] = temp[1]; minimizer[2] = temp[2]; minimizer[3] = temp[3];
                    }
                }
            }
        }
    }

    return min;
}

/// for average qcfs

void SpeciesTree::postorder_nodes(Node* root, std::vector<Node*>& out) {
    out.clear();
    if (!root) return;

    std::vector<Node*> st1;
    std::vector<Node*> st2;

    st1.push_back(root);
    while (!st1.empty()) {
        Node* u = st1.back();
        st1.pop_back();
        st2.push_back(u);
        for (Node* c : u->children) {
            if (c) st1.push_back(c);
        }
    }
    while (!st2.empty()) {
        out.push_back(st2.back());
        st2.pop_back();
    }
}

static inline double I_of_C(const std::array<std::array<std::size_t,4>,3>& C) {
    const auto& A = C[0];
    const auto& B = C[1];
    const auto& R = C[2];
    return
        (double)A[0]*B[1]*R[2]*R[3] + (double)A[1]*B[0]*R[2]*R[3] +
        (double)A[2]*B[3]*R[0]*R[1] + (double)A[3]*B[2]*R[0]*R[1] +
        (double)R[0]*B[1]*A[2]*A[3] + (double)R[1]*B[0]*A[2]*A[3] +
        (double)R[2]*B[3]*A[0]*A[1] + (double)R[3]*B[2]*A[0]*A[1] +
        (double)A[0]*R[1]*B[2]*B[3] + (double)A[1]*R[0]*B[2]*B[3] +
        (double)A[2]*R[3]*B[0]*B[1] + (double)A[3]*R[2]*B[0]*B[1];
}

// overload: buckets is std::array<vec,4>
static inline void build_bucket_of(std::vector<int8_t>& bucket_of,
                                   const std::array<std::vector<index_t>,4>& buckets,
                                   index_t num_taxa)
{
    bucket_of.assign((size_t)num_taxa, (int8_t)-1);
    for (int b = 0; b < 4; ++b) {
        for (index_t t : buckets[b]) {
            if ((size_t)t < bucket_of.size())
                bucket_of[(size_t)t] = (int8_t)b;
        }
    }
}

weight_t SpeciesTree::F_bucket_topology(const std::vector<Tree*>& input,
                         const std::array<std::vector<index_t>,4>& buckets,
                         index_t num_taxa,
                         const std::array<int,4>& roles)
{
    std::vector<int8_t> bucket_of;
    build_bucket_of(bucket_of, buckets, num_taxa);

    const std::size_t Sz[4] = {
        buckets[roles[0]].size(),
        buckets[roles[1]].size(),
        buckets[roles[2]].size(),
        buckets[roles[3]].size()
    };

    weight_t r = 0.0;

    std::vector<Node*> order;
    std::vector<std::array<std::size_t,4>> S;

    for (Tree* tree : input) {
        Node* root = tree->get_root();
        postorder_nodes(root, order);

        S.clear();

        for (Node* u : order) {
            if (u->children.size() == 0) {
                // leaf
                index_t t = u->index;
                std::array<std::size_t,4> c = {0,0,0,0};

                int8_t b = ((size_t)t < bucket_of.size()) ? bucket_of[(size_t)t] : (int8_t)-1;
                if (b == roles[0]) c[0] = 1;
                else if (b == roles[1]) c[1] = 1;
                else if (b == roles[2]) c[2] = 1;
                else if (b == roles[3]) c[3] = 1;

                S.push_back(c);
            } else {
                // binary assumed

                auto C2 = S.back(); S.pop_back();
                auto C1 = S.back(); S.pop_back();

                std::array<std::size_t,4> sum = {
                    C1[0] + C2[0],
                    C1[1] + C2[1],
                    C1[2] + C2[2],
                    C1[3] + C2[3]
                };

                std::array<std::size_t,4> rest = {
                    Sz[0] - sum[0],
                    Sz[1] - sum[1],
                    Sz[2] - sum[2],
                    Sz[3] - sum[3]
                };

                std::array<std::array<std::size_t,4>,3> C = {C1, C2, rest};
                r += I_of_C(C);

                S.push_back(sum);
            }
        }
    }

    if (input.empty()) return 0.0;

    return (r * 0.5) / (weight_t)input.size();
}

std::array<weight_t,3> SpeciesTree::freq_three_toplogies(const std::vector<Tree*>& input, Node* blob_node,
                             const std::array<index_t,4>& bucket_ids,
                                       index_t num_taxa)
{   
    index_t temp[4];
    for(index_t i = 0; i < 4; i++) {
        temp[i] = bucket_ids[i];
    }

    quartet_t q = join(temp);
    if (qCFs_average_cache.find(q) != qCFs_average_cache.end()) {
        return qCFs_average_cache[q];
    }

    
    const std::array<int,4> Q1 = {0,1,2,3}; // (0,1)|(2,3)
    const std::array<int,4> Q2 = {0,2,1,3}; // (0,2)|(1,3)
    const std::array<int,4> Q3 = {0,3,1,2}; // (0,3)|(1,2)
    
    std::array<std::vector<index_t>,4> buckets = {blob_node->multi_partitions[bucket_ids[0]],
                                                   blob_node->multi_partitions[bucket_ids[1]],
                                                   blob_node->multi_partitions[bucket_ids[2]],
                                                   blob_node->multi_partitions[bucket_ids[3]]};

    weight_t f1 = F_bucket_topology(input, buckets, num_taxa, Q1);
    weight_t f2 = F_bucket_topology(input, buckets, num_taxa, Q2);
    weight_t f3 = F_bucket_topology(input, buckets, num_taxa, Q3);

    weight_t m = f1 + f2 + f3;
    if (m <= 0.0) return {0.0, 0.0, 0.0};
    qCFs_average_cache[q] = {f1/m, f2/m, f3/m};
    return {f1/m, f2/m, f3/m};
}



// search quard with number of iterations 


weight_t SpeciesTree::search_quard_heuristic(std::vector<Tree *> &input,
                                            std::vector<std::vector<index_t>> &quad,
                                            unsigned long int iter_limit,
                                            index_t *minimizer,
                                            index_t branch_id,
                                            QCFWriter *qcf_writer) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;

    // basic guard
    if (quad.size() != 4) return (weight_t)-1;
    for (int d = 0; d < 4; ++d) if (quad[d].empty()) return (weight_t)-1;

    while (count < iter_limit) {
        // random start: one from each bucket
        indices[0] = quad[0][rand() % quad[0].size()];
        indices[1] = quad[1][rand() % quad[1].size()];
        indices[2] = quad[2][rand() % quad[2].size()];
        indices[3] = quad[3][rand() % quad[3].size()];

        weight_t old_min = min;

        count += neighbor_search_quard(
            input, quad, indices, &min, branch_id, qcf_writer
        );

        if (old_min < 0 || min < old_min) {
            minimizer[0] = indices[0];
            minimizer[1] = indices[1];
            minimizer[2] = indices[2];
            minimizer[3] = indices[3];
        }
    }

    return min;
}

weight_t SpeciesTree::search_quard_heuristic(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table,
                                            std::vector<std::vector<index_t>> &quad,
                                            unsigned long int iter_limit,
                                            index_t *minimizer) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;

    // basic guard
    if (quad.size() != 4) return (weight_t)-1;
    for (int d = 0; d < 4; ++d) if (quad[d].empty()) return (weight_t)-1;

    while (count < iter_limit) {
        // random start: one from each bucket
        indices[0] = quad[0][rand() % quad[0].size()];
        indices[1] = quad[1][rand() % quad[1].size()];
        indices[2] = quad[2][rand() % quad[2].size()];
        indices[3] = quad[3][rand() % quad[3].size()];

        weight_t old_min = min;

        count += neighbor_search_quard(qCFs_table, quad, indices, &min);

        if (old_min < 0 || min < old_min) {
            minimizer[0] = indices[0];
            minimizer[1] = indices[1];
            minimizer[2] = indices[2];
            minimizer[3] = indices[3];
        }
    }

    return min;
}



size_t SpeciesTree::neighbor_search_quard(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table,
                                         std::vector<std::vector<index_t>> &quad,
                                         index_t *current,
                                         weight_t *min) {
    auto [current_score, _] = get_pvalue_and_qCFs(qCFs_table,current);
    size_t k = 1;

    while (true) {
        index_t best_next[4] = { current[0], current[1], current[2], current[3] };
        weight_t best_next_score = current_score;

        // Try changing one coordinate at a time (within its bucket)
        for (int d = 0; d < 4; ++d) {
            index_t temp[4] = { current[0], current[1], current[2], current[3] };

            for (size_t i = 0; i < quad[d].size(); ++i) {
                index_t new_index = quad[d][i];
                if (new_index == current[d]) continue;

                temp[d] = new_index;
                auto [temp_score, _] = get_pvalue_and_qCFs(qCFs_table, temp);
                ++k;

                if (temp_score >= 0 && temp_score < best_next_score) {
                    best_next_score = temp_score;
                    best_next[0] = temp[0];
                    best_next[1] = temp[1];
                    best_next[2] = temp[2];
                    best_next[3] = temp[3];
                }

                temp[d] = current[d]; 
            }
        }

        // no improving neighbor
        if (best_next_score >= current_score) break;

        // take best move
        current_score = best_next_score;
        current[0] = best_next[0];
        current[1] = best_next[1];
        current[2] = best_next[2];
        current[3] = best_next[3];
    }

    if (*min < 0 || current_score < *min) *min = current_score;
    return k;
}



size_t SpeciesTree::neighbor_search_quard(std::vector<Tree *> &input,
                                         std::vector<std::vector<index_t>> &quad,
                                         index_t *current,
                                         weight_t *min,
                                         index_t branch_id,
                                         QCFWriter *qcf_writer) {
    auto [current_score, current_qcfs] = get_pvalue_and_qCFs(
        input, current, branch_id, qcf_writer
    );
    size_t k = 1;

    while (true) {
        index_t best_next[4] = { current[0], current[1], current[2], current[3] };
        weight_t best_next_score = current_score;

        // Try changing one coordinate at a time (within its bucket)
        for (int d = 0; d < 4; ++d) {
            index_t temp[4] = { current[0], current[1], current[2], current[3] };

            for (size_t i = 0; i < quad[d].size(); ++i) {
                index_t new_index = quad[d][i];
                if (new_index == current[d]) continue;

                temp[d] = new_index;
                auto [temp_score, temp_qcfs] = get_pvalue_and_qCFs(
                    input, temp, branch_id, qcf_writer
                );
                ++k;

                if (temp_score >= 0 && temp_score < best_next_score) {
                    best_next_score = temp_score;
                    best_next[0] = temp[0];
                    best_next[1] = temp[1];
                    best_next[2] = temp[2];
                    best_next[3] = temp[3];
                }

                temp[d] = current[d]; 
            }
        }

        // no improving neighbor
        if (best_next_score >= current_score) break;

        // take best move
        current_score = best_next_score;
        current[0] = best_next[0];
        current[1] = best_next[1];
        current[2] = best_next[2];
        current[3] = best_next[3];
    }

    if (*min < 0 || current_score < *min) *min = current_score;
    return k;
}





weight_t SpeciesTree::search_3f1a(std::vector<Tree *> &input, std::tuple<std::vector<Node *>, std::vector<Node *>, std::vector<Node *>, std::vector<Node *>> *quad, index_t* minimizer, index_t branch_id, QCFWriter *qcf_writer) {
    index_t i[4] = {0, 0, 0, 0};
    weight_t min = -1;
    
    auto& t = *quad;

    for (index_t alter = 0; alter < 4; alter ++) {
        
        const std::vector<Node*>* current_vec = nullptr;
        
        const std::vector<size_t> partitions_size = {std::get<0>(t).size(), std::get<1>(t).size(), std::get<2>(t).size(), std::get<3>(t).size()};
        

        switch (alter) {
            case 0: current_vec = &std::get<0>(t); break;
            case 1: current_vec = &std::get<1>(t); break;
            case 2: current_vec = &std::get<2>(t); break;
            case 3: current_vec = &std::get<3>(t); break;
        }
        
        index_t temp[4];
        
        for (int j = 0; j < 4; j++) {temp[j] = i[j];}

        for (temp[alter] = 0; temp[alter] < current_vec->size(); temp[alter] ++) {
            index_t cur_quart[4];
            cur_quart[0] = std::get<0>(t)[temp[0]]->index;
            cur_quart[1] = std::get<1>(t)[temp[1]]->index;
            cur_quart[2] = std::get<2>(t)[temp[2]]->index;
            cur_quart[3] = std::get<3>(t)[temp[3]]->index;
            // if it is not the altered partition, we use random index from the partition
            
            // for (int j = 0; j < 4; j++) {
            //     if (j == alter) {
            //         switch (j) {
            //             case 0: cur_quart[0] = std::get<0>(t)[temp[0]]->index; break;
            //             case 1: cur_quart[1] = std::get<1>(t)[temp[1]]->index; break;
            //             case 2: cur_quart[2] = std::get<2>(t)[temp[2]]->index; break;
            //             case 3: cur_quart[3] = std::get<3>(t)[temp[3]]->index; break;
            //         }
            //     } else {
            //         index_t rand_index = rand() % partitions_size[j];
            //         switch (j) {
            //             case 0: cur_quart[0] = std::get<0>(t)[rand_index]->index; break;
            //             case 1: cur_quart[1] = std::get<1>(t)[rand_index]->index; break;
            //             case 2: cur_quart[2] = std::get<2>(t)[rand_index]->index; break;
            //             case 3: cur_quart[3] = std::get<3>(t)[rand_index]->index; break;
            //         }
            //     }
            // }

            auto [score, qcfs] = get_pvalue_and_qCFs(
                input, cur_quart, branch_id, qcf_writer
            );
            
            
            if (min < 0 || score < min) {
                i[alter] = temp[alter];
                minimizer[0] = cur_quart[0];
                minimizer[1] = cur_quart[1];
                minimizer[2] = cur_quart[2];
                minimizer[3] = cur_quart[3];
                min = score;
            }
        }
    }
    return min;
}



bool SpeciesTree::is_match_with_split(const std::array<weight_t,3>& qcfs, index_t node_a1_id, index_t node_a2_id, index_t *indices) {
    // determine the split of the quartet based on qcf
    index_t temp[4];
    for (int i = 0; i < 4; ++i) temp[i] = indices[i];
    std::sort(temp, temp + 4);
    std::array<std::array<index_t, 4>, 2> top_2_frequent_topologies = computed_displayed_quartet_toplogy(temp, qcfs);
    std::array<index_t, 4> most_frequent_toplogy = top_2_frequent_topologies[0];
    if (most_frequent_toplogy[0] == node_a1_id && most_frequent_toplogy[1] == node_a2_id) {
        return true;
    } else if (most_frequent_toplogy[0] == node_a2_id && most_frequent_toplogy[1] == node_a1_id) {
        return true;
    } else if (most_frequent_toplogy[2] == node_a1_id && most_frequent_toplogy[3] == node_a2_id) {
        return true;
    } else if (most_frequent_toplogy[2] == node_a2_id && most_frequent_toplogy[3] == node_a1_id) {
        return true;
    } else {
        return false;
    }
}


weight_t SpeciesTree::search_2f2a(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t* minimizer, size_t &split_match_count, size_t &split_mismatch_count, index_t branch_id, QCFWriter *qcf_writer) {
    index_t i[4];
    weight_t min = -1;
    i[0] = 0; i[2] = 0;
        for (i[1] = 1; i[1] < A.size(); i[1] ++) {
                for (i[3] = 1; i[3] < B.size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = A[i[0]]->index;
                    temp[1] = A[i[1]]->index;
                    temp[2] = B[i[2]]->index;
                    temp[3] = B[i[3]]->index;
                    auto [score, qcfs] = get_pvalue_and_qCFs(
                        input, temp, branch_id, qcf_writer
                    );
                    if (is_match_with_split(qcfs, A[i[0]]->index, A[i[1]]->index, temp)) {
                        split_match_count++;
                    } else {
                        split_mismatch_count++;
                    }
                    
                    if (min < 0 || score < min) {
                        min = score;
                        minimizer[0] = temp[0]; minimizer[1] = temp[1]; minimizer[2] = temp[2]; minimizer[3] = temp[3];
                    }
                }
        }
    return min;
}

weight_t SpeciesTree::search_star(std::vector<Tree *> &input, 
                                  std::vector<Node *> &A,
                                  std::vector<Node *> &B) {
    index_t i[4];
    weight_t min = -1;
    size_t count = 0;
    for (i[0] = 0; i[0] < A.size(); i[0] ++) {
        for (i[1] = i[0] + 1; i[1] < A.size(); i[1] ++) {
            for (i[2] = 0; i[2] < B.size(); i[2] ++) {
                for (i[3] = i[2] + 1; i[3] < B.size(); i[3] ++) {
                    index_t temp[4];
                    temp[0] = A[i[0]]->index;
                    temp[1] = A[i[1]]->index;
                    temp[2] = B[i[2]]->index;
                    temp[3] = B[i[3]]->index;
                    count ++;
                    weight_t f[3];
                    weight_t score = get_pvalue_star(input, temp);
                    if (min < 0 || score > min) {
                        min = score;
                    }
                }
            }
        }
    }
    //std::cout << "naive iter: " << count << std::endl;
    return min;
}

//// fixed the bug in minimizer updates in the heuristic search
weight_t SpeciesTree::search(std::vector<Tree *> &input, 
                             std::vector<Node *> &A,
                             std::vector<Node *> &B,
                             size_t iter_limit,
                             index_t* minimizer,
                             index_t branch_id,
                             QCFWriter *qcf_writer) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;

    while (count < iter_limit) {
        indices[0] = A[rand() % A.size()]->index;
        do {indices[1] = A[rand() % A.size()]->index;} while (indices[0] == indices[1]);
        indices[2] = B[rand() % B.size()]->index;
        do {indices[3] = B[rand() % B.size()]->index;} while (indices[2] == indices[3]);
        
        weight_t old_min = min;

        count += neighbor_search(
            input, A, B, indices, &min, branch_id, qcf_writer
        );

        if (min < old_min || old_min < 0) {
            minimizer[0] = indices[0]; minimizer[1] = indices[1]; minimizer[2] = indices[2]; minimizer[3] = indices[3];
        }

        // std::cout << i << ' ' << min << std::endl;
    }
    
    //std::cout << "heuristic iter: " << count << std::endl;
    return min;
}

weight_t SpeciesTree::neighbor_search(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *current, weight_t *min, index_t branch_id, QCFWriter *qcf_writer) {
    weight_t current_f[3];
    weight_t current_score = get_pvalue_and_qCFs(
        input, current, branch_id, qcf_writer
    ).first;
    size_t k = 1;
    while (true) {
        index_t temp[4], next[4];
        weight_t temp_score, next_score = -1, temp_f[3], next_f[3];
        for (index_t j = 0; j < 4; j ++) 
            temp[j] = current[j];
        for (index_t i = 0; i < A.size(); i ++) {
            index_t new_index = A[i]->index;
            if (new_index == current[0] || new_index == current[1]) continue;
            temp[0] = new_index;
            temp_score = get_pvalue_and_qCFs(
                input, temp, branch_id, qcf_writer
            ).first; k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[0] = current[0];
            temp[1] = new_index;
            temp_score = get_pvalue_and_qCFs(
                input, temp, branch_id, qcf_writer
            ).first; k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[1] = current[1];
        }
        for (index_t i = 0; i < B.size(); i ++) {
            index_t new_index = B[i]->index;
            if (new_index == current[2] || new_index == current[3]) continue;
            temp[2] = new_index;
            temp_score = get_pvalue_and_qCFs(
                input, temp, branch_id, qcf_writer
            ).first; k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[2] = current[2];
            temp[3] = new_index;
            temp_score = get_pvalue_and_qCFs(
                input, temp, branch_id, qcf_writer
            ).first; k ++;
            if (next_score < 0 || temp_score < next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[3] = current[3];
        }
        if (next_score >= current_score) break;
        current_score = next_score;
        for (index_t j = 0; j < 4; j ++) 
            current[j] = next[j];
        current_f[0] = next_f[0]; current_f[1] = next_f[1]; current_f[2] = next_f[2];
    }
    if (*min < 0 || *min > current_score) {
        *min = current_score;
    }
    return k;
}

weight_t SpeciesTree::search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, size_t iter_limit) {
    index_t indices[4];
    weight_t min = -1;
    size_t count = 0;
    while (count < iter_limit) {
        indices[0] = A[rand() % A.size()]->index;
        do {indices[1] = A[rand() % A.size()]->index;} while (indices[0] == indices[1]);
        indices[2] = B[rand() % B.size()]->index;
        do {indices[3] = B[rand() % B.size()]->index;} while (indices[2] == indices[3]);
        count += neighbor_search_star(input, A, B, indices, &min);
        // std::cout << i << ' ' << min << std::endl;
    }
    //std::cout << "heuristic iter: " << count << std::endl;
    return min;
}

weight_t SpeciesTree::neighbor_search_star(std::vector<Tree *> &input, std::vector<Node *> &A, std::vector<Node *> &B, index_t *current, weight_t *min) {
    weight_t current_f[3];
    weight_t current_score = get_pvalue_star(input, current);
    size_t k = 1;
    while (true) {
        index_t temp[4], next[4];
        weight_t temp_score, next_score = -1, temp_f[3], next_f[3];
        for (index_t j = 0; j < 4; j ++) 
            temp[j] = current[j];
        for (index_t i = 0; i < A.size(); i ++) {
            index_t new_index = A[i]->index;
            if (new_index == current[0] || new_index == current[1]) continue;
            temp[0] = new_index;
            temp_score = get_pvalue_star(input, temp); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[0] = current[0];
            temp[1] = new_index;
            temp_score = get_pvalue_star(input, temp); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[1] = current[1];
        }
        for (index_t i = 0; i < B.size(); i ++) {
            index_t new_index = B[i]->index;
            if (new_index == current[2] || new_index == current[3]) continue;
            temp[2] = new_index;
            temp_score = get_pvalue_star(input, temp); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[2] = current[2];
            temp[3] = new_index;
            temp_score = get_pvalue_star(input, temp); k ++;
            if (next_score < 0 || temp_score > next_score) {
                next_score = temp_score;
                for (index_t j = 0; j < 4; j ++) 
                    next[j] = temp[j];
                next_f[0] = temp_f[0]; next_f[1] = temp_f[1]; next_f[2] = temp_f[2];
            }
            temp[3] = current[3];
        }
        if (next_score <= current_score) break;
        current_score = next_score;
        for (index_t j = 0; j < 4; j ++) 
            current[j] = next[j];
        current_f[0] = next_f[0]; current_f[1] = next_f[1]; current_f[2] = next_f[2];
    }
    if (*min < 0 || *min < current_score) {
        *min = current_score;
    }
    return k;
}

// the get_quartet will return the index = the index of sibling of the 0 index in indices - 1 in the indices array, 
// here we alway gonna passing the hybrid to the index of 0 in the indices array, 
// thus we can easily identify the sibilibing of hybrid in the displayed toplogy 
std::pair<weight_t, std::array<weight_t, 3>> SpeciesTree::get_pvalue_and_qCFs(
        std::vector<Tree *> &input,
        index_t *indices,
        index_t branch_id,
        QCFWriter *qcf_writer) {
    
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);

    if (pvalues.find(q) == pvalues.end() || qCFs_cache.find(q) == qCFs_cache.end()) {
        weight_t qCF[3] = {0, 0, 0};
        for (Tree *t : input) {
            index_t topology = t->get_quartet(temp);
            if (topology >= 0) qCF[topology] += 1;
        }

        if ((qCF[0] + qCF[1] + qCF[2]) == 0)
            pvalues[q] = 1.0;
        else
            pvalues[q] = pvalue(qCF);
        qCFs_cache[q] = {qCF[0], qCF[1], qCF[2]}; // main the global qCFs map
    }

    const std::array<weight_t, 3> &qcfs = qCFs_cache.at(q);
    if (qcf_writer != nullptr) {
        qcf_writer->write(branch_id, temp, qcfs);
    }
    return {pvalues.at(q), qcfs};
}

std::pair<weight_t, std::array<weight_t, 3>> SpeciesTree::get_pvalue_and_qCFs(std::unordered_map<quartet_t, std::array<weight_t, 3>> &qCFs_table, index_t *indices) {
    
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);


    std::array<weight_t, 3> qCFs;
    qCFs = qCFs_table[q];
    weight_t tmp[3] = {qCFs[0], qCFs[1], qCFs[2]};
    pvalues[q] = pvalue(tmp);
    qCFs_cache[q] = qCFs;
    return {pvalues[q], qCFs};
}


std::array<std::array<index_t, 4>, 2> SpeciesTree::computed_displayed_quartet_toplogy(index_t *indices,
                                                const std::array<weight_t,3>& qcf) {
    index_t temp[4];
    for (int i = 0; i < 4; ++i) temp[i] = indices[i];
    std::sort(temp, temp + 4);

    int best1 = 0;
    for (int t = 1; t < 3; ++t) {
        if (qcf[t] > qcf[best1] || (qcf[t] == qcf[best1] && t < best1)) best1 = t;
    }

    int best2 = (best1 == 0 ? 1 : 0);
    for (int t = 0; t < 3; ++t) {
        if (t == best1) continue;
        if (qcf[t] > qcf[best2] || (qcf[t] == qcf[best2] && t < best2)) best2 = t;
    }

    auto topo_to_quartet = [&](int topo) -> std::array<index_t, 4> {
        // pairs: (ret[0],ret[1]) and (ret[2],ret[3])
        switch (topo) {
            case 0: return { temp[0], temp[1], temp[2], temp[3] }; // 01|23
            case 1: return { temp[0], temp[2], temp[1], temp[3] }; // 02|13
            case 2: return { temp[0], temp[3], temp[1], temp[2] }; // 03|12
            default: return { temp[0], temp[1], temp[2], temp[3] };
        }
    };

    return { topo_to_quartet(best1), topo_to_quartet(best2) };
}



std::array<std::array<index_t, 4>, 2> SpeciesTree::computed_displayed_quartet_toplogy(index_t *indices) {
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);

    index_t best1 = 0, best2 = 1;
    for (index_t t = 0; t < 3; ++t) {
        if (qCFs_cache[t] > qCFs_cache[best1] || (qCFs_cache[t] == qCFs_cache[best1] && t < best1)) best1 = t;
    }
    best2 = (best1 == 0 ? 1 : 0);

    for (index_t t = 0; t < 3; ++t) {
        if (t == best1) continue;
        if (qCFs_cache[t] > qCFs_cache[best2] || (qCFs_cache[t] == qCFs_cache[best2] && t < best2)) best2 = t;
    }

    auto topo_to_quartet = [&](int topo) -> std::array<index_t, 4> {
        // Returns {a,b,c,d} meaning (a,b) and (c,d) are the two sibling pairs.
        switch (topo) {
            case 0: return { temp[0], temp[1], temp[2], temp[3] }; // 01|23
            case 1: return { temp[0], temp[2], temp[1], temp[3] }; // 02|13
            case 2: return { temp[0], temp[3], temp[1], temp[2] }; // 03|12
            default: return { temp[0], temp[1], temp[2], temp[3] };
        }
    };

    return {topo_to_quartet(best1), topo_to_quartet(best2)};
}

static inline index_t sibling_in_topology(const std::array<index_t,4> &topo, index_t taxon) {
    // pair1: topo[0] <-> topo[1], pair2: topo[2] <-> topo[3]
    if (topo[0] == taxon) return topo[1];
    if (topo[1] == taxon) return topo[0];
    if (topo[2] == taxon) return topo[3];
    if (topo[3] == taxon) return topo[2];
    std::cout << "Error: taxon " << taxon << " not found in topology" << std::endl;
    return (index_t)-1;
}

std::array<index_t, 2> SpeciesTree::siblings_in_two_best_topologies(const std::array<std::array<index_t,4>,2> &best2,
                                             index_t taxon) {
    std::array<index_t,2> sibs;
    sibs[0] = sibling_in_topology(best2[0], taxon);
    sibs[1] = sibling_in_topology(best2[1], taxon);
    return sibs;
}



weight_t SpeciesTree::get_pvalue(std::vector<Tree *> &input, index_t *indices) {
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);

    if (pvalues.find(q) == pvalues.end()) {
        weight_t qCF[3] = {0, 0, 0};
        for (Tree *t : input) {
            index_t topology = t->get_quartet(temp);
            if (topology >= 0) qCF[topology] += 1;
        }

        /*
        std::sort(qCF, qCF + 3);
        std::vector<weight_t> all = pvalue_all(indices);
        std::sort(all.begin(), all.begin() + 3);
        bool verification = qCF[0] == all[0] && qCF[1] == all[1] && qCF[2] == all[2];
        if (! verification) {
            std::cout << qCF[0] << ' ' << qCF[1] << ' ' << qCF[2] << std::endl;
            std::cout << all[0] << ' ' << all[1] << ' ' << all[2] << std::endl;
        }
        assert(verification);
        */

        if ((qCF[0] + qCF[1] + qCF[2]) == 0)
            pvalues[q] = 1.0;
        else
            pvalues[q] = pvalue(qCF);
        qCFs_cache[q] = {qCF[0], qCF[1], qCF[2]}; // main the global qCFs map
    }
    return pvalues[q];
}

weight_t SpeciesTree::get_pvalue_star(std::vector<Tree *> &input, index_t *indices) {
    index_t temp[4];
    for (index_t i = 0; i < 4; i ++) 
        temp[i] = indices[i];
    std::sort(temp, temp + 4);
    quartet_t q = join(temp);

    if (pvalues_star.find(q) == pvalues_star.end()) {
        weight_t qCF[3] = {0, 0, 0};
        for (Tree *t : input) {
            index_t topology = t->get_quartet(temp);
            if (topology >= 0) qCF[topology] += 1;
        }
        
        /*
        std::sort(qCF, qCF + 3);
        std::vector<weight_t> all = pvalue_all(indices);
        std::sort(all.begin(), all.begin() + 3);
        bool verification = qCF[0] == all[0] && qCF[1] == all[1] && qCF[2] == all[2];
        if (! verification) {
            std::cout << qCF[0] << ' ' << qCF[1] << ' ' << qCF[2] << std::endl;
            std::cout << all[0] << ' ' << all[1] << ' ' << all[2] << std::endl;
        }
        assert(verification);
        */
        if ((qCF[0] + qCF[1] + qCF[2]) == 0)
            pvalues_star[q] = -1.0;  // set to -1 so wouldn't impact search for max
        else
            pvalues_star[q] = pvalue_star(qCF);
    }
    return pvalues_star[q];
}

index_t Tree::get_quartet(index_t *indices) {
    Node *leaves[4];
    for (index_t i = 0; i < 4; i ++) {
        // add for absent tadxon in gene tree
        if (!index2node[indices[i]]) {
            return -1;
        }
        leaves[i] = index2node[indices[i]];
    }
    Node *lowest = NULL;
    index_t nodes[4] = {-1, -1, -1, -1};
    index_t lowest_count = 0;
    for (index_t i = 0; i < 4; i ++) {
        for (index_t j = i + 1; j < 4; j ++) {
            // Node *a = LCA_naive(leaves[i], leaves[j]);
            Node *a = LCA_via_rmq(leaves[i], leaves[j]);
            if (lowest == NULL || a->depth > lowest->depth) {
                lowest = a;
                lowest_count = 0;
                nodes[0] = i;
                nodes[1] = j;
            }
            if (a == lowest) lowest_count ++;
        }
    }
    if (lowest_count != 1) return -1;
    if (nodes[0] == 0) 
        return nodes[1] - 1;
    else 
        return 5 - nodes[0] - nodes[1];
}

void Tree::LCA_preprocessing() {
    LCA_preprocessing_with_ett_rmq_sparse_table();
}

void Tree::LCA_depth_first_search(Node *root, std::vector<Node *> &stack) {
    root->depth = stack.size() - 1;
    for (index_t i = 1; i < stack.size(); i *= 2) 
        root->ancestors.push_back(stack[stack.size() - i - 1]);
    if (root->children.size() != 0) {
        for (Node *child : root->children) {
            stack.push_back(child);
            LCA_depth_first_search(child, stack);
            stack.pop_back();
        }
    }
}

Node *Tree::LCA_fast(Node *x, Node *y) {
    if (x->depth > y->depth) {
        Node *temp = x;
        x = y;
        y = temp;
    }
    for (index_t i = y->ancestors.size(); i >= 0; i --) {
        if (i >= y->ancestors.size()) continue;
        if (y->ancestors[i]->depth >= x->depth) 
            y = y->ancestors[i];
    }
    if (x == y) return y;
    for (index_t i = y->ancestors.size(); i >= 0; i --) {
        if (i >= y->ancestors.size()) continue;
        if (y->ancestors[i] != x->ancestors[i]) {
            y = y->ancestors[i];
            x = x->ancestors[i];
        }
    }
    return y->ancestors[0];
}

Node *Tree::LCA_naive(Node *a, Node *b) {
    while (a->depth > b->depth) a = a->parent;
    while (b->depth > a->depth) b = b->parent;
    while (a != b) {a = a->parent; b = b->parent;}
    return a;
}
#endif  // ENABLE_TOB
