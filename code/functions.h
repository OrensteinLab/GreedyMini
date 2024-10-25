#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <unordered_set>
#include <string>
#include <tuple>
#include <cstdint>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>

// Function to save a vector to a file
void save_vector_to_file(const std::vector<uint64_t>& vec, const std::string& filename);

// Function to load a vector from a file
std::vector<uint64_t> load_vector_from_file(const std::string& filename);

// Function to check if a file exists
bool file_exists(const std::string& filename);

// Function to create initial lists and save them to files
void make_init_lists(uint32_t w, uint32_t k);

void ensure_init_lists(uint32_t w, uint32_t k);

// Function to load the initial lists from files
std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>> load_init_lists(uint32_t w, uint32_t k);

// Function to check if the initial lists are saved
bool check_init_lists(uint32_t w, uint32_t k);

// Recursive function to initialize gamechanger and non-gamechanger counts
void initdfs(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t lev, uint32_t w, uint32_t k,
    std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc);

// Function to initialize the gamechanger and non-gamechanger lists
std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>> init_lists(uint32_t w, uint32_t k);

// Function to choose a random k-mer with the best nongc/gc ratio up to a multiplicative error
uint64_t randcand(uint32_t k, double err, const std::vector<uint64_t>& gc, const std::vector<uint64_t>& nongc, uint64_t rank);

uint64_t randcand_v2(uint32_t k, const std::vector<uint64_t>& gc, const std::vector<uint64_t>& nongc);

uint64_t randcand_v3(uint32_t k, const std::vector<uint64_t>& gc, const std::vector<uint64_t>& nongc);

// Function to compute the next block
std::pair<uint64_t, uint32_t> nextblock(uint64_t num, uint32_t w, uint32_t k, uint64_t kmer, uint64_t n_kmers, uint64_t pos);

// Recursive depth-first search for gamechangers
uint64_t gc_dfs(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev);

// Helper function for gc_dfs when the previous preference could be minimum
uint64_t gc_dfs_prev_could_be_min(uint32_t w, uint64_t n_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev);

// Helper function for gc_dfs when the previous preference is not minimum
uint64_t gc_dfs_pref_is_not_min(uint32_t w, uint64_t n_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev);


// Wrapper for the rightdfs function
uint64_t rightdfs_wrapper(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t lev, uint32_t w, uint32_t k, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc);

// Right depth-first search for level zero
uint64_t rightdfs_lev_zero(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t w, uint64_t n_kmers_mask, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc);

// Right depth-first search for non-zero levels
uint64_t rightdfs(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t lev_left, uint64_t n_kmers_mask, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc);


uint64_t rightdfs_wrapper_v2(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t lev, uint32_t w, uint32_t k, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans, uint64_t current_rank);
uint64_t rightdfs_lev_zero_v2(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t w, uint64_t n_kmers_mask, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans, uint64_t current_rank);
uint64_t rightdfs_v2(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t lev_left, uint64_t n_kmers_mask, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans, uint64_t current_rank);



std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>,
    std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>, std::vector<uint64_t>, std::vector<uint64_t>>
    init_lists_specific(uint32_t w, uint32_t k, std::string path);

std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>,
    std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>, std::vector<uint64_t>, std::vector<uint64_t>>
    load_init_lists_specific(uint32_t w, uint32_t k, std::string name);

bool check_init_lists_specific(uint32_t w, uint32_t k, std::string name);

void ensure_init_lists_specific(uint32_t w, uint32_t k, std::string path, std::string name);

void make_init_lists_specific(uint32_t w, uint32_t k, std::string path, std::string name);

void update_specific_vectors(std::vector<uint64_t>& gc,
    std::vector<uint64_t>& nongc,
    std::vector<std::vector<uint64_t>>& kmer_to_gc_string_id,
    std::vector<std::vector<uint64_t>>& kmer_to_non_gc_string_id,
    std::vector<std::vector<uint64_t>>& string_id_to_non_gc_kmers,
    std::vector<uint64_t>& string_id_to_gc_prefix,
    std::vector<uint64_t>& string_id_to_gc_suffix,
    boost::dynamic_bitset<>& dead_windows,
    uint64_t nextbest);


std::pair<std::vector<uint64_t>, uint64_t> ecogreed_specific(uint32_t w, uint32_t k, double err, std::string path, std::string name, std::vector<uint64_t>& allSequences,
    std::vector<uint64_t> gc, std::vector<uint64_t> nongc, std::vector<uint64_t> ans,
    std::vector<std::vector<uint64_t>>kmer_to_gc_string_id, std::vector<std::vector<uint64_t>>kmer_to_non_gc_string_id, std::vector<std::vector<uint64_t>>& string_id_to_non_gc_kmers, std::vector<uint64_t>& string_id_to_gc_prefix,
    std::vector<uint64_t>& string_id_to_gc_suffix);

void update_counters(uint32_t w, uint32_t k, uint64_t nextbest, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans, uint64_t rank);

// Greedy choice of a good 2-ary order in small space
std::pair<std::vector<uint64_t>, uint64_t> ecogreed(uint32_t w, uint32_t k, double err, std::vector<uint64_t> gc, std::vector<uint64_t> nongc, std::vector<uint64_t> ans);


// Function to compute density of order for a specific w
uint64_t specific_wide_gc(uint32_t w, uint32_t k, const std::vector<uint64_t>& order);

// Function to compute the density of a specific order
void wide_gc_report(uint32_t w, uint32_t k, double err);

std::vector<uint64_t> get_random_order(uint32_t k);

std::vector<uint64_t> load_all_sequences_particular(uint32_t w, uint32_t k, std::string path);




#endif // FUNCTIONS_H
