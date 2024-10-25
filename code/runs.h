#ifndef RUNS_H
#define RUNS_H

#include <vector>
#include <cstdint>
#include "config.h"

// Structure to store results

struct Result {
    std::vector<uint64_t> ans;
    uint64_t gc_count;
    double actual_alpha; // since we randomize around the error
};

struct KeepTrackTop3 {
    std::vector<std::vector<uint64_t>> top_3_orders; // assume filled with 0 0 0
    std::vector<uint64_t> top_3_gc_counts;
    std::vector<uint64_t> all_gc_counts;
};

// Function to run ecogreed in a single thread multiple times
void run_ecogreed(uint32_t W, uint32_t K, Result& best_result, int iterations, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans, Config& config);

void run_ecogreed_specific(uint32_t W, uint32_t K, Result& best_result, int iterations, std::vector<uint64_t>& allSequences,
    std::vector<uint64_t> gc, std::vector<uint64_t> nongc, std::vector<uint64_t> ans,
    std::vector<std::vector<uint64_t>>& kmer_to_gc_string_id, std::vector<std::vector<uint64_t>>& kmer_to_non_gc_string_id, std::vector<std::vector<uint64_t>>& string_id_to_non_gc_kmers, std::vector<uint64_t>& string_id_to_gc_prefix,
    std::vector<uint64_t>& string_id_to_gc_suffix, Config& config);

// Function to run ecogreed in parallel
void parallel_run(Config &config, uint32_t W, uint32_t K, const bool save_results);
void parallel_run_specific(Config &config, uint32_t W, uint32_t K,  const bool save_results);
void single_run_swapper(uint32_t W, uint32_t K, double min_alpha, double max_alpha, const double max_time_seconds);


#endif // RUNS_H