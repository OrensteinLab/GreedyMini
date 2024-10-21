#ifndef RUNS_H
#define RUNS_H

#include <vector>
#include <cstdint>
#include "config.h"

// Structure to store results

struct Result {
    std::vector<uint64_t> ans;
    uint64_t gc_count;
    double actual_error; // since we randomize around the error
};

struct KeepTrackTop3 {
    std::vector<std::vector<uint64_t>> top_3_orders; // assume filled with 0 0 0
    std::vector<uint64_t> top_3_gc_counts;
    std::vector<uint64_t> all_gc_counts;
};

// Function to run ecogreed in a single thread multiple times
void run_ecogreed(uint32_t W, uint32_t K, double err, Result& best_result, int iterations, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans, Config& config);

//void run_ecogreed_top_3(uint32_t W, uint32_t K, double err, KeepTrackTop3& orders_tracker, int iterations, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans);


void run_ecogreed_specific(uint32_t W, uint32_t K, double err, Result& best_result, int iterations, const std::string& path, const std::string& name, std::vector<uint64_t>& allSequences,
    std::vector<uint64_t> gc, std::vector<uint64_t> nongc, std::vector<uint64_t> ans,
    std::vector<std::vector<uint64_t>>& kmer_to_gc_string_id, std::vector<std::vector<uint64_t>>& kmer_to_non_gc_string_id, std::vector<std::vector<uint64_t>>& string_id_to_non_gc_kmers, std::vector<uint64_t>& string_id_to_gc_prefix,
    std::vector<uint64_t>& string_id_to_gc_suffix, Config &config);

// Function to run ecogreed in parallel
void parallel_run(Config &config, uint32_t W, uint32_t K, double err, const int total_runs, const bool save_results);

//void parallel_run_top_3(uint32_t W, uint32_t K, double err, const int total_runs);

void parallel_run_specific(Config &config, uint32_t W, uint32_t K, double err, const int total_runs, const bool save_results, std::string path, std::string name);

void single_run(uint32_t W, uint32_t K, double err, const bool save_results);

void single_run_specific(uint32_t W, uint32_t K, double err, const bool save_results, std::string path, std::string name);


void single_run_swapper(uint32_t W, uint32_t K, double err, const double max_time_seconds);

// TODO: use this
//void multiple_runs_swapper(uint32_t W, uint32_t K, double err, const double max_time_seconds, const uint32_t n_repetitions);



#endif // RUNS_H