#pragma once
#ifndef DENSITY_AND_GC_H
#define DENSITY_AND_GC_H

#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
#include <string>
#include <tuple>
#include "window_types_manager.h"


using namespace boost::multiprecision;

uint64_t gc_iterative_particular_binary(uint32_t w, uint32_t k, std::vector<uint64_t> order, std::vector<uint64_t> numbers);

double density_particular_dna(uint32_t W, uint32_t K, std::string path, std::vector<uint64_t> order);

uint64_t prob_gc(uint32_t w, uint32_t k, const std::vector<uint64_t>& order);

uint64_t gc_dfs_extend_k(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev, uint32_t extended_k);

uint64_t gc_dfs_extend_k_prev_could_be_min(uint32_t w, uint32_t k, uint64_t n_extended_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev, uint32_t extended_k);

uint64_t gc_dfs_extend_k_pref_is_not_min(uint32_t w, uint32_t k, uint64_t n_extended_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev, uint32_t extended_k);

double density_expected_binary_extend_k(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint32_t extended_k);

std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>> get_gc_ruiners(uint32_t w, uint32_t k);
double get_special_case_expected_gc(uint32_t w, uint32_t k, WindowTypesManager window_types_manager);

uint64_t gc_dfs_dna(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t min_order, uint64_t min_indices, uint32_t lev, WindowTypesManager& window_types_manager, bool is_first_min, bool is_middle_min);

double density_expected_binary(uint32_t w, uint32_t k, const std::vector<uint64_t>& order);

double density_expected_dna(uint32_t w, uint32_t k, const std::vector<uint64_t>& order);

std::unordered_map<uint32_t, boost::multiprecision::cpp_int> dna_wide_gc(uint32_t maxw, uint32_t k, const std::vector<uint64_t>& order, bool verbose);

std::unordered_map<uint32_t, double> density_expected_dna_wide(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, bool verbose);

uint64_t gc_iterative_particular_dna(uint32_t w, uint32_t k, std::vector<uint64_t> order, std::vector<uint64_t> oddNumbers, std::vector<uint64_t> evenNumbers);

std::unordered_map<uint32_t, boost::multiprecision::cpp_int> wide_gc(uint32_t maxw, uint32_t k, const std::vector<uint64_t>& order, bool verbose);

std::unordered_map<uint32_t, double> density_expected_binary_wide(uint32_t maxw, uint32_t k, const std::vector<uint64_t>& order, bool verbose);

#endif // DENSITY_AND_GC_H

