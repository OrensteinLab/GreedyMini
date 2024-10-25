#pragma once
#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> 
#include <map>
#include <set>







void compute_and_store_densities(
    uint32_t w, uint32_t k, Config& config,
    std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map,
    std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map,
    std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map,
    std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map,
    std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map,
    std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map);

void temp_best_density(Config& config);

void expected_density_tests(Config& config, uint32_t max_w, uint32_t max_k);

void save_density_csvs(std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map);

void calculate_density_w_is_k_minus_one(const int max_time_seconds, Config& config, uint32_t max_k, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map);

void calculate_density_k_fixed_w_varying(const int max_time_seconds, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map);

void calculate_density_w_fixed_k_varying(const int max_time_seconds, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map);

void calculate_density_matrix(Config& config, uint32_t max_w, uint32_t max_k, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map);

void short_compute_and_store_densities(uint32_t w, uint32_t k, Config& config);

void big_w_extension(Config& config, uint32_t initial_w, uint32_t k, uint32_t max_w, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map);

void big_k_extension(Config& config, uint32_t w, uint32_t initial_k, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map);

void improve_from_previous_k(Config& config, uint32_t k, uint32_t w, double current_best_density, double max_swapping_time);

void particular_dna_tests(Config& config);

void dual_density_tests();

void save_particular_density_csvs(std::map<uint32_t, std::map<uint32_t, double>>& normal_order, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order, std::map<uint32_t, std::map<uint32_t, double>>& normal_order_on_expected, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order_on_expected);

void short_calculate_particular_density(uint32_t w, uint32_t k, Config& config);

void calculate_particular_density_k_fixed_w_varying(uint32_t k, const int max_time_in_seconds, Config& config, std::set<std::pair<uint32_t, uint32_t>>& computed_pairs, std::map<uint32_t, std::map<uint32_t, double>>& normal_order, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order, std::map<uint32_t, std::map<uint32_t, double>>& normal_order_on_expected, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order_on_expected);

void calculate_particular_density_w_fixed_k_varying(uint32_t w, const int max_time_in_seconds, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& normal_order, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order, std::map<uint32_t, std::map<uint32_t, double>>& normal_order_on_expected, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order_on_expected, std::set<std::pair<uint32_t, uint32_t>>& computed_pairs);

