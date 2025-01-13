#include <iostream>
#include <chrono>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <limits>
#include "functions.h" 
#include "runs.h"  
#include "tools.h"
#include "density_and_gc.h"
#include <algorithm> 
#include <cmath>
#include <random>      // For std::random_device and std::mt19937
#include <map>
#include <iomanip> // 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "tests.h"
#include <string>
#include "config.h"
#include <set>



void temp_best_density(Config& config) {
    std::map<uint32_t, std::map<uint32_t, double>> binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map;
    uint32_t w = 5;
    uint32_t k = 6;
    compute_and_store_densities(w, k, config, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);

    //w = 3;
    //k = 10;
    //compute_and_store_densities(w, k, config, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);


}

void expected_density_tests(Config& config, uint32_t max_w, uint32_t max_k) {
    print_to_both(config, "Starting generate orders and get densities\n");

    std::map<uint32_t, std::map<uint32_t, double>> binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map;
    calculate_density_matrix(config, max_w, max_k, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);

    const int max_time_seconds = config.max_mins_per_step * 60;
    calculate_density_w_fixed_k_varying(max_time_seconds, config, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);
    calculate_density_k_fixed_w_varying(max_time_seconds, config, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);
    //calculate_density_w_is_k_minus_one(max_time_seconds, config, max_k, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);


    save_density_csvs(binary_density_map, config, binary_density_after_swap_map, dna_density_map, dna_density_after_swap_map, binary_density_best_map, dna_density_best_map);

}

void particular_dna_tests(Config& config) {
    print_to_both(config, "\nRunning particular DNA tests \n");

    std::map<uint32_t, std::map<uint32_t, double>> normal_order, specifically_trained_order, normal_order_on_expected, specifically_trained_order_on_expected;

    const int max_time_in_seconds = config.max_mins_per_step * 60;
    std::set<std::pair<uint32_t, uint32_t>> computed_pairs; // To track computed cells

    uint32_t fixed_w = 12;
    calculate_particular_density_w_fixed_k_varying(fixed_w, max_time_in_seconds, config, normal_order, specifically_trained_order, normal_order_on_expected, specifically_trained_order_on_expected, computed_pairs);
    uint32_t fixed_k = 8;
    calculate_particular_density_k_fixed_w_varying(fixed_k, max_time_in_seconds, config, computed_pairs, normal_order, specifically_trained_order, normal_order_on_expected, specifically_trained_order_on_expected);

    save_particular_density_csvs(normal_order, config, specifically_trained_order, normal_order_on_expected, specifically_trained_order_on_expected);
}



void dual_density_tests() {
    uint32_t w = 12;
    uint32_t k = 8;
    uint32_t n_runs = 1000;
    std::vector<double> density_differences;
    std::vector<double> densities;
    std::vector<double> density_diff_to_density_ratio;

    for (uint32_t i = 0; i < n_runs; i++) {
        std::vector<uint64_t> order = get_random_order(k);
        std::vector<uint64_t> order_reversed = get_reversed_order(order);
        double binary_density = density_expected_binary(w, k, order);
        double binary_density_reversed = density_expected_binary(w, k, order_reversed);
        density_differences.push_back(binary_density - binary_density_reversed);
        densities.push_back(binary_density);
        densities.push_back(binary_density_reversed);
        double difference_abs = std::abs(binary_density - binary_density_reversed);
        density_diff_to_density_ratio.push_back(difference_abs / binary_density);
  //      std::cout << "normal density " << binary_density << "\t " << "reversed density " << binary_density_reversed << std::endl;

  //      // print the order and the reverse order
  //      std::cout << "Order: ";
  //      for (uint64_t val : order) {
		//	std::cout << val << " ";
		//}
  //      std::cout << std::endl;

		//std::cout << "Reverse order: ";
		//for (uint64_t val : order_reversed) {
		//	std::cout << val << " ";
		//}
		//std::cout << std::endl;

	}

    // make density differences positive
    for (double& diff : density_differences) {
        if (diff < 0) {
            diff = -diff;
        }
    }

    // calculate mean and standard deviation
    double mean = 0;
    for (double diff : density_differences) {
		mean += diff;
    }
	mean /= n_runs;

	double variance = 0;
	for (double diff : density_differences) {
        variance += (diff - mean) * (diff - mean);
	}
    variance /= n_runs;

    double std_dev = sqrt(variance);

    std::cout << "Mean: " << mean << std::endl;
    std::cout << "Standard deviation: " << std_dev << std::endl;


    // calculate max density difference
    double max_diff = 0;
    for (double diff : density_differences) {
        if (diff > max_diff) {
            max_diff = diff;
        }
    }

    std::cout << "Max density difference: " << max_diff << std::endl;


    // calculate mean and standard deviation for densities
    double mean_density = 0;
    for (double density : densities) {
		mean_density += density;
	}
    mean_density /= n_runs;

    double variance_density = 0;
    for (double density : densities) {
        variance_density += (density - mean_density) * (density - mean_density);
    }
    variance_density /= n_runs;

    double std_dev_density = sqrt(variance_density);

    std::cout << "Mean density: " << mean_density << std::endl;

    std::cout << "Standard deviation density: " << std_dev_density << std::endl;


    // calcualte the average ratio of the difference to the density
    double mean_diff_to_density_ratio = 0;
    for (double ratio : density_diff_to_density_ratio) {
		mean_diff_to_density_ratio += ratio;
	}
    // times 100 for percentage
    mean_diff_to_density_ratio *= 100;
    mean_diff_to_density_ratio /= n_runs;

    std::cout << "Mean difference to density ratio: " << mean_diff_to_density_ratio << "%" << std::endl;

}







void save_density_csvs(std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map)
{
    save_2d_to_csv(binary_density_map, "original_binary_density", "W", "K", config);
    save_2d_to_csv(binary_density_after_swap_map, "original_binary_density_after_swap", "W", "K", config);
    save_2d_to_csv(dna_density_map, "original_dna_density", "W", "K", config);
    save_2d_to_csv(dna_density_after_swap_map, "original_dna_density_after_swap", "W", "K", config);
    save_2d_to_csv(binary_density_best_map, "binary_density_best", "W", "K", config);
    save_2d_to_csv(dna_density_best_map, "dna_density_best", "W", "K", config);
}

void calculate_density_w_is_k_minus_one(const int max_time_seconds, Config& config, uint32_t max_k, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map)
{
    // Compute for w = k - 1 starting from w = 3
    auto start_time = std::chrono::steady_clock::now();
    for (uint32_t w = 3; w <= 40; ++w) {
        uint32_t k = w + 1;
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - start_time).count();
        if (elapsed_seconds >= max_time_seconds) {
            print_to_both(config, "Reached time limit for k = " + std::to_string(k) + ", stopping further iterations.\n");
            break;
        }
        if (k <= max_k) {
            if (binary_density_map.count(w) == 0 || binary_density_map[w].count(k) == 0) {
                compute_and_store_densities(w, k, config, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);
            }
        }
    }
}

void calculate_density_k_fixed_w_varying(const int max_time_seconds, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map)
{
    // Calculate for k fixed and w varying
    std::vector<uint32_t> k_values = { 8, 12 };

    for (uint32_t k : k_values) {
        auto start_time = std::chrono::steady_clock::now();

        uint32_t last_w_calculated = 3;

        // calculate until time limit
        for (uint32_t w = 3; w <= 40; ++w) {
            auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start_time).count();
            if (elapsed_seconds >= max_time_seconds) {
                print_to_both(config, "Reached time limit for k = " + std::to_string(k) + ", stopping further iterations.\n");
                break;
            }
            if (binary_density_map.count(w) == 0 || binary_density_map[w].count(k) == 0) {
                compute_and_store_densities(w, k, config, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);
            }
            last_w_calculated = w;
        }
        print_to_both(config, "Extending k from: W = " + std::to_string(last_w_calculated) + " K=" + std::to_string(k) + "\n");
        // afterwards, extend w until 200
        big_w_extension(config, last_w_calculated, k, 200, dna_density_best_map);
    }
}

void calculate_density_w_fixed_k_varying(const int max_time_seconds, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map)
{
    // Calculate for w fixed and k varying
    std::vector<uint32_t> w_values = { 8, 12 };
    for (uint32_t w : w_values) {
        auto start_time = std::chrono::steady_clock::now();

        uint32_t last_k_calculated = 3;

        // calculate until time limit
        for (uint32_t k = 3; k <= 40; ++k) {
            auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start_time).count();
            if (elapsed_seconds >= max_time_seconds) {
                print_to_both(config, "Reached time limit for k = " + std::to_string(k) + ", stopping further iterations.\n");
                break;
            }
            if (binary_density_map.count(w) == 0 || binary_density_map[w].count(k) == 0) {
                compute_and_store_densities(w, k, config, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);
            }
            last_k_calculated = k;
        }
        print_to_both(config, "Extending w from: W = " + std::to_string(w) + " K=" + std::to_string(last_k_calculated) + "\n");
        // afterwards, extend k until time limit
        big_k_extension(config, w, last_k_calculated, dna_density_best_map);
    }
}

void calculate_density_matrix(Config& config, uint32_t max_w, uint32_t max_k, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map)
{
    // First, compute for all combinations where 3 <= w <= max_w and 3 <= k <= max_k
    for (uint32_t w = 3; w <= max_w; ++w) {
        for (uint32_t k = 3; k <= max_k; ++k) {
            if (binary_density_map.count(w) == 0 || binary_density_map[w].count(k) == 0) {
                compute_and_store_densities(w, k, config, binary_density_map, binary_density_after_swap_map, binary_density_best_map, dna_density_map, dna_density_after_swap_map, dna_density_best_map);
            }
        }
    }
}


void short_compute_and_store_densities(uint32_t w, uint32_t k, Config& config) {
    print_to_both(config, "\nRunning for W = " + std::to_string(w) + " and K = " + std::to_string(k) + "\n");

    auto time_before_parallel_run = std::chrono::steady_clock::now();
    parallel_run(config, w, k, true);
    auto elapsed_seconds_parallel_run = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - time_before_parallel_run).count();

    std::vector<uint64_t> best_order = load_order(w, k, config.min_alpha, config.max_alpha, false);
    std::vector<uint64_t> best_order_reversed = get_reversed_order(best_order);
    double binary_density = density_expected_binary(w, k, best_order);
    double binary_density_reversed = density_expected_binary(w, k, best_order_reversed);

    if (binary_density_reversed < binary_density) {
        best_order = best_order_reversed;
        save_order(w, k, config.min_alpha, config.max_alpha, best_order_reversed, false);
    }

    // this indicates same swapping time as the parallel run
    if (config.max_mins_per_step == std::numeric_limits<uint32_t>::max()) {
        single_run_swapper(w, k, config.min_alpha, config.max_alpha, elapsed_seconds_parallel_run);
    }
    else {
        single_run_swapper(w, k, config.min_alpha, config.max_alpha, config.max_swapper_time_minutes*60);
    }

    std::vector<uint64_t> best_order_after_swap = load_order(w, k, config.min_alpha, config.max_alpha, true);
    std::vector<uint64_t> best_order_after_swap_reversed = get_reversed_order(best_order_after_swap);



    double binary_density_after_swap = density_expected_binary(w, k, best_order_after_swap);
    double binary_density_after_swap_reversed = density_expected_binary(w, k, best_order_after_swap_reversed);

    double dna_density_upper_bound = (binary_density_after_swap + binary_density_after_swap_reversed) / 2;


    //double dna_density = density_expected_dna(w, k, best_order);
    //double dna_density_after_swap = density_expected_dna(w, k, best_order_after_swap);

    //print_to_both(config, "Binary density before swap: " + std::to_string(binary_density) + "\n");
    print_to_both(config, "Binary density: " + std::to_string(binary_density_after_swap) + "\n");
    //print_to_both(config, "DNA density before swap: " + std::to_string(dna_density) + "\n");
    //print_to_both(config, "DNA density after swap: " + std::to_string(dna_density_after_swap) + "\n");

    print_to_both(config, "DNA density upper bound: " + std::to_string(dna_density_upper_bound) + "\n");



    
}


void compute_and_store_densities(uint32_t w, uint32_t k, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& binary_density_best_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_after_swap_map, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map) {
    print_to_both(config, "\nRunning for W = " + std::to_string(w) + " and K = " + std::to_string(k) + "\n");

    auto time_before_parallel_run = std::chrono::steady_clock::now();
    parallel_run(config, w, k, true);
    auto elapsed_seconds_parallel_run = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - time_before_parallel_run).count();

    std::vector<uint64_t> best_order = load_order(w, k, config.min_alpha, config.max_alpha, false);
    std::vector<uint64_t> best_order_reversed = get_reversed_order(best_order);
    double binary_density = density_expected_binary(w, k, best_order);
    double binary_density_reversed = density_expected_binary(w, k, best_order_reversed);

    if (binary_density_reversed < binary_density) {
        print_to_both(config, "Reverse order is better (" +
            std::to_string(binary_density_reversed) + " < " + std::to_string(binary_density) + ")\n");
        best_order = best_order_reversed;
        save_order(w, k, config.min_alpha, config.max_alpha, best_order_reversed, false);
    }

    single_run_swapper(w, k, config.min_alpha, config.max_alpha, elapsed_seconds_parallel_run);
    std::vector<uint64_t> best_order_after_swap = load_order(w, k, config.min_alpha, config.max_alpha, true);
    double binary_density_after_swap = density_expected_binary(w, k, best_order_after_swap);
    double dna_density = density_expected_dna(w, k, best_order);
    double dna_density_after_swap = density_expected_dna(w, k, best_order_after_swap);

    binary_density_map[w][k] = binary_density;
    binary_density_after_swap_map[w][k] = binary_density_after_swap;
    dna_density_map[w][k] = dna_density;
    dna_density_after_swap_map[w][k] = dna_density_after_swap;

    print_to_both(config, "Binary density before swap: " + std::to_string(binary_density) + "\n");
    print_to_both(config, "Binary density after swap: " + std::to_string(binary_density_after_swap) + "\n");
    print_to_both(config, "DNA density before swap: " + std::to_string(dna_density) + "\n");
    print_to_both(config, "DNA density after swap: " + std::to_string(dna_density_after_swap) + "\n");


    // if no previous k -> skip but still save it as "best"
    if (binary_density_after_swap_map[w].count(k - 1) == 0) {
        binary_density_best_map[w][k] = binary_density_after_swap;
        dna_density_best_map[w][k] = dna_density_after_swap;

        // get final time
        auto final_time = std::chrono::steady_clock::now();
        print_to_both(config, "Total time taken: " + std::to_string(std::chrono::duration<double>(final_time - time_before_parallel_run).count()) + " seconds\n");
        return;
    }

    // otherwise try to extend from previous k if it is better
    double previous_k_binary_density = binary_density_best_map[w][k - 1];
    double current_k_binary_density = binary_density_after_swap_map[w][k];

    if (previous_k_binary_density < current_k_binary_density) {
		print_to_both(config, "Previous k is better (" +
            std::to_string(previous_k_binary_density) + " < " + std::to_string(current_k_binary_density) + ")\n");
		improve_from_previous_k(config, k, w, current_k_binary_density, elapsed_seconds_parallel_run);
	}

    // best between the order and the previous k's extension
    best_order = load_order(w, k, config.min_alpha, config.max_alpha, true);
    double best_binary_density = density_expected_binary(w, k, best_order);
    double best_dna_density = density_expected_dna(w, k, best_order);

    binary_density_best_map[w][k] = best_binary_density;
    dna_density_best_map[w][k] = best_dna_density;


    // get final time
    auto final_time = std::chrono::steady_clock::now();
    print_to_both(config, "Total time taken: " + std::to_string(std::chrono::duration<double>(final_time - time_before_parallel_run).count()) + " seconds\n");
}

void big_w_extension(Config& config, uint32_t initial_w, uint32_t k, uint32_t max_w, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map) {
    std::map<uint32_t, double> binary_density_map;
    std::map<uint32_t, double> binary_density_map_reversed;
    print_to_both(config, "\nRunning binary W extension for K = " + std::to_string(k) + ", starting W = " + std::to_string(initial_w) + "\n");

    // Load the best order for the initial w and k
    std::vector<uint64_t> best_order = load_order(initial_w, k, config.min_alpha, config.max_alpha, true);
    std::vector<uint64_t> reversed_order = get_reversed_order(best_order);
    std::unordered_map<uint32_t, double> w_to_density = density_expected_binary_wide(max_w, k, best_order, true);
    std::unordered_map<uint32_t, double> w_to_density_reversed = density_expected_binary_wide(max_w, k, reversed_order, true);

    for (auto& [w, density] : w_to_density) {
        print_to_both(config, "Density for W = " + std::to_string(w) + " is " + std::to_string(density) + "\n");
        print_to_both(config, "Density for W = " + std::to_string(w) + " reversed is " + std::to_string(w_to_density_reversed[w]) + "\n");
        binary_density_map[w] = density;
        binary_density_map_reversed[w] = w_to_density_reversed[w];
        // For w values that are new - use the lower bound
        if (w > initial_w) {
            dna_density_best_map[w][k] = (density + w_to_density_reversed[w]) / 2;
		}
    }
}

void big_k_extension(Config& config, uint32_t w, uint32_t initial_k, std::map<uint32_t, std::map<uint32_t, double>>& dna_density_best_map) {
    std::map<uint32_t, double> binary_density_map;
    std::map<uint32_t, double> binary_density_map_reversed;
    print_to_both(config, "\nRunning binary K extension for W = " + std::to_string(w) + ", starting K = " + std::to_string(initial_k) + "\n");

    // load the best order for w and the initial k
    std::vector<uint64_t> best_order = load_order(w, initial_k, config.min_alpha, config.max_alpha, true);
    std::vector<uint64_t> reversed_order = get_reversed_order(best_order);

    auto starting_time = std::chrono::high_resolution_clock::now();

    for (int extended_k = initial_k + 1; extended_k < 31; extended_k++) {
        auto time_passed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - starting_time).count();
        if (time_passed >= (config.max_mins_per_step * 60)) {
            print_to_both(config, "Reached " + std::to_string(config.max_mins_per_step) + " minutes for W = " + std::to_string(w) + ", stopping further K iterations.\n");
            break;
        }
        print_to_both(config, "\nRunning for W = " + std::to_string(w) + " and K = " + std::to_string(initial_k) + " and extended K = " + std::to_string(extended_k) + "\n");
        double extended_density = density_expected_binary_extend_k(w, initial_k, best_order, extended_k);
        double extended_density_reversed = density_expected_binary_extend_k(w, initial_k, reversed_order, extended_k);
        binary_density_map[extended_k] = extended_density;
        binary_density_map_reversed[extended_k] = extended_density_reversed;
        // For k values that are new - use the lower bound (all are new)
        dna_density_best_map[w][extended_k] = (extended_density + extended_density_reversed) / 2;


        print_to_both(config, "Density for extended K = " + std::to_string(extended_k) + " is " + std::to_string(extended_density) + "\n");
        print_to_both(config, "Density for extended K = " + std::to_string(extended_k) + " reversed is " + std::to_string(extended_density_reversed) + "\n");
    }
}

void improve_from_previous_k(Config& config,  uint32_t k, uint32_t w, double current_best_density, double max_swapping_time) {
    std::vector<uint64_t> order = load_order(w, k-1, config.min_alpha, config.max_alpha, true);
    print_to_both(config, "Extending k and swapping for  (w,k) = (" + std::to_string(w) + ", " + std::to_string(k) + ")\n");

    std::vector<uint64_t> explicit_k_extension_order = get_explicitly_extended_order(order);



    save_order(w, k, 2 + config.min_alpha, config.max_alpha, explicit_k_extension_order, false);
    single_run_swapper(w, k, 2 + config.min_alpha, config.max_alpha, max_swapping_time);
    std::vector<uint64_t> explicit_k_extension_order_after_swap = load_order(w, k, 2 + config.min_alpha, config.max_alpha, true);
    double binary_density = density_expected_binary(w, k, explicit_k_extension_order_after_swap);


    // check if we improved
    if (binary_density > current_best_density) {
        print_to_both(config, "Extending and swapping didn't improve binary density" + std::to_string(binary_density) + " > " + std::to_string(current_best_density) + "\n");
        return;
	}

    // if improved, save the order instead of the previous one (As if it was a normal order)
    print_to_both(config, "Improved binary density from " + std::to_string(current_best_density) + " to " + std::to_string(binary_density) + "\n");
    save_order(w, k, config.min_alpha, config.max_alpha, explicit_k_extension_order_after_swap, true);
}

void save_particular_density_csvs(std::map<uint32_t, std::map<uint32_t, double>>& normal_order, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order, std::map<uint32_t, std::map<uint32_t, double>>& normal_order_on_expected, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order_on_expected)
{
    save_2d_to_csv(normal_order, "normal_order_dna_density", "W", "K", config);
    save_2d_to_csv(specifically_trained_order, "specifically_trained_order_dna_density", "W", "K", config);
    save_2d_to_csv(normal_order_on_expected, "normal_order_dna_density_on_expected", "W", "K", config);
    save_2d_to_csv(specifically_trained_order_on_expected, "specifically_trained_order_on_expected", "W", "K", config);
}



void short_calculate_particular_density(uint32_t w, uint32_t k, Config& config) {
    auto start_time = std::chrono::steady_clock::now();


    print_to_both(config, "\nRunning for W = " + std::to_string(w) + " and K = " + std::to_string(k) + "\n");

    parallel_run_specific(config, w, k, true);

    // Load the best orders
    std::vector<uint64_t> best_specifically_trained_order = load_order_specific(w, k, config.min_alpha, config.max_alpha, false, config.name);
    std::vector<uint64_t> best_specifically_trained_order_reversed = get_reversed_order(best_specifically_trained_order);

    // Calculate densities
    double particular_dna_density_specific = density_particular_dna(w, k, config, best_specifically_trained_order);
    double particular_dna_density_specific_reversed = density_particular_dna(w, k, config, best_specifically_trained_order_reversed);
    if (particular_dna_density_specific_reversed < particular_dna_density_specific) {
        print_to_both(config, "Particular reverse order is better (" + std::to_string(particular_dna_density_specific_reversed) + " < " + std::to_string(particular_dna_density_specific) + ")\n");
        best_specifically_trained_order = best_specifically_trained_order_reversed;
        particular_dna_density_specific = particular_dna_density_specific_reversed;
        save_order_specific(w, k, config.min_alpha, config.max_alpha, best_specifically_trained_order, false, config.name);
    }


    print_to_both(config, "Density for best specifically trained order: " + std::to_string(particular_dna_density_specific) + "\n");



    auto current_time = std::chrono::steady_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
    print_to_both(config, "Total time taken: " + std::to_string(elapsed_seconds) + " seconds\n");
}




void calculate_particular_density_w_fixed_k_varying(uint32_t w, const int max_time_in_seconds, Config& config, std::map<uint32_t, std::map<uint32_t, double>>& normal_order, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order, std::map<uint32_t, std::map<uint32_t, double>>& normal_order_on_expected, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order_on_expected, std::set<std::pair<uint32_t, uint32_t>>& computed_pairs)
{
    // First loop over W, varying K
    auto start_time = std::chrono::steady_clock::now();
    for (uint32_t k = 3; k < 17; k++) {
        auto current_time = std::chrono::steady_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

        if (elapsed_seconds >= max_time_in_seconds) {
            print_to_both(config, "Reached " + std::to_string(config.max_mins_per_step) + " minutes for W = " + std::to_string(w) + ", stopping further K iterations.\n");
            break;
        }

        print_to_both(config, "\nRunning for W = " + std::to_string(w) + " and K = " + std::to_string(k) + "\n");

        // Check if the order already exists; if not, create it
        if (!does_order_exist(w, k, config.min_alpha, config.max_alpha, false)) {
            print_to_both(config, "Order does not exist for W = " + std::to_string(w) + " and K = " + std::to_string(k) + "\n");
            print_to_both(config, "Creating order\n");
            parallel_run(config, w, k, false);
        }

        // Run specific training
        parallel_run_specific(config, w, k, true);

        // Load the best orders
        std::vector<uint64_t> best_order = load_order(w, k, config.min_alpha, config.max_alpha, false);
        std::vector<uint64_t> best_specifically_trained_order = load_order_specific(w, k, config.min_alpha, config.max_alpha, false, config.name);
        std::vector<uint64_t> best_specifically_trained_order_reversed = get_reversed_order(best_specifically_trained_order);

        // Calculate densities
        double particular_dna_density = density_particular_dna(w, k, config, best_order);
        double particular_dna_density_specific = density_particular_dna(w, k, config, best_specifically_trained_order);
        double particular_dna_density_specific_reversed = density_particular_dna(w, k, config, best_specifically_trained_order_reversed);
        if (particular_dna_density_specific_reversed < particular_dna_density_specific) {
            print_to_both(config, "Particular reverse order is better (" + std::to_string(particular_dna_density_specific_reversed) + " < " + std::to_string(particular_dna_density_specific) + ")\n");
            best_specifically_trained_order = best_specifically_trained_order_reversed;
            particular_dna_density_specific = particular_dna_density_specific_reversed;
            save_order_specific(w, k, config.min_alpha, config.max_alpha, best_specifically_trained_order, false, config.name);
        }
        normal_order[w][k] = particular_dna_density;
        specifically_trained_order[w][k] = particular_dna_density_specific;

        print_to_both(config, "Density for best order: " + std::to_string(particular_dna_density) + "\n");
        print_to_both(config, "Density for best specifically trained order: " + std::to_string(particular_dna_density_specific) + "\n");

        double dna_density_expected = density_expected_dna(w, k, best_order);
        double dna_density_expected_specific = density_expected_dna(w, k, best_specifically_trained_order);
        normal_order_on_expected[w][k] = dna_density_expected;
        specifically_trained_order_on_expected[w][k] = dna_density_expected_specific;

        print_to_both(config, "Density for best order on expected: " + std::to_string(dna_density_expected) + "\n");
        print_to_both(config, "Density for best specifically trained order on expected: " + std::to_string(dna_density_expected_specific) + "\n");

        // Mark this cell as computed
        computed_pairs.insert(std::make_pair(w, k));
    }
}

void calculate_particular_density_k_fixed_w_varying(uint32_t k, const int max_time_in_seconds, Config& config, std::set<std::pair<uint32_t, uint32_t>>& computed_pairs, std::map<uint32_t, std::map<uint32_t, double>>& normal_order, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order, std::map<uint32_t, std::map<uint32_t, double>>& normal_order_on_expected, std::map<uint32_t, std::map<uint32_t, double>>& specifically_trained_order_on_expected)
{
    auto start_time = std::chrono::steady_clock::now();
    for (uint32_t w = 3; w < 22; w++) {
        auto current_time = std::chrono::steady_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

        if (elapsed_seconds >= max_time_in_seconds) {
            print_to_both(config, "Reached " + std::to_string(config.max_mins_per_step) + " minutes for K = " + std::to_string(k) + ", stopping further W iterations.\n");
            break;
        }

        // Skip already computed cells
        if (computed_pairs.count(std::make_pair(w, k)) > 0) {
            continue;
        }

        print_to_both(config, "\nRunning for W = " + std::to_string(w) + " and K = " + std::to_string(k) + "\n");

        // Check if the order already exists; if not, create it
        if (!does_order_exist(w, k, config.min_alpha, config.max_alpha, false)) {
            print_to_both(config, "Order does not exist for W = " + std::to_string(w) + " and K = " + std::to_string(k) + "\n");
            // temp: STOP if no order
            print_to_both(config, "STOPPING \n");
            break;

            print_to_both(config, "Creating order\n");
            parallel_run(config, w, k, true);
        }

        // Run specific training
        parallel_run_specific(config, w, k, true);

        // Load the best orders
        std::vector<uint64_t> best_order = load_order(w, k, config.min_alpha, config.max_alpha, false);
        std::vector<uint64_t> best_specifically_trained_order = load_order_specific(w, k, config.min_alpha, config.max_alpha, false, config.name);
        std::vector<uint64_t> best_specifically_trained_order_reversed = get_reversed_order(best_specifically_trained_order);

        // Calculate densities
        double particular_dna_density = density_particular_dna(w, k, config, best_order);
        double particular_dna_density_specific = density_particular_dna(w, k, config, best_specifically_trained_order);
        double particular_dna_density_specific_reversed = density_particular_dna(w, k, config, best_specifically_trained_order_reversed);
        if (particular_dna_density_specific_reversed < particular_dna_density_specific) {
            print_to_both(config, "Particular reverse order is better (" + std::to_string(particular_dna_density_specific_reversed) + " < " + std::to_string(particular_dna_density_specific) + ")\n");
            best_specifically_trained_order = best_specifically_trained_order_reversed;
            particular_dna_density_specific = particular_dna_density_specific_reversed;
            save_order_specific(w, k, config.min_alpha, config.max_alpha, best_specifically_trained_order, false, config.name);

        }

        normal_order[w][k] = particular_dna_density;
        specifically_trained_order[w][k] = particular_dna_density_specific;

        print_to_both(config, "Density for best order: " + std::to_string(particular_dna_density) + "\n");
        print_to_both(config, "Density for best specifically trained order: " + std::to_string(particular_dna_density_specific) + "\n");

        double dna_density_expected = density_expected_dna(w, k, best_order);
        double dna_density_expected_specific = density_expected_dna(w, k, best_specifically_trained_order);
        normal_order_on_expected[w][k] = dna_density_expected;
        specifically_trained_order_on_expected[w][k] = dna_density_expected_specific;

        print_to_both(config, "Density for best order on expected: " + std::to_string(dna_density_expected) + "\n");
        print_to_both(config, "Density for best specifically trained order on expected: " + std::to_string(dna_density_expected_specific) + "\n");

        // Mark this cell as computed
        computed_pairs.insert(std::make_pair(w, k));
    }
}



void swaps_only(Config& config) {
    single_run_swapper_v2(config);
 
}