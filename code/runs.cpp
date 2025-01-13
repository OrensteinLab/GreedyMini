#include <iostream>
#include <chrono>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <limits>
#include <fstream>
#include <filesystem>  // For creating directories
#include "functions.h" 
#include "runs.h"
#include "tools.h"
#include "swapper.h"
#include "results_manager.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "density_and_gc.h"
#include "functions.h"
#include "config.h"
#include "preprocess_particular.h"


// Mutex for protecting shared resources
std::mutex mtx;

// Function to run ecogreed in a single thread
void run_ecogreed(uint32_t W, uint32_t K, Result& best_result, int iterations, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans, Config &config) {
    for (int i = 0; i < iterations; ++i) {
        // TODO: PICK WHICH VERSION IS BEST
        double alpha = sample_alpha(config);
        //double actual_error = err;
        auto [new_ans, gc_count] = ecogreed(W, K, alpha,gc ,nongc ,ans );

        std::lock_guard<std::mutex> guard(mtx);
        if (gc_count < best_result.gc_count) {
            best_result.ans = new_ans;
            best_result.gc_count = gc_count;
            best_result.actual_alpha = alpha;
        }
    }
}
//
//void run_ecogreed_top_3(uint32_t W, uint32_t K, double err, KeepTrackTop3& orders_tracker, int iterations, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans) {
//    for (int i = 0; i < iterations; ++i) {
//        auto [new_ans, gc_count] = ecogreed(W, K, err, gc, nongc, ans);
//
//        std::lock_guard<std::mutex> guard(mtx);
//        orders_tracker.all_gc_counts.push_back(gc_count);
//
//        // Combine top_3_orders and top_3_gc_counts into a vector of pairs
//        std::vector<std::pair<uint64_t, std::vector<uint64_t>>> top_3_pairs;
//        for (int j = 0; j < 3; ++j) {
//            top_3_pairs.emplace_back(orders_tracker.top_3_gc_counts[j], orders_tracker.top_3_orders[j]);
//        }
//
//        // Sort the pairs by the gc_count (first element in the pair)
//        std::sort(top_3_pairs.begin(), top_3_pairs.end(), [](const auto& a, const auto& b) {
//            return a.first < b.first; // Sort based on gc_count
//            });
//
//        // Update the top 3 gc counts and orders after sorting
//        for (int j = 0; j < 3; ++j) {
//            orders_tracker.top_3_gc_counts[j] = top_3_pairs[j].first;
//            orders_tracker.top_3_orders[j] = top_3_pairs[j].second;
//        }
//
//        // Now check if the new gc_count is smaller than the largest (last) of the top 3
//        if (gc_count < orders_tracker.top_3_gc_counts[2]) {
//            orders_tracker.top_3_gc_counts[2] = gc_count;
//            orders_tracker.top_3_orders[2] = new_ans;
//
//            // Resort the top 3 after adding the new gc count and order
//            top_3_pairs.clear();
//            for (int j = 0; j < 3; ++j) {
//                top_3_pairs.emplace_back(orders_tracker.top_3_gc_counts[j], orders_tracker.top_3_orders[j]);
//            }
//
//            std::sort(top_3_pairs.begin(), top_3_pairs.end(), [](const auto& a, const auto& b) {
//                return a.first < b.first;
//                });
//
//            for (int j = 0; j < 3; ++j) {
//                orders_tracker.top_3_gc_counts[j] = top_3_pairs[j].first;
//                orders_tracker.top_3_orders[j] = top_3_pairs[j].second;
//            }
//        }
//    }
//}
//
//



void run_ecogreed_specific(uint32_t W, uint32_t K, Result& best_result, int iterations, std::vector<uint64_t>& allSequences,
    std::vector<uint64_t> gc, std::vector<uint64_t> nongc, std::vector<uint64_t> ans,
    std::vector<std::vector<uint64_t>>& kmer_to_gc_string_id, std::vector<std::vector<uint64_t>>& kmer_to_non_gc_string_id, std::vector<std::vector<uint64_t>>& string_id_to_non_gc_kmers, std::vector<uint64_t>& string_id_to_gc_prefix,
    std::vector<uint64_t>& string_id_to_gc_suffix, Config &config) {


    for (int i = 0; i < iterations; ++i) {
        double alpha = sample_alpha(config);
        auto [final_ans, gc_count] = ecogreed_specific(W, K, alpha, config.path, config.name, allSequences, gc, nongc, ans, kmer_to_gc_string_id, kmer_to_non_gc_string_id, string_id_to_non_gc_kmers,string_id_to_gc_prefix, string_id_to_gc_suffix);
        std::lock_guard<std::mutex> guard(mtx);
        if (gc_count <  best_result.gc_count) {
            best_result.ans = final_ans;
            best_result.gc_count = gc_count;
            best_result.actual_alpha = alpha;
        }
    }
}

void parallel_run(Config &config, uint32_t W, uint32_t K, const bool save_results) {
    ensure_init_lists(W, K);
 

    auto [gc, nongc, ans] = load_init_lists(W, K);

    unsigned int num_cores = config.n_cores;
    std::vector<std::thread> threads;
    Result best_result = { {}, std::numeric_limits<uint64_t>::max(), 0 };

    //TODO: REMOVE LATER
    //num_cores /= 2;



    int total_iterations = config.greedy_mini_runs + (num_cores - (config.greedy_mini_runs % num_cores))%num_cores;

    // Calculate the number of iterations per thread
    int iterations_per_thread = total_iterations / num_cores;


    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Launch threads to run ecogreed in parallel
    for (unsigned int i = 0; i < num_cores; ++i) {
        threads.emplace_back(run_ecogreed, W, K, std::ref(best_result), iterations_per_thread,std::ref(gc), std::ref(nongc), std::ref(ans), std::ref(config));
    }

    // Wait for all threads to finish
    for (auto& t : threads) {
        t.join();
    }


    // Print the actual error of the best result
    std::cout << "Alpha of the best result: " << best_result.actual_alpha << std::endl;


    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_duration = end - start;

    // Save results if the flag is set
    if (save_results) {
        save_order(W, K, config.min_alpha, config.max_alpha, best_result.ans, false);
    }

  
    std::cout << "GreedyMini time: " << total_duration.count() << " seconds, using " << num_cores << " threads to run " << total_iterations << " iterations" << std::endl;

    return;
}

//void parallel_run_top_3(uint32_t W, uint32_t K, double err, const int total_runs) {
//    ensure_init_lists(W, K);
//
//    auto [gc, nongc, ans] = load_init_lists(W, K);
//
//    unsigned int num_cores = std::thread::hardware_concurrency();
//    std::vector<std::thread> threads;
//    KeepTrackTop3 orders_tracker = {
//        {std::vector<uint64_t>{}, std::vector<uint64_t>{}, std::vector<uint64_t>{}}, // top_3_orders (three empty vectors)
//        {std::numeric_limits<uint64_t>::max(),std::numeric_limits<uint64_t>::max(),std::numeric_limits<uint64_t>::max()}, // top_3_gc_counts (initialized to max values)
//        {} // all_gc_counts (empty vector)
//    };
//
//    //TODO: REMOVE LATER
//    num_cores /= 2;
//
//
//
//    int total_iterations = total_runs + (num_cores - (total_runs % num_cores)) % num_cores;
//
//    std::cout << "Total iterations: " << total_iterations << std::endl;
//
//    // Calculate the number of iterations per thread
//    int iterations_per_thread = total_iterations / num_cores;
//
//
//
//    // Launch threads to run ecogreed in parallel
//    for (unsigned int i = 0; i < num_cores; ++i) {
//        threads.emplace_back(run_ecogreed_top_3, W, K, err, std::ref(orders_tracker), iterations_per_thread, std::ref(gc), std::ref(nongc), std::ref(ans));
//    }
//
//    // Wait for all threads to finish
//    for (auto& t : threads) {
//        t.join();
//    }
//
//
//    // print top 3 orders and gc counts
//    for (int i = 0; i < 3; ++i) {
//        std::cout << "Top " << i + 1 << " GC count: " << orders_tracker.top_3_gc_counts[i] << std::endl;
//        std::cout << "Top " << i + 1 << " order: ";
//        std::cout << "[";
//        for (int j = 0; j < orders_tracker.top_3_orders[i].size(); ++j) {
//            std::cout << orders_tracker.top_3_orders[i][j];
//            if (j != orders_tracker.top_3_orders[i].size() - 1) {
//                std::cout << ", ";
//            }
//        }
//        std::cout << "]" << std::endl;
//        std::cout << std::endl;
//    }
//
//    //  print all gc counts 
//    std::cout << "All GC counts: ";
//    std::cout << "[";
//    for (int i = 0; i < orders_tracker.all_gc_counts.size(); ++i) {
//        std::cout << orders_tracker.all_gc_counts[i];
//        if (i != orders_tracker.all_gc_counts.size() - 1) {
//            std::cout << ", ";
//        }
//    }
//    std::cout << "]" << std::endl;
//    return;
//}

//void single_run(uint32_t W, uint32_t K, double err, const bool save_results)
//{
//    ensure_init_lists(W, K);
//    auto [gc, nongc, ans] = load_init_lists(W, K);
//
//
//
//    // Start timing
//    auto start = std::chrono::high_resolution_clock::now();
//
//    // Call the ecogreed function
//    auto [new_ans, gc_count] = ecogreed(W, K, err,gc ,nongc ,ans );
//
//    // Stop timing
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> duration = end - start;
//
//     // Save results if the flag is set
//    if (save_results) {
//
//        // Save log file as CSV
//        std::ofstream log_file("logs/single_" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(err) + ".csv");
//
//        // Write the header (titles)
//        log_file << "Total running time w/o preprocessing (seconds),Best GC count\n";
//
//        // Write the corresponding values
//        log_file << duration.count() << "," << gc_count << "\n";
//
//        // Close the log file
//        log_file.close();
//
//        // Save result to a file
//        save_order(W, K, err, new_ans, false);
//        //save_vector_to_file(ans, "output/ordering/" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(err) + ".bin");
//    }
//
//
//    // Print the gamechanger count
//    std::cout << "Gamechanger count: " << gc_count << std::endl;
//
//    // Print the time taken
//    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
//    return;
//
//}

//
//void single_run_specific(uint32_t W, uint32_t K, double err, const bool save_results, std::string path, std::string name) {
//    ensure_init_lists_specific(W, K, path, name);
//
//    std::vector<uint64_t> allSequences = load_all_sequences_particular(W, K, path);
//    auto [gc, nongc, initial_ans, kmer_to_gc_string_id, kmer_to_non_gc_string_id, string_id_to_non_gc_kmers, string_id_to_gc_prefix, string_id_to_gc_suffix] = load_init_lists_specific(W, K, name);
//
//    // Start timing
//    auto start = std::chrono::high_resolution_clock::now();
//
//    // Call the ecogreed function
//    auto [ans, gc_count] = ecogreed_specific(W, K, err, path, name, allSequences, gc, nongc, initial_ans, kmer_to_gc_string_id, kmer_to_non_gc_string_id, string_id_to_non_gc_kmers, string_id_to_gc_prefix, string_id_to_gc_suffix);
//
//    // Stop timing
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> duration = end - start;
//
//    // Save results if the flag is set
//    if (save_results) {
//        // make subfolders for 'name'
//        std::filesystem::create_directories("output/ordering/" + name);
//        std::filesystem::create_directories("logs/" + name);
//
//        // Save log file as CSV
//        std::ofstream log_file("logs/" + name + "/single_" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(err) + ".csv");
//
//        // Write the header (titles)
//        log_file << "Total running time w/o preprocessing (seconds),Best GC count\n";
//
//        // Write the corresponding values
//        log_file << duration.count() << "," << gc_count << "\n";
//
//        // Close the log file
//        log_file.close();
//
//        // Save result to a file
//        save_order_specific(W, K, err, ans, false, name);
//    }
//
//
//
//
//    // Print the time taken
//    std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
//    return;
//
//}

void parallel_run_specific(Config &config, uint32_t W, uint32_t K, const bool save_results) {

    ensure_sequence_is_processed(config, W, K);
    
    ensure_init_lists_specific(W, K, config);

    unsigned int num_cores = config.n_cores;
    std::vector<std::thread> threads;
    Result best_result = { {}, std::numeric_limits<uint64_t>::max() };

    //TODO: REMOVE LATER
    //num_cores /= 2;


    int total_iterations = config.greedy_mini_runs + (num_cores - (config.greedy_mini_runs % num_cores)) % num_cores;

  

    // Calculate the number of iterations per thread
    int iterations_per_thread = total_iterations / num_cores;


    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<uint64_t> allSequences = load_all_sequences_particular(W, K, config);
    auto [gc, nongc, initial_ans, kmer_to_gc_string_id, kmer_to_non_gc_string_id, string_id_to_non_gc_kmers, string_id_to_gc_prefix, string_id_to_gc_suffix] = load_init_lists_specific(W, K, config);

    // Launch threads to run ecogreed in parallel
    for (unsigned int i = 0; i < num_cores; ++i) {
        threads.emplace_back(run_ecogreed_specific, W, K, std::ref(best_result), iterations_per_thread, std::ref(allSequences), std::ref(gc), std::ref(nongc), std::ref(initial_ans), std::ref(kmer_to_gc_string_id), std::ref(kmer_to_non_gc_string_id), std::ref(string_id_to_non_gc_kmers), std::ref(string_id_to_gc_prefix), std::ref(string_id_to_gc_suffix), std::ref(config));
    }

    // Wait for all threads to finish
    for (auto& t : threads) {
        t.join();
    }

    // Print the actual alpha of the best result
    std::cout << "Alpha of the best result: " << best_result.actual_alpha << std::endl;



    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_duration = end - start;

    // Save results if the flag is set
    if (save_results) {
        //  make subfolders for 'name'
        std::filesystem::create_directories("output/minimizers/" + config.name);

        // Save result to a file
        save_order_specific(W, K, config.min_alpha, config.max_alpha, best_result.ans, false, config.name);
    }


    std::cout << "GreedyMiniParticular time: " << total_duration.count() << " seconds, using " << num_cores << " threads to run " << total_iterations << " iterations" << std::endl;

}

void single_run_swapper(uint32_t W, uint32_t K, double min_alpha, double max_alpha, const double max_time_seconds)
{
    // Load original order from file
    std::vector<uint64_t> order = load_order(W, K, min_alpha, max_alpha, false);

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Number of threads to use
    size_t num_threads = std::thread::hardware_concurrency();
   
    // TODO: chjange later
    num_threads /= 2;

    // Vector to store results from each thread
    std::vector<std::pair<std::vector<uint64_t>, uint64_t>> results(num_threads);

    // Vector of threads
    std::vector<std::thread> threads;


    for (size_t i = 0; i < num_threads; ++i) {
        threads.emplace_back([&, i]() {
            // Each thread needs its own copy of order
            std::vector<uint64_t> order_copy = order;

            bool verbose;
            if (i == 0) {
                // during testing it was true for this only
                verbose = true;
            }
            else {
                verbose = false;
            }

            // Run swapper_f_v5
            auto result = swapper_f(W, K, order_copy, max_time_seconds, verbose);

            // Store result in results[i]
            results[i] = result;

            // Optional: Print thread completion
            //std::cout << "Thread " << i + 1 << " completed with GC count: " << result.second << std::endl;
            });
    }

    // Wait for all threads to finish
    for (auto& t : threads) {
        t.join();
    }

    // Find the result with the lowest swapped_gc_count
    auto best_it = std::min_element(results.begin(), results.end(),
        [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

    auto worst_it = std::max_element(results.begin(), results.end(),
		[](const auto& a, const auto& b) {
			return a.second < b.second;
		});

    std::vector<uint64_t> swapped_ans = best_it->first;
    uint64_t swapped_gc_count = best_it->second;
    uint64_t worst_swapped_gc_count = worst_it->second;

    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Swapper time: " << duration.count() << " seconds. " << std::endl;
    save_order(W, K, min_alpha, max_alpha, swapped_ans, true);

}



void single_run_swapper_v2(Config& config)
{

    uint64_t max_time_seconds = 60;
    if (config.max_swapper_time_minutes != std::numeric_limits<uint32_t>::max()) {
        max_time_seconds = config.max_swapper_time_minutes * 60;
        print_to_both(config, "Running swapper with max time: " + std::to_string(max_time_seconds) + " seconds\n");

    }
    else {
        print_to_both(config, "Running swapper with 1 minute max time limit (default)\n");
    }







    // Load original order from file
    std::vector<uint64_t> order = load_order_path(config.path);

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Number of threads to use
    size_t num_threads = std::thread::hardware_concurrency();

    // TODO: chjange later
    num_threads /= 2;


    // Vector to store results from each thread
    std::vector<std::pair<std::vector<uint64_t>, uint64_t>> results(num_threads);

    // Vector of threads
    std::vector<std::thread> threads;


    for (size_t i = 0; i < num_threads; ++i) {
        threads.emplace_back([&, i]() {
            // Each thread needs its own copy of order
            std::vector<uint64_t> order_copy = order;

            bool verbose;
            if (i == 0) {
                // during testing it was true for this only
                verbose = true;
                
            }
            else {
                verbose = false;
            }

            // Run swapper_f_v5
            auto result = swapper_f_v2(config.w, config.k, order_copy, max_time_seconds, verbose, (i==0));

            // Store result in results[i]
            results[i] = result;

            // Optional: Print thread completion
            //std::cout << "Thread " << i + 1 << " completed with GC count: " << result.second << std::endl;
            });
    }

    // Wait for all threads to finish
    for (auto& t : threads) {
        t.join();
    }

    // Find the result with the lowest swapped_gc_count
    auto best_it = std::min_element(results.begin(), results.end(),
        [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

    auto worst_it = std::max_element(results.begin(), results.end(),
        [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

    std::vector<uint64_t> swapped_ans = best_it->first;
    uint64_t swapped_gc_count = best_it->second;
    uint64_t worst_swapped_gc_count = worst_it->second;

    double density = calc_density(swapped_gc_count, config.k, config.w);
    std::cout << "gc_count: " << swapped_gc_count << std::endl;
    print_to_both(config, "Density after swap:" + std::to_string(density) + "\n");

    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Swapper time: " << duration.count() << " seconds. " << std::endl;


    // output path is same as input path but with _swapped added after .bin
    std::string output_path = config.path;
    output_path.insert(output_path.find(".bin"), "_swapped");



    save_order_path(output_path, swapped_ans);

}





//void multiple_runs_swapper(uint32_t W, uint32_t K, double err, const double max_time_seconds, const uint32_t n_repetitions)
//{
//    // Load original order from file
//    std::vector<uint64_t> order = load_order(W, K, err, false);
//
//    // Start timing
//    auto start = std::chrono::high_resolution_clock::now();
//
//    // Number of threads to use
//    size_t num_threads = std::thread::hardware_concurrency();
//
//    // TODO: change later
//    num_threads /= 2;
//
//    // Vector to store results from each thread: {best_result, worst_result}
//    std::vector<std::pair<std::pair<std::vector<uint64_t>, uint64_t>, std::pair<std::vector<uint64_t>, uint64_t>>> results(num_threads);
//
//    // Vector of threads
//    std::vector<std::thread> threads;
//
//    for (size_t i = 0; i < num_threads; ++i) {
//        threads.emplace_back([&, i]() {
//            // Each thread needs its own copy of order
//            std::vector<uint64_t> order_copy = order;
//
//            // Initialize best and worst GC counts
//            uint64_t best_gc_count = UINT64_MAX;  // Best starts high, we want to minimize
//            uint64_t worst_gc_count = 0;          // Worst starts low, we want to maximize
//
//            std::vector<uint64_t> best_order;
//            std::vector<uint64_t> worst_order;
//
//            for (uint32_t j = 0; j < n_repetitions; ++j) {
//                // Run swapper_f_v5
//                auto result = swapper_f_v5(W, K, order_copy, max_time_seconds);
//
//                // Check if the current result is better (lower GC count)
//                if (result.second < best_gc_count) {
//                    best_gc_count = result.second;
//                    best_order = result.first;
//                }
//
//                // Check if the current result is worse (higher GC count)
//                if (result.second > worst_gc_count) {
//                    worst_gc_count = result.second;
//                    worst_order = result.first;
//                }
//            }
//
//            // Store both the best and worst result for this thread
//            results[i] = { {best_order, best_gc_count}, {worst_order, worst_gc_count} };
//            });
//    }
//
//    // Wait for all threads to finish
//    for (auto& t : threads) {
//        t.join();
//    }
//
//    // Find the overall best result across threads (lowest GC count)
//    auto best_it = std::min_element(results.begin(), results.end(),
//        [](const auto& a, const auto& b) {
//            return a.first.second < b.first.second;  // Compare the best GC counts
//        });
//
//    // Find the overall worst result across threads (highest GC count)
//    auto worst_it = std::max_element(results.begin(), results.end(),
//        [](const auto& a, const auto& b) {
//            return a.second.second < b.second.second;  // Compare the worst GC counts
//        });
//
//    std::vector<uint64_t> swapped_ans = best_it->first.first;
//    uint64_t swapped_gc_count = best_it->first.second;
//    uint64_t worst_swapped_gc_count = worst_it->second.second;
//
//    // Stop timing
//    auto end = std::chrono::high_resolution_clock::now();
//
//    // Calculate the time taken
//    std::chrono::duration<double> duration = end - start;
//
//    // Print the time taken and the best/worst GC counts
//    std::cout << "Swapper time: " << duration.count() << " seconds for max time of " << max_time_seconds << " seconds" << std::endl;
//    std::cout << "Best Swapper GC count: " << swapped_gc_count << std::endl;
//    std::cout << "Worst Swapper GC count: " << worst_swapped_gc_count << std::endl;
//
//    // Save the swapped order (best result)
//    save_order(W, K, err, swapped_ans, true);
//
//    // Save the log file
//    std::ofstream log_file("logs/swapper_" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(err) + ".csv");
//
//    // Write the header (titles)
//    log_file << "Total running time w/o preprocessing (seconds),Best GC count,Worst GC count\n";
//
//    // Write the corresponding values
//    log_file << duration.count() << "," << swapped_gc_count << "," << worst_swapped_gc_count << "\n";
//}