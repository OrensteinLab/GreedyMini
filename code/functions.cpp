#include <unordered_set>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <unordered_map>
#include <tuple>
#include <numeric> 
#include <chrono>
#include "functions.h" 
#include "tools.h"
#include <boost/multiprecision/cpp_int.hpp>  // Include Boost Multiprecision for arbitrary-precision integers
#include <boost/multiprecision/cpp_dec_float.hpp> // Include multiprecision floating point
//#include "gc_counter.h"
#include "density_and_gc.h"
#include <chrono>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <random>  

std::mutex load_mutex;




void update_counters(uint32_t w, uint32_t k, uint64_t nextbest, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc, std::vector<uint64_t>& ans, uint64_t rank) {
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;
    uint64_t max_pos = 1ULL << (w + k);

    auto [block, dep] = nextblock(0, w, k, nextbest, n_kmers, max_pos); 
    while (block < (1ULL << (w + k))) {  
        std::unordered_set<uint64_t> kmers;
        bool skip = false;
        for (uint32_t i = 0; i < dep; ++i) {  
            uint64_t leftkmer = (block >> (w - i)) & n_kmers_mask;
            skip = (ans[leftkmer] < rank);
            if (skip) {  
                std::tie(block, dep) = nextblock(((block >> (w - i)) + 1) << (w - i), w, k, nextbest, n_kmers, max_pos);  // increase leftkmer by 1 and find next block after
                break;
            }
            else {
                kmers.insert(leftkmer);
            }
        }
        if (!skip) {
            uint64_t subt = rightdfs_wrapper(nextbest, kmers, dep, w, k, gc, nongc);
            uint64_t pref = block >> w;
            gc[pref] -= subt;
            if (!kmers.empty()) {
                kmers.erase(pref);
                for (uint64_t km : kmers) {
                    nongc[km] -= subt;
                }
            }
            std::tie(block, dep) = nextblock(block + (1ULL << (w - dep)), w, k, nextbest, n_kmers, max_pos);  // finding next block
        }
    }
}


std::pair<std::vector<uint64_t>, uint64_t> ecogreed(uint32_t w, uint32_t k, double err, std::vector<uint64_t> gc, std::vector<uint64_t> nongc, std::vector<uint64_t> ans) {
    uint64_t n_kmers = 1ULL << k;


    for (uint64_t rank = 0; rank < n_kmers; ++rank) {
        uint64_t nextbest = randcand(k, err, gc, nongc, rank);
        //uint64_t nextbest = randcand_v2(k, gc, nongc);
        //uint64_t nextbest = randcand_v3(k, gc, nongc);
        //std::cout << "WARNING, USING RANDCAND_V3" << std::endl;

        if (nextbest == n_kmers) {
            break;
        }
        ans[nextbest] = rank;
        update_counters(w, k, nextbest, gc, nongc, ans, rank);
    }
    uint64_t gc_count = prob_gc(w, k, ans);
    return { ans, gc_count };

}

void make_init_lists(uint32_t w, uint32_t k) {
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    auto [gc, nongc, ans] = init_lists(w, k);


    // Save the lists to files
    save_vector_to_file(gc, "temp/gc_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    save_vector_to_file(nongc, "temp/nongc_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    save_vector_to_file(ans, "temp/ans_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");

    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_duration = end - start;

    // Print the time taken
    std::cout << "Preprocessing time for W: " << w << " and K: " << k << " is " << total_duration.count() << " seconds" << std::endl;

    // Save the time to a log file
    std::ofstream log_file("logs/preprocessing_time" + std::to_string(w) + "_" + std::to_string(k) + ".txt");
    log_file <<  total_duration.count(); // in seconds
    log_file.close();
}

std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>> load_init_lists(uint32_t w, uint32_t k) {
    std::lock_guard<std::mutex> guard(load_mutex);
    std::vector<uint64_t> gc = load_vector_from_file("temp/gc_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    std::vector<uint64_t> nongc = load_vector_from_file("temp/nongc_" + std::to_string(w) + "_" + std::to_string(k) +  ".bin");
    std::vector<uint64_t> ans = load_vector_from_file("temp/ans_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");

    return { gc, nongc, ans };
}

bool check_init_lists(uint32_t w, uint32_t k) {
    return file_exists("temp/gc_" + std::to_string(w) + "_" + std::to_string(k)  + ".bin") &&
        file_exists("temp/nongc_" + std::to_string(w) + "_" + std::to_string(k)  + ".bin") &&
        file_exists("temp/ans_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
}

void ensure_init_lists(uint32_t w, uint32_t k) {
	if (!check_init_lists(w, k)) {
		make_init_lists(w, k);
	}
}

void initdfs(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t lev, uint32_t w, uint32_t k,
    std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc) {
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;
    bool flag = kmers.find(kmer) == kmers.end();

    if (lev == w) {
        if (flag) {
            gc[kmer]++;
        }
    }
    else {
        if (flag) {
            nongc[kmer] += 1ULL << (w - lev);
            kmers.insert(kmer);
        }
        uint64_t nextkmer = (kmer << 1) & n_kmers_mask;
        initdfs(nextkmer, kmers, lev + 1, w, k, gc, nongc);
        initdfs(nextkmer + 1, kmers, lev + 1, w, k, gc, nongc);
        if (flag) {
            kmers.erase(kmer);
        }
    }
}

std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>> init_lists(uint32_t w, uint32_t k) {
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;

    std::vector<uint64_t> gc(n_kmers, 0);  // (live windows) becoming gamechangers if i is now chosen
    std::vector<uint64_t> nongc(n_kmers, 0); // (live windows) becoming non-gamechangers if i is now chosen
    std::vector<uint64_t> ans(n_kmers, n_kmers); // current partial order

    for (uint64_t pref = 0; pref < n_kmers; ++pref) { // fill the arrays gc and nongc
        gc[pref] += 1ULL << w;
        std::unordered_set<uint64_t> kmers = { pref };
        uint64_t k1 = (2 * pref) & n_kmers_mask;
        initdfs(k1, kmers, 1, w, k, gc, nongc);
        initdfs(k1 + 1, kmers, 1, w, k, gc, nongc);
    }

    return { gc, nongc, ans };
}

uint64_t randcand(uint32_t k, double initial_err, const std::vector<uint64_t>& gc, const std::vector<uint64_t>& nongc, uint64_t rank) {
    uint64_t n_kmers = 1ULL << k;


    // Check for freeride case
    for (uint64_t i = 0; i < n_kmers; ++i) {
        if ((gc[i] == 0) && (nongc[i] > 0)) {
            return i;
        }
    }

    // Find the best nongc/gc ratio
    double best_ratio = 0.0;
    for (uint64_t i = 0; i < n_kmers; ++i) {
        if (gc[i] > 0) {
            double ratio = static_cast<double>(nongc[i]) / gc[i];
            if (ratio > best_ratio) {
                best_ratio = ratio;
            }
        }
    }

    // TODO: decide which version is better
    // calculate error, it starts at initial_err and ends at 1-e^(-5)
    double err = 1 + (initial_err - 1) * std::exp(-5.0 * (static_cast<double>(rank) * k / n_kmers));
    //double err = initial_err;

    // If best_ratio is greater than 0, form a candidate list
    if (best_ratio > 0.0) {
        std::vector<uint64_t> candidates;
        for (uint64_t i = 0; i < n_kmers; ++i) {
            if (gc[i] > 0) {
                double ratio = static_cast<double>(nongc[i]) / gc[i];
                if (ratio >= best_ratio * err) {
                    candidates.push_back(i);
                }
            }
        }


        // Random number generator
        std::random_device rd;  // Non-deterministic seed (can also use a fixed seed for reproducibility)
        std::mt19937 gen(rd()); // Mersenne Twister engine seeded with random device



        if (!candidates.empty()) {
            std::uniform_int_distribution<uint64_t> dist(0, candidates.size() - 1);



            return candidates[dist(gen)];
        }
    }



    // If no valid candidates, return n_kmers to indicate stopping
    return n_kmers;
}

uint64_t randcand_v2(uint32_t k, const std::vector<uint64_t>& gc, const std::vector<uint64_t>& nongc) {
    uint64_t n_kmers = 1ULL << k;

    // Check for freeride case
    for (uint64_t i = 0; i < n_kmers; ++i) {
        if ((gc[i] == 0) && (nongc[i] > 0)) {
            return i;
        }
    }

    // Find the top two nongc/gc ratios
    double best_ratio = 0.0;
    double second_best_ratio = 0.0;
    uint64_t best_candidate = n_kmers;
    uint64_t second_best_candidate = n_kmers;

    for (uint64_t i = 0; i < n_kmers; ++i) {
        if (gc[i] > 0) {
            double ratio = static_cast<double>(nongc[i]) / gc[i];
            if (ratio > best_ratio) {
                second_best_ratio = best_ratio;
                second_best_candidate = best_candidate;
                best_ratio = ratio;
                best_candidate = i;
            }
            else if (ratio > second_best_ratio) {
                second_best_ratio = ratio;
                second_best_candidate = i;
            }
        }
    }

    // If we have valid best and second-best candidates, pick one at random
    if (best_candidate < n_kmers && second_best_candidate < n_kmers) {
        // Random number generator
        std::random_device rd;  // Non-deterministic seed
        std::mt19937 gen(rd()); // Mersenne Twister engine seeded with random device

        // 50-50 chance to pick either best or second best
        std::uniform_int_distribution<uint64_t> dist(0, 1);
        uint64_t selected_candidate = dist(gen) == 0 ? best_candidate : second_best_candidate;

        return selected_candidate;
    }

    // If no valid candidates, return n_kmers to indicate stopping
    return n_kmers;
}


uint64_t randcand_v3(uint32_t k, const std::vector<uint64_t>& gc, const std::vector<uint64_t>& nongc) {
    uint64_t n_kmers = 1ULL << k;


    // Check for freeride case
    for (uint64_t i = 0; i < n_kmers; ++i) {
        if ((gc[i] == 0) && (nongc[i] > 0)) {
            return i;
        }
    }

    // If first time selecting - > select a random kmer
    bool has_zero_gc = false;
    for (uint64_t i = 0; i < n_kmers; ++i) {
        if (gc[i] == 0) {
            has_zero_gc = true;
            break;
        }
    }

    if (!has_zero_gc) {
        // Random number generator
        std::random_device rd;  // Non-deterministic seed
        std::mt19937 gen(rd()); // Mersenne Twister engine seeded with random device

        std::uniform_int_distribution<uint64_t> dist(0, n_kmers - 1);
        return dist(gen); // Return a random k-mer
    }

    // Find the top two nongc/gc ratios
    double best_ratio = 0.0;
    double second_best_ratio = 0.0;
    uint64_t best_candidate = n_kmers;
    uint64_t second_best_candidate = n_kmers;

    for (uint64_t i = 0; i < n_kmers; ++i) {
        if (gc[i] > 0) {
            double ratio = static_cast<double>(nongc[i]) / gc[i];
            if (ratio > best_ratio) {
                second_best_ratio = best_ratio;
                second_best_candidate = best_candidate;
                best_ratio = ratio;
                best_candidate = i;
            }
            else if (ratio > second_best_ratio) {
                second_best_ratio = ratio;
                second_best_candidate = i;
            }
        }
    }

    // If we have valid best and second-best candidates, pick one at random
    if (best_candidate < n_kmers && second_best_candidate < n_kmers) {
        // Random number generator
        std::random_device rd;  // Non-deterministic seed
        std::mt19937 gen(rd()); // Mersenne Twister engine seeded with random device

        // 50-50 chance to pick either best or second best
        std::uniform_int_distribution<uint64_t> dist(0, 1);
        uint64_t selected_candidate = dist(gen) == 0 ? best_candidate : second_best_candidate;

        return selected_candidate;
    }

    // If no valid candidates, return n_kmers to indicate stopping
    return n_kmers;
}





// if processed                                         0000011100KMER0000
// next will process number that comes after            0000011100KMER0000 + 10000
// since we already processed all the subtree below it, we can just skip it
std::pair<uint64_t, uint32_t> nextblock(uint64_t num, uint32_t w, uint32_t k, uint64_t kmer, uint64_t n_kmers, uint64_t pos) {
    uint32_t depth = w + k;
    // pos is initially 2^(w+k)

    for (uint32_t i = 0; i <= w; ++i) {
        uint32_t cur_shift = w + k - i;
        uint64_t cur = (num & ((1ULL << cur_shift) - 1)) >> (w - i);

        uint64_t temp;
        if (cur < kmer) {
            temp = ((num >> (w - i)) + kmer - cur) << (w - i);
        }
        else if (cur > kmer) {
            temp = ((num >> (w - i)) + kmer + n_kmers - cur) << (w - i);
        }
        else {
            temp = num;
        }

        if (temp < pos) {
            pos = temp;
            depth = i;
        }
    }

    return { pos, depth };
}

uint64_t rightdfs_wrapper(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t lev, uint32_t w, uint32_t k, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc) {
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;

    if (lev == 0) {
        return rightdfs_lev_zero(kmer, kmers, w, n_kmers_mask, gc, nongc);
    }
    else {
        uint32_t lev_left = w - lev;
        return rightdfs(kmer, kmers, lev_left, n_kmers_mask, gc, nongc);
    }
}


// This function is used when we treat level 0 only - we need to update the gc[] sicne it doesnt update in the main function, and update nongc for all varying kmers as well
uint64_t rightdfs_lev_zero(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t w, uint64_t n_kmers_mask, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc) {
    if (nongc[kmer] == 0 && gc[kmer] == 0) {
        return 0;  // Stop as no more live windows contain kmer
    }
    else {
        uint64_t nextkmer = (2 * kmer) & n_kmers_mask;
        uint64_t sub = 0;
        kmers.insert(kmer);
        sub = rightdfs(nextkmer, kmers, w - 1, n_kmers_mask, gc, nongc) +
            rightdfs(nextkmer + 1, kmers, w - 1, n_kmers_mask, gc, nongc);
        //gc[kmer] -= sub;
        kmers.erase(kmer);
        return sub;
    }
}
// This function is used for all levels except 0
// it's job is to update gc when we reach a node -> in case it was a suffix gc
// Updates nonegc along the way  while keeping tack of the kmers we have seen
uint64_t rightdfs(uint64_t kmer, std::unordered_set<uint64_t>& kmers, uint32_t lev_left, uint64_t n_kmers_mask, std::vector<uint64_t>& gc, std::vector<uint64_t>& nongc) {
    if (nongc[kmer] == 0 && gc[kmer] == 0) {
        return 0;  // Stop as no more live windows contain kmer
    }
    else {
        bool flag = kmers.find(kmer) == kmers.end();

        if (lev_left == 0) {
            gc[kmer] -= flag;
            return 1;
        }
        else {
            uint64_t nextkmer = (2 * kmer) & n_kmers_mask;

            uint64_t sub = 0;
            // only remove from nongc first time we see a kmer (all the subtree below it)
            if (flag) {
                kmers.insert(kmer);
                sub = rightdfs(nextkmer, kmers, lev_left - 1, n_kmers_mask, gc, nongc) +
                    rightdfs(nextkmer + 1, kmers, lev_left - 1, n_kmers_mask, gc, nongc);
                nongc[kmer] -= sub;
                kmers.erase(kmer);
            }
            else {
                sub = rightdfs(nextkmer, kmers, lev_left - 1, n_kmers_mask, gc, nongc) +
                    rightdfs(nextkmer + 1, kmers, lev_left - 1, n_kmers_mask, gc, nongc);
            }

            return sub;
        }
    }
}






// TODO:
// Just keep reference from each kmer to which strings it is in (doesnt matter what it is in them)
// And keep from each string what non-gc and what gc kmers it has
std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>,
    std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>, std::vector<uint64_t>, std::vector<uint64_t>>
    init_lists_specific(uint32_t w, uint32_t k, Config& config) {

    auto [odd, even] = load_odd_even_pair_from_file(config, w, k);
    uint64_t n_odd = odd.size();

    // for GC/nongc we add 1 more for empty kmer which we ignore later
    std::vector<uint64_t> gc = std::vector<uint64_t>((1ULL << k) + 1, 0);
    std::vector<uint64_t> nongc = std::vector<uint64_t>((1ULL << k) + 1, 0);
    std::vector<uint64_t> ans = std::vector<uint64_t>(1ULL << k, 1ULL << k);

    // Initialize kmer_to_gc_string_id as a vector of size 2^k, each containing a vector of size 10
    std::vector<std::vector<uint64_t>> kmer_to_gc_string_id =
        std::vector<std::vector<uint64_t>>(1ULL << k);

    std::vector<std::vector<uint64_t>> kmer_to_non_gc_string_id =
        std::vector<std::vector<uint64_t>>(1ULL << k);

    // Initialize string_id_to_non_gc_kmers as a vector of size 10, each containing a vector of size w
    std::vector<std::vector<uint64_t>> string_id_to_non_gc_kmers =
        std::vector<std::vector<uint64_t>>(n_odd);

    // 0-2^k-1 to actual kmers, if 2^k then it means empty
    std::vector<uint64_t> string_id_to_gc_prefix = std::vector<uint64_t>(n_odd, 1ULL << k);
    std::vector<uint64_t> string_id_to_gc_suffix = std::vector<uint64_t>(n_odd, 1ULL << k);


    // loop over all odd numbers and update the vectors
    for (uint64_t i = 0; i < n_odd; ++i) {
		uint64_t oddNumber = odd[i];

        uint64_t firstKmer = oddNumber >> w;
        
        
        // update prefix GC
        gc[firstKmer] += 1;
        kmer_to_gc_string_id[firstKmer].push_back(i);
        string_id_to_gc_prefix[i] = firstKmer;


        // create a set of all seen kmers
        std::unordered_set<uint64_t> kmers = { firstKmer };

        // check if last kmer is unique and add all non_gc_kmers to the list
        for (uint32_t j = 1; j < w; ++j) {

            uint64_t kmer_at_j = (oddNumber >> (w - j)) & ((1ULL << k) - 1);

            // in case of a new kmer
            if (kmers.find(kmer_at_j) == kmers.end()) {
                nongc[kmer_at_j] += 1;
                string_id_to_non_gc_kmers[i].push_back(kmer_at_j);
                kmer_to_non_gc_string_id[kmer_at_j].push_back(i);
                kmers.insert(kmer_at_j);
            }
		}

        // If last kmer is unique -> it is a GC, if not we handles it previously
        uint64_t lastKmer = oddNumber & ((1ULL << k) - 1);

        if (kmers.find(lastKmer) == kmers.end()) {
			gc[lastKmer] += 1;
            kmer_to_gc_string_id[lastKmer].push_back(i);
            string_id_to_gc_suffix[i] = lastKmer;
        }	
	}

    return { gc, nongc, ans, kmer_to_gc_string_id,kmer_to_non_gc_string_id, string_id_to_non_gc_kmers, string_id_to_gc_prefix, string_id_to_gc_suffix };
}

std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>,
    std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>, std::vector<uint64_t>, std::vector<uint64_t>>
    load_init_lists_specific(uint32_t w, uint32_t k, Config& config) {

    std::lock_guard<std::mutex> guard(load_mutex);

    std::vector<uint64_t> gc = load_vector_from_file("temp/" + config.name + "/specific_gc_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    std::vector<uint64_t> nongc = load_vector_from_file("temp/" + config.name + "/specific_nongc_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    std::vector<uint64_t> ans = load_vector_from_file("temp/" + config.name + "/specific_ans_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");

    std::vector<uint64_t> string_id_to_gc_prefix = load_vector_from_file("temp/" +
        config.name + "/specific_string_id_to_gc_prefix_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    std::vector<uint64_t> string_id_to_gc_suffix = load_vector_from_file("temp/" +
        config.name + "/specific_string_id_to_gc_suffix_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");

    // Load the vectors of vectors
    std::vector<std::vector<uint64_t>> kmer_to_gc_string_id =
        load_vector_of_vectors_from_file("temp/" + config.name + "/specific_kmer_to_gc_string_id_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");

    std::vector<std::vector<uint64_t>> kmer_to_non_gc_string_id =
        load_vector_of_vectors_from_file("temp/" + config.name + "/specific_kmer_to_non_gc_string_id_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");

    std::vector<std::vector<uint64_t>> string_id_to_non_gc_kmers =
        load_vector_of_vectors_from_file("temp/" + config.name + "/specific_string_id_to_non_gc_kmers_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");




    return { gc, nongc, ans, kmer_to_gc_string_id,kmer_to_non_gc_string_id, string_id_to_non_gc_kmers, string_id_to_gc_prefix, string_id_to_gc_suffix };
}

bool check_init_lists_specific(uint32_t w, uint32_t k, Config& config) {
    std::string base_path = "temp/" + config.name + "/";
    std::string suffix = "_" + std::to_string(w) + "_" + std::to_string(k) + ".bin";

    return file_exists(base_path + "specific_gc" + suffix) &&
        file_exists(base_path + "specific_nongc" + suffix) &&
        file_exists(base_path + "specific_ans" + suffix) &&
        file_exists(base_path + "specific_kmer_to_gc_string_id" + suffix) &&
        file_exists(base_path + "specific_kmer_to_non_gc_string_id" + suffix) &&
        file_exists(base_path + "specific_string_id_to_non_gc_kmers" + suffix) && 
        file_exists(base_path + "specific_string_id_to_gc_prefix" + suffix) &&
        file_exists(base_path + "specific_string_id_to_gc_suffix" + suffix);

       
        

}

void ensure_init_lists_specific(uint32_t w, uint32_t k, Config& config) {
	if (!check_init_lists_specific(w, k, config)) {
		make_init_lists_specific(w, k, config);
	}
}

void make_init_lists_specific(uint32_t w, uint32_t k, Config& config) {
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Initialize the lists
    auto [gc, nongc, ans, kmer_to_gc_string_id, kmer_to_non_gc_string_id, string_id_to_non_gc_kmers, string_id_to_gc_prefix, string_id_to_gc_suffix] = init_lists_specific(w, k, config);

    // Create a directory in temp for 'name'
    std::filesystem::create_directory("temp/" + config.name);

    // Save the lists to files
    save_vector_to_file(gc, "temp/" + config.name + "/specific_gc_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    save_vector_to_file(nongc, "temp/" + config.name + "/specific_nongc_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    save_vector_to_file(ans, "temp/" + config.name + "/specific_ans_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    save_vector_to_file(string_id_to_gc_prefix, "temp/" + config.name + "/specific_string_id_to_gc_prefix_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    save_vector_to_file(string_id_to_gc_suffix, "temp/" + config.name + "/specific_string_id_to_gc_suffix_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");


    // Save vectors of vectors
    save_vector_of_vectors_to_file(kmer_to_gc_string_id, "temp/" + config.name + "/specific_kmer_to_gc_string_id_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    save_vector_of_vectors_to_file(kmer_to_non_gc_string_id, "temp/" + config.name + "/specific_kmer_to_non_gc_string_id_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    save_vector_of_vectors_to_file(string_id_to_non_gc_kmers, "temp/" + config.name + "/specific_string_id_to_non_gc_kmers_" + std::to_string(w) + "_" + std::to_string(k) + ".bin");
    


    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_duration = end - start;

    // Print the time taken
    std::cout << "Preprocessing time for W: " << w << " and K: " << k << " is " << total_duration.count() << " seconds" << std::endl;

    // Save the time to a log file
    std::ofstream log_file("logs/preprocessing_time/" + config.name + "_" + std::to_string(w) + "_" + std::to_string(k) + ".txt");
    log_file << total_duration.count(); // in seconds
    log_file.close();
}


// TODO: some code is duplicate, we can probably store kmer_to_string_id instead of to gc and non_gc
// updates the vectors after a new kmer is chosen
void update_specific_vectors(std::vector<uint64_t>& gc,
    std::vector<uint64_t>& nongc,
    std::vector<std::vector<uint64_t>>& kmer_to_gc_string_id,
    std::vector<std::vector<uint64_t>>& kmer_to_non_gc_string_id,
    std::vector<std::vector<uint64_t>>& string_id_to_non_gc_kmers,
    std::vector<uint64_t>& string_id_to_gc_prefix,
    std::vector<uint64_t>& string_id_to_gc_suffix,
    boost::dynamic_bitset<>& dead_windows,
    uint64_t nextbest) {
    

    // handle windows in which it was a GC kmer
    for (uint64_t string_id : kmer_to_gc_string_id[nextbest]) {
        // skip dead windows
        if (dead_windows.test(string_id)) {
            continue;
        }

		for (uint64_t non_gc_kmer : string_id_to_non_gc_kmers[string_id]) {
			nongc[non_gc_kmer]--;
		}



        uint64_t prefix_kmer = string_id_to_gc_prefix[string_id];
        uint64_t suffix_kmer = string_id_to_gc_suffix[string_id];
        // if is 2^k -> no kmer we still have space for it but we ignore the value so it's fine
        gc[prefix_kmer]--;
        gc[suffix_kmer]--;

        // mark window as dead
        dead_windows.set(string_id);
	
    }

    // handle windows in which it was a non_GC kmer
    for (uint64_t string_id : kmer_to_non_gc_string_id[nextbest]) {
        // skip dead windows
        if (dead_windows.test(string_id)) {
			continue;
		}
        for (uint64_t non_gc_kmer : string_id_to_non_gc_kmers[string_id]) {
            nongc[non_gc_kmer]--;
        }
        uint64_t prefix_kmer = string_id_to_gc_prefix[string_id];
        uint64_t suffix_kmer = string_id_to_gc_suffix[string_id];
        gc[prefix_kmer]--;
        gc[suffix_kmer]--;
        dead_windows.set(string_id);
    }

    return;
}

std::vector<uint64_t> load_all_sequences_particular(uint32_t w, uint32_t k, Config& config) {
    std::lock_guard<std::mutex> guard(load_mutex);
    auto [oddNumbers, evenNumbers] = load_odd_even_pair_from_file(config, w, k);
    return oddNumbers;
}


std::pair<std::vector<uint64_t>, uint64_t> ecogreed_specific(uint32_t w, uint32_t k, double err, std::string path, std::string name, std::vector<uint64_t>& allSequences,
    std::vector<uint64_t> gc, std::vector<uint64_t> nongc, std::vector<uint64_t> ans,
    std::vector<std::vector<uint64_t>>kmer_to_gc_string_id, std::vector<std::vector<uint64_t>>kmer_to_non_gc_string_id, std::vector<std::vector<uint64_t>>& string_id_to_non_gc_kmers, std::vector<uint64_t>& string_id_to_gc_prefix,
    std::vector<uint64_t>& string_id_to_gc_suffix) {

    uint64_t n_kmers = 1ULL << k;
    uint64_t n_sequences = allSequences.size();
    boost::dynamic_bitset<> dead_windows(n_sequences);

    for (uint64_t rank = 0; rank < n_kmers; ++rank) {
        uint64_t nextbest = randcand(k, err, gc, nongc, rank);
        if (nextbest == n_kmers) {
            break;
        }
        ans[nextbest] = rank;

        update_specific_vectors(gc, nongc, kmer_to_gc_string_id, kmer_to_non_gc_string_id, string_id_to_non_gc_kmers, string_id_to_gc_prefix, string_id_to_gc_suffix, dead_windows, nextbest);

     }
        uint64_t gc_count = gc_iterative_particular_binary(w, k, ans, allSequences);

    return { ans, gc_count };
}


void wide_gc_report(uint32_t w, uint32_t k, double err) {


    // start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Load the ordering vector from a file
    std::vector<uint64_t> order = load_vector_from_file("output/ordering/" + std::to_string(w) + "_" + std::to_string(k) + "_" + std::to_string(err) + ".bin");

    // Compute the wide_gc report
    std::unordered_map<uint32_t, boost::multiprecision::cpp_int> answers_per_w = wide_gc(400, k, order, true);
        
        //wide_gc_multiprecision(400, k, order); // since w+1

    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the total time taken and save it to a log file
    std::chrono::duration<double> total_duration = end - start;

    std::ofstream log_file("logs/wide_gc_time" + std::to_string(w) + "_" + std::to_string(k) + "_" + std::to_string(err) + ".txt");
    log_file << total_duration.count(); // in seconds
    log_file.close();

    // Open the result file as a CSV file
    std::ofstream result_file("output/densities/" + std::to_string(w) + "_" + std::to_string(k) + "_" + std::to_string(err) + ".csv");

    // Write the header for the CSV
    result_file << "w,gc,density\n";

    using boost::multiprecision::cpp_dec_float_100; // 50 digits precision

    // Extract and sort keys
    std::vector<uint32_t> keys;
    for (const auto& [key, _] : answers_per_w) {
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());

    // Write the results in CSV format with sorted keys
    for (const uint32_t& key : keys) {
        cpp_dec_float_100 density = cpp_dec_float_100(answers_per_w[key]) / cpp_dec_float_100(std::pow(2, key + k));
        result_file << key << "," << answers_per_w[key] << "," << density.str(10) << "\n"; // Save density with 10 decimal places
    }

    // Close the result file
    result_file.close();
}



std::vector<uint64_t> get_random_order(uint32_t k) {
    std::vector<uint64_t> order;
    uint32_t n_kmers = 1ULL << k;
    for (uint64_t i = 0; i < n_kmers; i++) {
        order.push_back(i);
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(order.begin(), order.end(), g);
    return order;
}