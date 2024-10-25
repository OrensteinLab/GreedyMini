#include <vector>
#include <cstdint>
#include <iostream>
#include <chrono>
#include <bitset>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>
#include "tools.h"
#include "density_and_gc.h"
#include "window_types_manager.h"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <sstream>  // Include this for istringstream
#include <iostream>
#include <fstream>
#include "functions.h"
#include "config.h"

using namespace boost::multiprecision;

// FUNCTIONS FOR CALCULATING DEWSITY FOR BINARY - NORMAL

uint64_t prob_gc(uint32_t w, uint32_t k, const std::vector<uint64_t>& order) {
    uint64_t total = 0;
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;

    for (uint64_t pref = 0; pref < n_kmers; ++pref) {
        uint64_t nextkmer = (pref << 1) & n_kmers_mask;
        total += gc_dfs(w, k, order, nextkmer, order[pref], 1) +
            gc_dfs(w, k, order, nextkmer + 1, order[pref], 1);
    }

    return total;
}

uint64_t gc_dfs(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev) {
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;
    return gc_dfs_prev_could_be_min(w, n_kmers_mask, order, kmer, minorder, lev);
}

uint64_t gc_dfs_prev_could_be_min(uint32_t w, uint64_t n_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev) {
    if (lev == w) {
        return 1;
    }
    else {
        uint64_t nextkm = (kmer << 1) & n_kmers_mask;
        if (order[kmer] >= minorder) {
            return gc_dfs_prev_could_be_min(w, n_kmers_mask, order, nextkm, minorder, lev + 1) +
                gc_dfs_prev_could_be_min(w, n_kmers_mask, order, nextkm + 1, minorder, lev + 1);
        }
        else {
            return gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm, order[kmer], lev + 1) +
                gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm + 1, order[kmer], lev + 1);
        }
    }
}

uint64_t gc_dfs_pref_is_not_min(uint32_t w, uint64_t n_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev) {
    if (lev == w) {
        return order[kmer] < minorder;
    }
    else {
        uint64_t nextkm = (kmer << 1) & n_kmers_mask;
        if (order[kmer] >= minorder) {
            return gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm, minorder, lev + 1) +
                gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm + 1, minorder, lev + 1);
        }
        else {
            return gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm, order[kmer], lev + 1) +
                gc_dfs_pref_is_not_min(w, n_kmers_mask, order, nextkm + 1, order[kmer], lev + 1);
        }
    }
}

double density_expected_binary(uint32_t w, uint32_t k, const std::vector<uint64_t>& order) {
	uint64_t total = prob_gc(w, k, order);
	double density = calc_density(total, k, w);
	return density;
}

// FUNCTIONS FOR CALCULATING DENSITY FOR BINARY - WIDE

std::unordered_map<uint32_t, boost::multiprecision::cpp_int> wide_gc(uint32_t maxw, uint32_t k, const std::vector<uint64_t>& order, bool verbose) {
    uint64_t n_kmers = 1ULL << k;
    std::unordered_map<uint32_t, boost::multiprecision::cpp_int> answers_per_w;

    uint64_t size = n_kmers - std::count(order.begin(), order.end(), n_kmers);
    std::vector<uint64_t> adjusted_order(order.size());
    std::transform(order.begin(), order.end(), adjusted_order.begin(), [size](uint64_t i) { return i < size ? i : size; });

    // prefix gc count for each prefix, whereas the minimum order in the entire window is a specified value
    std::vector<std::vector<boost::multiprecision::cpp_int>> prefix_gc_count_0(n_kmers, std::vector<boost::multiprecision::cpp_int>(size + 1, 0));
    // not gc count for each prefix, whereas the minimum order in the entire window is a specified value
    std::vector<std::vector<boost::multiprecision::cpp_int>> not_gc_count_0(n_kmers, std::vector<boost::multiprecision::cpp_int>(size + 1, 0));
    // suffix gc count for each suffix
    std::vector<boost::multiprecision::cpp_int> suffix_gc_count_0(n_kmers, 0);

    for (uint64_t suffix = 0; suffix < n_kmers; ++suffix) {
        prefix_gc_count_0[suffix][adjusted_order[suffix]] = 1;
    }

    std::vector<std::vector<boost::multiprecision::cpp_int>> prefix_gc_count_1(n_kmers, std::vector<boost::multiprecision::cpp_int>(size + 1, 0));
    std::vector<std::vector<boost::multiprecision::cpp_int>> not_gc_count_1(n_kmers, std::vector<boost::multiprecision::cpp_int>(size + 1, 0));
    std::vector<boost::multiprecision::cpp_int> suffix_gc_count_1 = suffix_gc_count_0;

    for (uint32_t w = 0; w < maxw; ++w) {
        if (verbose) {
            std::cout << "Calculating for W = " << w + 1 << std::endl;
        }
        for (uint64_t suffix = 0; suffix < n_kmers; ++suffix) {
            suffix_gc_count_1[suffix] = 0;

            uint64_t suff1 = suffix >> 1;
            uint64_t suff2 = suff1 + (1ULL << (k - 1));

            for (uint64_t min_order_in_window = 0; min_order_in_window <= adjusted_order[suffix]; ++min_order_in_window) {
                prefix_gc_count_1[suffix][min_order_in_window] = prefix_gc_count_0[suff1][min_order_in_window] + prefix_gc_count_0[suff2][min_order_in_window];
                not_gc_count_1[suffix][min_order_in_window] = not_gc_count_0[suff1][min_order_in_window] + not_gc_count_0[suff2][min_order_in_window];
            }

            for (uint64_t min_order_in_window = adjusted_order[suffix] + 1; min_order_in_window <= size; ++min_order_in_window) {
                suffix_gc_count_1[suffix] += prefix_gc_count_0[suff1][min_order_in_window] + not_gc_count_0[suff1][min_order_in_window] +
                    prefix_gc_count_0[suff2][min_order_in_window] + not_gc_count_0[suff2][min_order_in_window];
                prefix_gc_count_1[suffix][min_order_in_window] = 0;
                not_gc_count_1[suffix][min_order_in_window] = 0;
            }

            if (adjusted_order[suffix] < adjusted_order[suff1]) {
                suffix_gc_count_1[suffix] += suffix_gc_count_0[suff1];
            }
            else {
                not_gc_count_1[suffix][adjusted_order[suff1]] += suffix_gc_count_0[suff1];
            }

            if (adjusted_order[suffix] < adjusted_order[suff2]) {
                suffix_gc_count_1[suffix] += suffix_gc_count_0[suff2];
            }
            else {
                not_gc_count_1[suffix][adjusted_order[suff2]] += suffix_gc_count_0[suff2];
            }
        }

        boost::multiprecision::cpp_int ans = std::accumulate(suffix_gc_count_1.begin(), suffix_gc_count_1.end(), boost::multiprecision::cpp_int(0)) +
            std::accumulate(prefix_gc_count_1.begin(), prefix_gc_count_1.end(), boost::multiprecision::cpp_int(0),
                [](boost::multiprecision::cpp_int sum, const std::vector<boost::multiprecision::cpp_int>& vec) { return sum + std::accumulate(vec.begin(), vec.end(), boost::multiprecision::cpp_int(0)); });

        answers_per_w[w + 1] = ans;

        std::swap(prefix_gc_count_0, prefix_gc_count_1);
        std::swap(suffix_gc_count_0, suffix_gc_count_1);
        std::swap(not_gc_count_0, not_gc_count_1);
    }

    return answers_per_w;
}

std::unordered_map<uint32_t, double> density_expected_binary_wide(uint32_t maxw, uint32_t k, const std::vector<uint64_t>& order, bool verbose) {
	std::unordered_map<uint32_t, boost::multiprecision::cpp_int> answers_per_w = wide_gc(maxw, k, order, verbose);
	std::unordered_map<uint32_t, double> densities_per_w;

	for (auto const& x : answers_per_w) {
		double density = calc_density(x.second, k, x.first);
		densities_per_w[x.first] = density;
	}

	return densities_per_w;
}

// FUNCTIONS FOR CALCULATING DENSITY FOR EXTENDED ORDERS

uint64_t get_extended_order(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint32_t extended_k, uint64_t extended_kmer) {
    uint64_t first_k_bits = extended_kmer & ((1ULL << k) - 1);
    uint64_t last_bits = extended_kmer >> k;
    uint64_t order_of_first_k_bits = order[first_k_bits];
    uint64_t diff_in_k = extended_k - order_of_first_k_bits;
    // shift order_of_first_k_bits by k
    uint64_t order_of_first_k_bits_shifted = order_of_first_k_bits << k;
    // add last_bits to order_of_first_k_bits_shifted
    uint64_t new_order = order_of_first_k_bits_shifted + last_bits;
    return new_order;
}

uint64_t prob_gc_extend_k(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint32_t extended_k) {
    uint64_t total = 0;
    uint64_t n_extended_kmers = 1ULL << extended_k;
    uint64_t n_extended_kmers_mask = n_extended_kmers - 1;

    for (uint64_t pref = 0; pref < n_extended_kmers; ++pref) {
        uint64_t nextkmer = (pref << 1) & n_extended_kmers_mask;
        uint64_t order_pref = get_extended_order(w, k, order, extended_k, pref);
        total += gc_dfs_extend_k(w, k, order, nextkmer, order_pref, 1, extended_k) +
            gc_dfs_extend_k(w, k, order, nextkmer + 1, order_pref, 1, extended_k);
    }

    return total;
}

uint64_t gc_dfs_extend_k(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev, uint32_t extended_k) {
    uint64_t n_extended_kmers = 1ULL << extended_k;
    uint64_t n_extended_kmers_mask = n_extended_kmers - 1;
    return gc_dfs_extend_k_prev_could_be_min(w, k, n_extended_kmers_mask, order, kmer, minorder, lev, extended_k);
}

uint64_t gc_dfs_extend_k_prev_could_be_min(uint32_t w, uint32_t k, uint64_t n_extended_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev, uint32_t extended_k) {
    if (lev == w) {
        return 1;
    }
    else {
        uint64_t order_kmer = get_extended_order(w, k, order, extended_k, kmer);
        uint64_t nextkm = (kmer << 1) & n_extended_kmers_mask;

        if (order_kmer >= minorder) {
            return gc_dfs_extend_k_prev_could_be_min(w, k, n_extended_kmers_mask, order, nextkm, minorder, lev + 1, extended_k) +
                gc_dfs_extend_k_prev_could_be_min(w, k, n_extended_kmers_mask, order, nextkm + 1, minorder, lev + 1, extended_k);
        }
        else {
            return gc_dfs_extend_k_pref_is_not_min(w, k, n_extended_kmers_mask, order, nextkm, order_kmer, lev + 1, extended_k) +
                gc_dfs_extend_k_pref_is_not_min(w, k, n_extended_kmers_mask, order, nextkm + 1, order_kmer, lev + 1, extended_k);
        }
    }
}

uint64_t gc_dfs_extend_k_pref_is_not_min(uint32_t w, uint32_t k, uint64_t n_extended_kmers_mask, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t minorder, uint32_t lev, uint32_t extended_k) {
    uint64_t order_kmer = get_extended_order(w, k, order, extended_k, kmer);
    if (lev == w) {
        return order_kmer < minorder;
    }
    else {
        uint64_t nextkm = (kmer << 1) & n_extended_kmers_mask;
        if (order_kmer >= minorder) {
            return gc_dfs_extend_k_pref_is_not_min(w, k, n_extended_kmers_mask, order, nextkm, minorder, lev + 1, extended_k) +
                gc_dfs_extend_k_pref_is_not_min(w, k, n_extended_kmers_mask, order, nextkm + 1, minorder, lev + 1, extended_k);
        }
        else {
            return gc_dfs_extend_k_pref_is_not_min(w, k, n_extended_kmers_mask, order, nextkm, order_kmer, lev + 1, extended_k) +
                gc_dfs_extend_k_pref_is_not_min(w, k, n_extended_kmers_mask, order, nextkm + 1, order_kmer, lev + 1, extended_k);
        }
    }
}

double density_expected_binary_extend_k(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint32_t extended_k) {
	uint64_t total = prob_gc_extend_k(w, k, order, extended_k);
	double density = calc_density(total, w, extended_k);
	return density;
}

// FUCNTIONS FOR EXPECTED DENSITY ON DNA - NORMAL

std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<uint64_t>> get_gc_ruiners(uint32_t w, uint32_t k) {
    int count_checked = 100000000;
    uint64_t total_sequences = 1ULL << (w + k);


    // Initialize vectors to store the results
    std::vector<uint64_t> prefix_ruiners(total_sequences, 0);
    std::vector<uint64_t> suffix_ruiners(total_sequences, 0);
    std::vector<uint64_t> both_ruiners(total_sequences, 0);

    if (total_sequences <= count_checked) {
        // If 2^(w+k) <= 100,000,000, use all possible binary sequences of length w+k
        for (uint64_t sequence = 0; sequence < total_sequences; ++sequence) {
            uint64_t prefix_ruined_indices = 0;
            uint64_t suffix_ruined_indices = 0;

            // Extract the first and last kmers
            uint64_t first_kmer = sequence >> w;
            uint64_t last_kmer = sequence & ((1ULL << k) - 1);

            // Iterate over each window position in the sequence starting from the second position and without the last position
            for (uint32_t j = 1; j < w; ++j) {
                uint64_t kmer_at_j = (sequence >> (w - j)) & ((1ULL << k) - 1);

                // Compare kmer_at_j with first_kmer lexicographically
                if (kmer_at_j < first_kmer) {
                    prefix_ruined_indices |= (1ULL << j);
                }
                // Compare kmer_at_j with last_kmer lexicographically
                if (kmer_at_j <= last_kmer) {
                    suffix_ruined_indices |= (1ULL << j);
                }
            }

            // Store results in respective vectors
            prefix_ruiners[sequence] = prefix_ruined_indices;
            suffix_ruiners[sequence] = suffix_ruined_indices;
            // Both is the intersection of prefix and suffix
            both_ruiners[sequence] = prefix_ruined_indices & suffix_ruined_indices;
        }
    }
    else {
        // Random number generator for creating random uint64_t sequences
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint64_t> dis(0, total_sequences - 1);

        // Create count_checked random uint64_t sequences
        for (int i = 0; i < count_checked; ++i) {
            uint64_t prefix_ruined_indices = 0;
            uint64_t suffix_ruined_indices = 0;
            uint64_t sequence = dis(gen);

            // Extract the first kmer
            uint64_t first_kmer = sequence >> w;
            uint64_t last_kmer = sequence & ((1ULL << k) - 1);

            // Iterate over each window position in the sequence starting from the second position and without the last position
            for (uint32_t j = 1; j < w; ++j) {
                uint64_t kmer_at_j = (sequence >> (w - j)) & ((1ULL << k) - 1);

                // Compare kmer_at_j with first_kmer lexicographically
                if (kmer_at_j < first_kmer) {
                    prefix_ruined_indices |= (1ULL << j);
                }
                // Compare kmer_at_j with last_kmer lexicographically
                if (kmer_at_j < last_kmer) {
                    suffix_ruined_indices |= (1ULL << j);
                }
            }

            // Store results in respective vectors
            prefix_ruiners[i] = prefix_ruined_indices;
            suffix_ruiners[i] = suffix_ruined_indices;
            // Both is the intersection of prefix and suffix
            both_ruiners[i] = prefix_ruined_indices & suffix_ruined_indices;
        }
    }

    // Return the tuple of vectors
    return { prefix_ruiners, suffix_ruiners, both_ruiners };
}

double get_proportion_gc_per_window_type(uint64_t w, uint64_t k, std::vector<uint64_t>& prefix_ruiners, std::vector<uint64_t>& suffix_ruiners, std::vector<uint64_t>& both_ruiners, uint64_t window_type) {
    // check if window has 1's in first bit
    bool has_first_bit = (window_type & 1) == 1;

    // check if window has 1's in the w'th bit
    bool has_last_bit = (window_type & (1ULL << (w))) == (1ULL << (w));

    // print both bools
    //std::cout << "First bit: " << has_first_bit << ", Last bit: " << has_last_bit << std::endl;


    // if both bits are 1's:
    if (has_first_bit && has_last_bit) {
		uint64_t not_ruined_count = 0;
        // go over all both_ruiners and add 1 if their bitwise AND is at least 1
        for (uint64_t i = 0; i < both_ruiners.size(); ++i) {
			if ((both_ruiners[i] & window_type) == 0) {
                not_ruined_count += 1;
			}
		}
        return static_cast<double>(not_ruined_count) / static_cast<double>(both_ruiners.size());
	}
	// if only the first bit is 1:
    else if (has_first_bit) {
        uint64_t not_ruined_count = 0;
        // go over all prefix_ruiners and add 1 if their bitwise AND is at least 1
        for (uint64_t i = 0; i < prefix_ruiners.size(); ++i) {
            if ((prefix_ruiners[i] & window_type) == 0) {
                not_ruined_count += 1;
            }
        }
        return static_cast<double>(not_ruined_count) / static_cast<double>(prefix_ruiners.size());
    }
    else {
        uint64_t not_ruined_count = 0;
		// go over all suffix_ruiners and add 1 if their bitwise AND is at least 1
		for (uint64_t i = 0; i < suffix_ruiners.size(); ++i) {
            if ((suffix_ruiners[i] & window_type) == 0) {
                not_ruined_count += 1;
			}
		}
		return static_cast<double>(not_ruined_count) / static_cast<double>(suffix_ruiners.size());
	}
}

double get_special_case_expected_gc(uint32_t w, uint32_t k, WindowTypesManager window_types_manager) {
    // generate prefix/suffix/both ruiners
    auto [prefix_ruiners, suffix_ruiners, both_ruiners] = get_gc_ruiners(w, k);

    double expected_gc = 0;

    // go over all window types
    for (auto const& x : window_types_manager.vector_counts) {
		expected_gc += get_proportion_gc_per_window_type(w, k, prefix_ruiners, suffix_ruiners, both_ruiners, x.first) * x.second;
	}
    
	return expected_gc;
}

uint64_t gc_dfs_dna(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t kmer, uint64_t min_order, uint64_t min_indices, uint32_t lev, WindowTypesManager& window_types_manager, bool is_first_min, bool is_middle_min) {
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;

    
    if (lev == w) {
        bool is_last_min = false;
        if (order[kmer] < min_order) {
            return 1; // in case suffix gc
        }
        else if (order[kmer] == min_order) {
            min_indices += (1ULL << lev); // might be an interesting case
            is_last_min = true;
        }


        if ((is_first_min || is_last_min)) {
            // check if interesting case or not
            if (is_middle_min) {
                window_types_manager.add_vector(min_indices);
                return 0;
            }
            else {
                return 1; // suffix/prefix gc
            }
        }
		else {
			return 0;
		}
    }
    else {

        // Here if not on leaves
        uint64_t next_kmer = (kmer << 1) & n_kmers_mask;

        if (order[kmer] < min_order) {
            min_indices = (1ULL << lev);
            return gc_dfs_dna(w, k, order, next_kmer, order[kmer], min_indices, lev + 1, window_types_manager, false, true) +
				gc_dfs_dna(w, k, order, next_kmer + 1, order[kmer], min_indices, lev + 1, window_types_manager, false, true);

        }
        else if (order[kmer] == min_order) {
            min_indices += (1ULL << lev);
            return gc_dfs_dna(w, k, order, next_kmer, min_order, min_indices, lev + 1, window_types_manager, is_first_min, true) +
                gc_dfs_dna(w, k, order, next_kmer + 1, min_order, min_indices, lev + 1, window_types_manager, is_first_min, true);
        }
		else {
			return gc_dfs_dna(w, k, order, next_kmer, min_order, min_indices, lev + 1, window_types_manager, is_first_min, is_middle_min) +
				gc_dfs_dna(w, k, order, next_kmer + 1, min_order, min_indices, lev + 1, window_types_manager, is_first_min, is_middle_min);
		}

    }

}

uint64_t prob_gc_dna(uint32_t w, uint32_t k, const std::vector<uint64_t>& order) {
    WindowTypesManager window_types_manager = WindowTypesManager();
    uint64_t total = 0;
    uint64_t n_kmers = 1ULL << k;
    uint64_t n_kmers_mask = n_kmers - 1;

    for (uint64_t pref = 0; pref < n_kmers; ++pref) {
        uint64_t next_kmer = (pref << 1) & n_kmers_mask;
        total += gc_dfs_dna(w, k, order, next_kmer, order[pref],1, 1, window_types_manager, true, false) +
            gc_dfs_dna(w, k, order, next_kmer + 1, order[pref], 1, 1, window_types_manager, true, false);
    }

    // This counts number of game-changers that are gc regardless of 2nd order
    total *= 1ULL << (w + k);

    // This approximates/counts the number of gamechangers that might be gc/no based on 2nd order
    double special_case_expected_gc = get_special_case_expected_gc(w, k, window_types_manager);
    special_case_expected_gc *= 1ULL << (w + k);

    // add to total
    total += static_cast<uint64_t>(special_case_expected_gc);

    return total;
}

double density_expected_dna(uint32_t w, uint32_t k, const std::vector<uint64_t>& order) {
	uint64_t total = prob_gc_dna(w, k, order);
	double density = calc_density(total, 2*w, 2 * k);
	return density;
}

// FUNCTIONS FOR EXPECTED DENSITY ON DNA - WIDE

std::unordered_map<uint32_t, boost::multiprecision::cpp_int> dna_wide_gc(uint32_t maxw, uint32_t k, const std::vector<uint64_t>& order, bool verbose = false) {
    uint64_t n_kmers = 1ULL << k;
    std::unordered_map<uint32_t, boost::multiprecision::cpp_int> answers_per_w;

    uint64_t size = n_kmers - std::count(order.begin(), order.end(), n_kmers);
    std::vector<uint64_t> adjusted_order(order.size());
    std::transform(order.begin(), order.end(), adjusted_order.begin(), [size](uint64_t i) { return i < size ? i : size; });

    // prefix gc count for each prefix, whereas the minimum order in the entire window is a specified value, the last value is the amount of time minimum order appears
    // NOTE: number of occurrences is stored as -1 of it as we always ahve at least 1
    std::vector<std::vector<std::vector<boost::multiprecision::cpp_int>>> prefix_gc_count_0(
        n_kmers,
        std::vector<std::vector<boost::multiprecision::cpp_int>>(size + 1, std::vector<boost::multiprecision::cpp_int>(maxw + k, 0))
    );
    // not gc count for each prefix, whereas the minimum order in the entire window is a specified value, the last value is the amount of time minimum order appears
    std::vector<std::vector<std::vector<boost::multiprecision::cpp_int>>> not_gc_count_0(
        n_kmers,
        std::vector<std::vector<boost::multiprecision::cpp_int>>(size + 1, std::vector<boost::multiprecision::cpp_int>(maxw + k, 0))
    );
    // suffix gc count for each suffix, 2nd term is the amount of time minimum order appears 
    // TODO: might need to check if now 2 elements are the min order with one being last if to count it as a gc and then divide by 2?
    // Need to count number of occurrences, as now if we have multiple times we can can still have a gc
    std::vector<std::vector<std::vector<boost::multiprecision::cpp_int>>> suffix_gc_count_0(
        n_kmers,
        std::vector<std::vector<boost::multiprecision::cpp_int>>(size + 1, std::vector<boost::multiprecision::cpp_int>(maxw + k, 0))
    );

    std::vector<std::vector<std::vector<boost::multiprecision::cpp_int>>> prefix_and_suffix_gc_count_0(
        n_kmers,
        std::vector<std::vector<boost::multiprecision::cpp_int>>(size + 1, std::vector<boost::multiprecision::cpp_int>(maxw + k, 0))
    );

    for (uint64_t suffix = 0; suffix < n_kmers; ++suffix) {
        prefix_and_suffix_gc_count_0[suffix][adjusted_order[suffix]][0] = 1;
    }

    uint64_t max_possible_order = *std::max_element(adjusted_order.begin(), adjusted_order.end());


    std::vector<std::vector<std::vector<boost::multiprecision::cpp_int>>> prefix_gc_count_1(
        n_kmers,
        std::vector<std::vector<boost::multiprecision::cpp_int>>(size + 1, std::vector<boost::multiprecision::cpp_int>(maxw + k, 0))
    );
    std::vector<std::vector<std::vector<boost::multiprecision::cpp_int>>> not_gc_count_1(
        n_kmers,
        std::vector<std::vector<boost::multiprecision::cpp_int>>(size + 1, std::vector<boost::multiprecision::cpp_int>(maxw + k, 0))
    );
    std::vector<std::vector<std::vector<boost::multiprecision::cpp_int>>> suffix_gc_count_1(
        n_kmers,
        std::vector<std::vector<boost::multiprecision::cpp_int>>(size + 1, std::vector<boost::multiprecision::cpp_int>(maxw + k, 0))
    );
    std::vector<std::vector<std::vector<boost::multiprecision::cpp_int>>> prefix_and_suffix_gc_count_1(
        n_kmers,
        std::vector<std::vector<boost::multiprecision::cpp_int>>(size + 1, std::vector<boost::multiprecision::cpp_int>(maxw + k, 0))
    );




    for (uint32_t w = 1; w <= maxw; ++w) {
        if (verbose) {
            std::cout << "Window size: " << w << std::endl;
        }
        for (uint64_t suffix = 0; suffix < n_kmers; ++suffix) {
            
            // reset the values in the temporary arrays
            for (uint64_t min_order_in_window = 0; min_order_in_window <= max_possible_order; ++min_order_in_window) {
                for (uint64_t min_apperances_in_window = 0; min_apperances_in_window < w+2; ++min_apperances_in_window) {
                    suffix_gc_count_1[suffix][min_order_in_window][min_apperances_in_window] = 0;
                    prefix_gc_count_1[suffix][min_order_in_window][min_apperances_in_window] = 0;
                    not_gc_count_1[suffix][min_order_in_window][min_apperances_in_window] = 0;
                    prefix_and_suffix_gc_count_1[suffix][min_order_in_window][min_apperances_in_window] = 0;
                }
			}

            // generate the preceding kmers for the given suffix
            uint64_t suff1 = suffix >> 1;
            uint64_t suff2 = suff1 + (1ULL << (k - 1));


            // case min order in window is strictly below the order of the suffix -> can only be a prefix GC or a non-gc
            for (uint64_t min_order_in_window = 0; min_order_in_window < adjusted_order[suffix]; ++min_order_in_window) {
                for (uint64_t min_apperances_in_window = 0; min_apperances_in_window < w +1; ++min_apperances_in_window) {
                    prefix_gc_count_1[suffix][min_order_in_window][min_apperances_in_window] += prefix_gc_count_0[suff1][min_order_in_window][min_apperances_in_window] + prefix_gc_count_0[suff2][min_order_in_window][min_apperances_in_window];
                    prefix_gc_count_1[suffix][min_order_in_window][min_apperances_in_window] += prefix_and_suffix_gc_count_0[suff1][min_order_in_window][min_apperances_in_window] + prefix_and_suffix_gc_count_0[suff2][min_order_in_window][min_apperances_in_window];

                    not_gc_count_1[suffix][min_order_in_window][min_apperances_in_window] += not_gc_count_0[suff1][min_order_in_window][min_apperances_in_window] + not_gc_count_0[suff2][min_order_in_window][min_apperances_in_window];
                    not_gc_count_1[suffix][min_order_in_window][min_apperances_in_window] += suffix_gc_count_0[suff1][min_order_in_window][min_apperances_in_window] + suffix_gc_count_0[suff2][min_order_in_window][min_apperances_in_window];
                }
            }
            // case min order is same order as suffix -> will always be a gamechanger, but can be either a suffix or a prefix and suffix gamechanger
            for (uint64_t min_apperances_in_window = 0; min_apperances_in_window < w +1; ++min_apperances_in_window) {
                prefix_and_suffix_gc_count_1[suffix][adjusted_order[suffix]][min_apperances_in_window+1] += prefix_and_suffix_gc_count_0[suff1][adjusted_order[suffix]][min_apperances_in_window] + prefix_and_suffix_gc_count_0[suff2][adjusted_order[suffix]][min_apperances_in_window];
                prefix_and_suffix_gc_count_1[suffix][adjusted_order[suffix]][min_apperances_in_window + 1] += prefix_gc_count_0[suff1][adjusted_order[suffix]][min_apperances_in_window] + prefix_gc_count_0[suff2][adjusted_order[suffix]][min_apperances_in_window];

                suffix_gc_count_1[suffix][adjusted_order[suffix]][min_apperances_in_window + 1] += suffix_gc_count_0[suff1][adjusted_order[suffix]][min_apperances_in_window] + suffix_gc_count_0[suff2][adjusted_order[suffix]][min_apperances_in_window];
                suffix_gc_count_1[suffix][adjusted_order[suffix]][min_apperances_in_window + 1] += not_gc_count_0[suff1][adjusted_order[suffix]][min_apperances_in_window] + not_gc_count_0[suff2][adjusted_order[suffix]][min_apperances_in_window];
            }

            // case min order is larger than the suffix -> can only be a suffix GC
            for (uint64_t min_order_in_window = adjusted_order[suffix] + 1; min_order_in_window <= size; ++min_order_in_window) {
                for (uint64_t min_apperances_in_window = 0; min_apperances_in_window < w + 1 ; ++min_apperances_in_window) {
                    suffix_gc_count_1[suffix][adjusted_order[suffix]][0] += prefix_gc_count_0[suff1][min_order_in_window][min_apperances_in_window] + prefix_gc_count_0[suff2][min_order_in_window][min_apperances_in_window];
                    suffix_gc_count_1[suffix][adjusted_order[suffix]][0] += suffix_gc_count_0[suff1][min_order_in_window][min_apperances_in_window] + suffix_gc_count_0[suff2][min_order_in_window][min_apperances_in_window];
                    suffix_gc_count_1[suffix][adjusted_order[suffix]][0] += prefix_and_suffix_gc_count_0[suff1][min_order_in_window][min_apperances_in_window] + prefix_and_suffix_gc_count_0[suff2][min_order_in_window][min_apperances_in_window];
                    suffix_gc_count_1[suffix][adjusted_order[suffix]][0] += not_gc_count_0[suff1][min_order_in_window][min_apperances_in_window] + not_gc_count_0[suff2][min_order_in_window][min_apperances_in_window];
                }
            }
        }
        using boost::multiprecision::cpp_dec_float_100;

        cpp_dec_float_100 ans = 0;


        for (uint64_t suffix = 0; suffix < n_kmers; ++suffix) {
            for (uint64_t min_order_in_window = 0; min_order_in_window <= size; ++min_order_in_window) {
                for (uint64_t min_apperances_in_window = 0; min_apperances_in_window < w + k; ++min_apperances_in_window) {
                    ans += cpp_dec_float_100(suffix_gc_count_1[suffix][min_order_in_window][min_apperances_in_window]) / cpp_dec_float_100(min_apperances_in_window+1);
                    ans += cpp_dec_float_100(prefix_gc_count_1[suffix][min_order_in_window][min_apperances_in_window]) / cpp_dec_float_100(min_apperances_in_window+1);
                    ans += (2 * cpp_dec_float_100(prefix_and_suffix_gc_count_1[suffix][min_order_in_window][min_apperances_in_window]) / cpp_dec_float_100(min_apperances_in_window+1));
                }
            }
        }
        answers_per_w[w] = boost::multiprecision::cpp_int(ans);

        std::swap(prefix_gc_count_0, prefix_gc_count_1);
        std::swap(suffix_gc_count_0, suffix_gc_count_1);
        std::swap(not_gc_count_0, not_gc_count_1);
        std::swap(prefix_and_suffix_gc_count_0, prefix_and_suffix_gc_count_1);
    }

    return answers_per_w;
}

std::unordered_map<uint32_t, double> density_expected_dna_wide(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, bool verbose) {
	auto answers_per_w = dna_wide_gc(w, k, order, verbose);
    std::unordered_map<uint32_t, double> densities_per_w;
    for (auto const& x : answers_per_w) {
		densities_per_w[x.first] = calc_density(x.second, k, x.first);
    }
	return densities_per_w;
}

// FUNCTIONS FOR PARTICULAR DENSITY ON DNA

uint64_t gc_iterative_particular_dna(uint32_t w, uint32_t k, std::vector<uint64_t> order, std::vector<uint64_t> oddNumbers, std::vector<uint64_t> evenNumbers) {
    uint64_t gc_count = 0;
    uint32_t oddNumbersLength = static_cast<uint32_t>(oddNumbers.size());
    for (size_t i = 0; i < oddNumbersLength; i++)
    {
        // for each odd sequence, we want to check whether or not it is a GC/not GC/ perhaps a GC

        uint64_t oddSequence = oddNumbers[i];

        uint64_t first_kmer = oddSequence >> w;
        uint64_t minOrder = order[first_kmer];

        // Stores the indices of the kmers that have the minimum order in binary, that is 0010100 means that the 3rd and 5th kmers have the minimum order
        uint64_t minIndices = 1;

        bool isFirstMin = true;
        bool isMiddleMin = false;

        // checks for all kmers except the last one
        for (size_t j = 1; j < w; j++)
        {
            uint64_t kmer_at_j = (oddSequence >> (w - j)) & ((1ULL << k) - 1);
            if (order[kmer_at_j] < minOrder) {
                isFirstMin = false;
                minOrder = order[kmer_at_j];
                minIndices = 1ULL << j;
                isMiddleMin = true;
            }
            else if (order[kmer_at_j] == minOrder) {
                minIndices |= 1ULL << j;
                isMiddleMin = true;
			}
        }

        // checks for the last kmer
        bool isLastMin = false;
        uint64_t kmer_at_w = oddSequence & ((1ULL << k) - 1);
        if (order[kmer_at_w] < minOrder) {
            // will always be a suffix GC
			gc_count += 1;
            continue;
		}
		else if (order[kmer_at_w] == minOrder) {
			minIndices |= 1ULL << w;
            isLastMin = true;
		}

        // Finally we handle the cases where it COULD be a GC
        if (isFirstMin || isLastMin) {
            if (!isMiddleMin) {
                gc_count += 1;
                continue;
            }
            else {
                // Here if this is an interesting case so we need to use the even sequence to check if it is a GC
                uint64_t evenSequence = evenNumbers[i];

                uint64_t evenFirstKmer = evenSequence >> w;
                uint64_t evenLastKmer = evenSequence & ((1ULL << k) - 1);
                uint64_t evenFirstOrder = order[evenFirstKmer];
                uint64_t evenLastOrder = order[evenLastKmer];

                bool canBePrefixGC = isFirstMin;
                bool canBeSuffixGC = isLastMin;

                for (size_t l = 1; l < w; l++)
                {
                    // skip indexes that aren't the minimum in the odd sequence
                    if ((minIndices & (1ULL << l)) == 0) {
						continue;
					}
                    // Early stopping
                    if (!canBePrefixGC && !canBeSuffixGC) {
                        break;
                    }


                    uint64_t kmer_at_l = (evenSequence >> (w - l)) & ((1ULL << k) - 1);
                    uint64_t order_at_l = order[kmer_at_l];

                    // check if this ruins prefix/suffix GC
                    if (order_at_l < evenFirstOrder) {
						canBePrefixGC = false;
					}
                    if (order_at_l <= evenLastOrder) {
                        canBeSuffixGC = false;
                    }
                }
                if (canBePrefixGC || canBeSuffixGC) {
					gc_count += 1;
				}
            }
			

        }
			
    }
    return gc_count;
}

uint64_t gc_iterative_particular_binary(uint32_t w, uint32_t k, std::vector<uint64_t> order, std::vector<uint64_t> numbers) {
    uint64_t gc_count = 0;
    uint32_t numbersLength = static_cast<uint32_t>(numbers.size());

    for (size_t i = 0; i < numbersLength; i++)
    {
        uint64_t sequence = numbers[i];
        uint64_t first_kmer = sequence >> w;
        uint64_t minOrder = order[first_kmer];

        bool isFirstMin = true;

        // check for all kmers except the last one
        for (size_t j = 1; j < w; j++)
		{
			uint64_t kmer_at_j = (sequence >> (w - j)) & ((1ULL << k) - 1);
			if (order[kmer_at_j] < minOrder) {
				isFirstMin = false;
				minOrder = order[kmer_at_j];
			}
		}

 

        // Suffix GC
        uint64_t kmer_at_w = sequence & ((1ULL << k) - 1);
        if (order[kmer_at_w] < minOrder) {
			gc_count += 1;
			continue;
		}

        // Prefix GC
        if (isFirstMin) {
            gc_count += 1;
            continue;
        }


    }

    return gc_count;
}

double density_particular_dna(uint32_t W, uint32_t K, std::string path, std::vector<uint64_t> order)
{
    auto [oddNumbers, evenNumbers] = load_odd_even_pair_from_file(path, W,K);
    uint64_t gcCount = gc_iterative_particular_dna(W, K, order, oddNumbers, evenNumbers);
    uint64_t kmers_picked = gcCount + 1;
    uint32_t oddNumbersLength = static_cast<uint32_t>(oddNumbers.size());
    uint32_t kmers_possible = oddNumbersLength + W; // since each context has w+1 kmers;

    double density = calc_density_specific<uint64_t>(kmers_picked, K, W, kmers_possible);
    return density;
}

double density_dna(uint32_t W, uint32_t K, double min_alpha, double max_alpha, bool swapped) {
    std::vector<uint64_t> order = load_order(W, K, min_alpha, max_alpha, swapped);
    uint64_t gc_count = prob_gc_dna(W, K, order);
    double density = calc_density<uint64_t>(gc_count, 2 * W, 2 * K);
    return density;
}


