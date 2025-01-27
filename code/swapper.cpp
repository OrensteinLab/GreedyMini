#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>
#include <chrono>
#include "swapper.h"
#include "functions.h"
//#include "gc_counter.h"
#include "density_and_gc.h"
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <random>
#include <chrono>
#include <algorithm>
#include "tools.h"

uint64_t cost(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t u, uint64_t v, uint64_t y) {
    std::vector<uint64_t> kmer(w + 1, 0);   // kmers in active nodes of DFS
    std::vector<uint8_t> has_v(w + 1, 0);   // flags "has v"
    std::vector<uint8_t> visits(w + 1, 0);  // num of previous visits to a node
    std::vector<uint64_t> curcost(w + 1, 0);// cost over subtree
    int32_t lev = 0;                        // use signed int for level control
    kmer[0] = u;

    while (visits[0] < 2) { // until the last return to the root
        if (visits[lev] == 0) {
            if (order[kmer[lev]] < y) {
                lev--;
                visits[lev]++;
            }
            else {
                if (lev == static_cast<int32_t>(w)) {
                    if (has_v[lev]) {
                        curcost[w - 1]++;
                    }
                    lev--;
                    visits[lev]++;
                }
                else {
                    if (kmer[lev] == v) {
                        has_v[lev] = 1;
                    }
                    lev++;
                    kmer[lev] = (kmer[lev - 1] << 1) % (1ULL << k);
                    visits[lev] = 0;
                    curcost[lev] = 0;
                    has_v[lev] = has_v[lev - 1];
                }
            }
        }
        else if (visits[lev] == 1) {
            lev++;
            kmer[lev]++;
            visits[lev] = 0;
            curcost[lev] = 0;
            has_v[lev] = has_v[lev - 1];
        }
        else {
            lev--;
            visits[lev]++;
            curcost[lev] += curcost[lev + 1];
        }
    }

    uint64_t tempcost = curcost[0];
    // print
    //std::cout << "Original prefix cost: " << tempcost << std::endl;

    visits[0] = 0;
    curcost[0] = 0;
    lev = 0;

    // Same for suffix tree
    while (visits[0] < 2) {
        if (visits[lev] == 0) {
            if (order[kmer[lev]] < y || (lev > 0 && kmer[lev] == u)) {
                lev--;
                visits[lev]++;
            }
            else {
                if (lev == static_cast<int32_t>(w)) {
                    if (has_v[w] && kmer[w] != v) {
                        curcost[w - 1]++;
                    }
                    lev--;
                    visits[lev]++;
                }
                else {
                    if (kmer[lev] == v) {
                        has_v[lev] = 1;
                    }
                    lev++;
                    kmer[lev] = kmer[lev - 1] >> 1;
                    visits[lev] = 0;
                    curcost[lev] = 0;
                    has_v[lev] = has_v[lev - 1];
                }
            }
        }
        else if (visits[lev] == 1) {
            lev++;
            kmer[lev] += 1ULL << (k - 1);
            visits[lev] = 0;
            curcost[lev] = 0;
            has_v[lev] = has_v[lev - 1];
        }
        else {
            lev--;
            visits[lev]++;
            curcost[lev] += curcost[lev + 1];
        }
    }
    // print
    //std::cout << "Original suffix cost: " << curcost[0] << std::endl;
    return tempcost + curcost[0];
}


uint64_t cost_prefix(uint32_t lev, uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t u, uint64_t v, uint64_t y, uint64_t cur_kmer, bool has_v) {
    if ((lev == w)) {
        if (order[cur_kmer] < y) {
			return 0;
		}
        return (has_v ? 1 : 0);
    }

    uint64_t next_kmer = (cur_kmer << 1) % (1ULL << k);
    bool next_has_v = has_v || (cur_kmer == v);

    uint64_t result = 0;
    if (order[cur_kmer] >= y) {
        result += cost_prefix(lev + 1, w, k, order, u, v, y, next_kmer, next_has_v);
        result += cost_prefix(lev + 1, w, k, order, u, v, y, next_kmer + 1, next_has_v);
    }
    return result;
}

uint64_t cost_suffix(uint32_t lev, uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t u, uint64_t v, uint64_t y, uint64_t cur_kmer, bool has_v) {
// cur_kmer==u check is because if cur_kmer ==u it means its a prefix GC and not a suffix
    if (lev == w) {
        if (order[cur_kmer] < y || cur_kmer == u) {
			return 0;
		}
        return ((has_v) && (cur_kmer != v)) ? 1 : 0;
    }

    uint64_t next_kmer1 = cur_kmer >> 1;
    uint64_t next_kmer2 = next_kmer1 + (1ULL << (k - 1));
    bool next_has_v = has_v || (cur_kmer == v);

    uint64_t result = 0;

    if (lev == 0) {
        result += cost_suffix(lev + 1, w, k, order, u, v, y, next_kmer1, next_has_v);
        result += cost_suffix(lev + 1, w, k, order, u, v, y, next_kmer2, next_has_v);
    }else if ((order[cur_kmer] >= y) && (cur_kmer != u)) {
        result += cost_suffix(lev + 1, w, k, order, u, v, y, next_kmer1, next_has_v);
        result += cost_suffix(lev + 1, w, k, order, u, v, y, next_kmer2, next_has_v);
    }

    return result;
}

uint64_t cost_recursive(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t u, uint64_t v, uint64_t y) {
    uint64_t prefix_cost = cost_prefix(0, w, k, order, u, v, y, u, false);
    uint64_t suffix_cost = cost_suffix(0, w, k, order, u, v, y, u, false);

    // print the prefix and suffix cost
    //std::cout << "Prefix cost: " << prefix_cost << std::endl;
    //std::cout << "Suffix cost: " << suffix_cost << std::endl;
    return prefix_cost + suffix_cost;
}


std::pair<std::vector<uint64_t>, uint64_t> swapper_f(uint32_t w, uint32_t k, std::vector<uint64_t>& order, double max_time_seconds, bool verbose) {
    uint64_t size = (1ULL << k) - std::count(order.begin(), order.end(), 1ULL << k);


    uint64_t GCcount = prob_gc(w, k, order);
    std::vector<uint64_t> perm(size);

    for (uint64_t i = 0; i < size; ++i) {
        perm[i] = std::distance(order.begin(), std::find(order.begin(), order.end(), i));
    }

    uint64_t tries = 0;
    std::vector<uint64_t> o1 = order;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<uint64_t> dis(0, size - 2);



    auto start_time = std::chrono::high_resolution_clock::now();

    uint32_t n_swaps = 0;
    uint32_t n_good_swaps = 0;

    while (true) {

        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - start_time;

        if (elapsed.count() >= max_time_seconds) {
            if (verbose) {
                std::cout << "Time limit reached, swapped " << n_swaps << " times out of " << tries << " times (on main thread)" << std::endl;
            }
            break;
        }

        uint64_t y = dis(gen);
        // how many prefix/suffix gamechangers where perm[y] is picked would be ruined if perm[y+1] was smaller rank
        int64_t gap = cost_recursive(w, k, o1, perm[y + 1], perm[y], y) - cost_recursive(w, k, o1, perm[y], perm[y + 1], y);


        if (gap < 0) {  // swap is profitable, always swap
            GCcount += gap;
            std::swap(o1[perm[y]], o1[perm[y + 1]]);
            std::swap(perm[y], perm[y + 1]);
            n_swaps++;
            n_good_swaps++;
        }
        else if (gap == 0 && tries % 2 == 1) {  // swap is neutral, we swap on odd tries and skip even tries
            std::swap(o1[perm[y]], o1[perm[y + 1]]);
            std::swap(perm[y], perm[y + 1]);
            n_swaps++;

        }
        tries++;
    }

    return { o1, GCcount };
}




uint64_t dp_gc_count(uint32_t w, uint32_t k, uint32_t rank, std::vector<uint64_t>& order, std::vector<uint64_t>& table_0, std::vector<uint64_t>& table_1) {

    std::vector<uint64_t> old_vec = std::vector<uint64_t>(1ULL << k, 0);
    std::vector<uint64_t> new_vec = std::vector<uint64_t>(1ULL << k, 0);


    uint64_t max_rank = 1 << k;

    // calculate prefix GC count

    // setuo
    old_vec[table_0[rank]] = 1;
    old_vec[table_1[rank]] = 1;

    // dp
    for (uint32_t iteration = 1; iteration < w; iteration++)
    {
        for (uint64_t following_rank = rank; following_rank < max_rank; following_rank++)
        {
            new_vec[table_0[following_rank]] += old_vec[following_rank];
            new_vec[table_1[following_rank]] += old_vec[following_rank];
            old_vec[following_rank] = 0;
        }
        std::swap(old_vec, new_vec);
    }
    // Calculate the sum of old_vec from index 'rank' to the end
    uint64_t prefix_gc = std::accumulate(old_vec.begin() + rank, old_vec.end(), 0ULL);



    // cleaning up and initializing for suffix GC count
    for (uint64_t previous_rank = 0; previous_rank <= rank; previous_rank++)
    {
        old_vec[previous_rank] = 0;
        new_vec[previous_rank] = 0;
    }
    for (uint64_t following_rank = rank+1; following_rank < max_rank; following_rank++)
    {
        old_vec[following_rank] = 1;
	}

    // dp
    for (uint32_t iteration = 0; iteration < w; iteration++)
    {
        for (uint64_t following_rank = rank+1; following_rank < max_rank; following_rank++)
        {
            new_vec[table_0[following_rank]] += old_vec[following_rank];
            new_vec[table_1[following_rank]] += old_vec[following_rank];
            old_vec[following_rank] = 0;
        }
        old_vec[rank] = 0;
        std::swap(old_vec, new_vec);
    }

    uint64_t suffix_gc = old_vec[rank];

    return prefix_gc + suffix_gc;
}









std::pair<std::vector<uint64_t>, uint64_t> swapper_f_v2(uint32_t w, uint32_t k, std::vector<uint64_t>& order, double max_time_seconds, bool verbose, bool is_first) {
    uint64_t size = (1ULL << k) - std::count(order.begin(), order.end(), 1ULL << k);
    uint64_t n_kmers = (1ULL << k);
    uint64_t kmer_mask = (1ULL << k) - 1;
    uint64_t left_bit = 1ULL << (k - 1);


    // make sure ranks are unique
    uint64_t temp_rank = size;
    for (uint64_t kmer = 0; kmer < n_kmers; kmer++)
    {
        if(order[kmer] == n_kmers)
		{
			order[kmer] = temp_rank;
            temp_rank++;
		}
    }

    // rank to the kmer it belongs to
    std::vector<uint64_t> perm(n_kmers);
    for (uint64_t i = 0; i < n_kmers; i++)
    {
        // Find rank i in order:
        auto it = std::find(order.begin(), order.end(), i);
        if (it == order.end())
        {
            perm[i] = n_kmers;  // or handle error
        }
        else
        {
            perm[i] = static_cast<uint64_t>(it - order.begin());
        }
    }



    // initialize the tables with zeros
    std::vector<uint64_t> table_0(n_kmers, 0);
    std::vector<uint64_t> table_1(n_kmers, 0);


    // fill the tables
    for (uint64_t rank = 0; rank < n_kmers; rank++)
    {
        uint64_t next_kmer = (perm[rank] << 1) & kmer_mask;
        table_0[rank] = order[next_kmer];
        table_1[rank] = order[next_kmer + 1];
    }


    uint64_t gc_count = 0;
    for (uint64_t i = 0; i < size; i++)
    {
        gc_count += dp_gc_count(w, k, i, order, table_0, table_1);
    }

    // print gc count, TODO: remove
    //if (is_first) {
    //    std::cout << "GC count: " << gc_count << std::endl;
    //    double density = calc_density(gc_count, k, w);
    //    std::cout << "Density: " << density << std::endl;
    //}


    uint64_t tries = 0;
    uint64_t n_swaps = 0;
    uint64_t n_good_swaps = 0;
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<uint64_t> dis(0, size - 2);



    auto start_time = std::chrono::high_resolution_clock::now();

    while (true) {
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - start_time;

        if (elapsed.count() >= max_time_seconds) {
            if (verbose) {
                std::cout << "Time limit reached, swapped " << n_swaps << " times out of " << tries << " times (on main thread)" << std::endl;
            }
            break;
        }



        uint64_t y = dis(gen);

        // count how many charged contexts we have right now which are charged by rank y or y+1 kmers
        uint64_t charged_before_swap = dp_gc_count(w, k, y, order, table_0, table_1) + dp_gc_count(w, k, y + 1, order, table_0, table_1);


        // swap the kmers with ranks y and y+1
        uint64_t a = perm[y];
        uint64_t b = perm[y + 1];

        uint64_t c = (a >> 1);
        uint64_t d = (b >> 1);

        if (a & 1) {
            table_1[order[c]] = y + 1;
            table_1[order[c + left_bit]] = y + 1;
        }
        else {
            table_0[order[c]] = y + 1;
            table_0[order[c + left_bit]] = y + 1;
        }

        if (b & 1) {
            table_1[order[d]] = y;
            table_1[order[d + left_bit]] = y;
        }
        else {
            table_0[order[d]] = y;
            table_0[order[d + left_bit]] = y;
        }

        std::swap(table_0[y], table_0[y + 1]);
        std::swap(table_1[y], table_1[y + 1]);

        // charged_before_swap
        uint64_t charged_after_swap = dp_gc_count(w, k, y, order, table_0, table_1) + dp_gc_count(w, k, y + 1, order, table_0, table_1);



        // if profitable or 50% when neutral
        if (charged_before_swap > charged_after_swap) {
            uint64_t gap = charged_before_swap - charged_after_swap;
            gc_count -= gap;
            //std::cout << "GC count: " << gc_count << std::endl;
            std::swap(order[perm[y]], order[perm[y + 1]]);
            std::swap(perm[y], perm[y + 1]);
            n_swaps++;
            n_good_swaps++;
        }
        else if ((charged_before_swap == charged_after_swap && tries % 2 == 1)) {
            std::swap(order[perm[y]], order[perm[y + 1]]);
            std::swap(perm[y], perm[y + 1]);
            n_swaps++;
        }
        else {
            // undo the changes
            std::swap(table_0[y], table_0[y + 1]);
            std::swap(table_1[y], table_1[y + 1]);

            if (a & 1) {
				table_1[order[c]] = y;
				table_1[order[c + left_bit]] = y;
			}
			else {
				table_0[order[c]] = y;
				table_0[order[c + left_bit]] = y;
			}

            if (b & 1) {
                table_1[order[d]] = y + 1;
                table_1[order[d + left_bit]] = y + 1;
            }
            else {
                table_0[order[d]] = y + 1;
				table_0[order[d + left_bit]] = y + 1;
            }
        }
        tries++;
    }

    // restore the UHS form of the order
    for (uint64_t kmer = 0; kmer < n_kmers; kmer++)
    {
        if (order[kmer] >= size) {
            order[kmer] = 1ULL << k;
        }
    }

    return { order, gc_count };




}


