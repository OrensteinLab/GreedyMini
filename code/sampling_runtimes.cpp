#include <iostream>
#include <chrono>
#include <vector>
#include <thread>
#include <mutex>
#include <limits>
#include "tools.h"
#include <algorithm> 
#include <cmath>
#include <random>     
#include <map>
#include <iomanip> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "config.h"
#include <iostream>
//#include <cstdlib>  
#include <ctime>    
//#include <cstdlib>  
#include <stdexcept>
#include <unordered_set>
#include <queue>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <utility> 
#include "sampling_runtimes.h"



std::vector<bool> get_bits(uint64_t n_bits) {
	std::vector<bool> bits(n_bits);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0, 1);

	for (size_t i = 0; i < n_bits; ++i) {
		bits[i] = dist(gen); // Either 0 or 1
	}

	return bits;

}






double test_runtime_xor_hash(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {
	uint64_t k_mask = (1 << k) - 1;
	uint64_t dna_k_mask = (1 << (2 * k)) - 1;

	uint64_t random_number = 5215;


	uint64_t sum = 0;

	// start timing
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t current_kmer = 0;
	for (uint64_t i = 0; i < N; i += 2) {
		current_kmer = (current_kmer << 1) & dna_k_mask;
		current_kmer |= (bits[i] << 1) | bits[i + 1];
		uint64_t kmer_rank = current_kmer ^ random_number;
		sum += kmer_rank;
	}



	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken xor hash (ms):  " << elapsed_us / 1000.0 << "\n";

	return elapsed_us / 1000.0;

}




double test_runtime_gm(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {

	uint64_t k_mask = (1 << k) - 1;

	// 1. Create a vector "order" of size 2^k
	std::vector<std::uint64_t> order(1ULL << k);

	// 2. Set up random number generation
	std::random_device rd;              // Non-deterministic seed
	std::mt19937_64 gen(rd());          // Mersenne Twister 64-bit
	std::uniform_int_distribution<std::uint64_t> dist(
		0, std::numeric_limits<std::uint64_t>::max()
	);

	// 3. Fill the vector with random data
	for (auto& elem : order) {
		elem = dist(gen);
	}




	//std::vector<uint64_t> order = load_gm_order(w, k);

	// multiply all ranks in order by 2^k
	for (uint64_t i = 0; i < order.size(); i++) {
		order[i] = order[i] << k;
	}


	uint64_t sum = 0;




	// start timing 
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t odd_kmer = 0;
	uint64_t even_kmer = 0;

	for (uint64_t i = 0; i < N; i += 2) {
		even_kmer = (even_kmer << 1) & k_mask;
		even_kmer |= bits[i];
		odd_kmer = (odd_kmer << 1) & k_mask;
		odd_kmer |= bits[i + 1];
		uint64_t kmer_rank = order[odd_kmer] + even_kmer;
		sum += kmer_rank;
	}


	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken gm (ms):  " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}


double test_runtime_gm_no_memory_access(uint32_t w, uint32_t k, uint64_t N, const std::vector<bool>& bits) {

	uint64_t k_mask = (1 << k) - 1;
	uint64_t sum = 0;




	// start timing 
	auto start = std::chrono::high_resolution_clock::now();

	uint64_t odd_kmer = 0;
	uint64_t even_kmer = 0;

	for (uint64_t i = 0; i < N; i += 2) {
		even_kmer = (even_kmer << 1) & k_mask;
		even_kmer |= bits[i];
		odd_kmer = (odd_kmer << 1) & k_mask;
		odd_kmer |= bits[i + 1];
		uint64_t kmer_rank = odd_kmer + even_kmer;
		sum += kmer_rank;
	}


	auto end = std::chrono::high_resolution_clock::now();

	std::cout << "sum:\t" << sum << std::endl;

	// Measure elapsed time in nanoseconds
	auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	// Or microseconds
	auto elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

	std::cout << "Time taken gm no memory (ms):  " << elapsed_us / 1000.0 << "\n";
	return elapsed_us / 1000.0;
}



void perform_all_sampling_tests() {
	const uint64_t N = 6'000'000'000; // 3 billion nuc
	uint32_t w = 10;

	std::vector<bool> bits = get_bits(N);



	// Choose an output file name
	std::string filename = "sampling_runtime.csv";

	// Open file in write (trunc) mode initially. If you want to append,
	// change std::ios_base::out | std::ios_base::trunc to std::ios_base::app
	std::ofstream ofs(filename, std::ios_base::out | std::ios_base::trunc);
	if (!ofs.is_open()) {
		std::cerr << "Error: Failed to open file: " << filename << std::endl;
		return;
	}

	// Write the CSV header
	ofs << "k,time_random,time_gm_no_mem,time_gm\n";


	for (uint32_t k = 3; k <= 31; ++k)
	{
		std::cout << "\n\n\nRunning tests for k = " << k << std::endl;
		double time_random = test_runtime_xor_hash(w, k, N, bits);
		double time_gm_no_mem = test_runtime_gm_no_memory_access(w, k, N, bits);
		double time_gm = test_runtime_gm(w, k, N, bits);


		// Write a single row of CSV output
		ofs << k << ","
			<< time_random << ","
			<< time_gm_no_mem << ","
			<< time_gm << "\n";
	}

	// Close the file
	ofs.close();

	std::cout << "Results have been saved to " << filename << std::endl;

}