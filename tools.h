#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <string>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdint> 
#include <iostream>
#include <fstream>
#include <map>
#include "config.h"


using namespace boost::multiprecision;


void save_vector_to_file(const std::vector<uint64_t>& vec, const std::string& filename);
std::vector<uint64_t> load_vector_from_file(const std::string& filename);
bool file_exists(const std::string& filename);


template <typename GCType>
double calc_density_specific(GCType kmers_picked, uint32_t k, uint32_t w, uint32_t n_possible_kmers) {
    using boost::multiprecision::cpp_dec_float_100;

    cpp_dec_float_100 denominator = cpp_dec_float_100(n_possible_kmers);

    // Calculate the density as GC_count / (2^(w + k))
    cpp_dec_float_100 density = cpp_dec_float_100(kmers_picked) / denominator;

    // Convert density to double and return
    return density.convert_to<double>();
}

template <typename GCType>
double calc_density(GCType GC_count, uint32_t k, uint32_t w) {
    using boost::multiprecision::cpp_dec_float_100;

    // Calculate 2^(w + k)
    cpp_dec_float_100 denominator = cpp_dec_float_100(std::pow(2, w + k));

    // Calculate the density as GC_count / (2^(w + k))
    cpp_dec_float_100 density = cpp_dec_float_100(GC_count) / denominator;

    // Convert density to double and return
    return density.convert_to<double>();
}


std::vector<uint64_t> get_reversed_order(const std::vector<uint64_t>& order);
std::vector<std::string> split(const std::string& str, char delimiter);
std::vector<uint64_t> get_explicitly_extended_order(const std::vector<uint64_t>& order);
void save_vector_of_vectors_to_file(const std::vector<std::vector<uint64_t>>& vec_of_vecs, const std::string& filename);
std::vector<std::vector<uint64_t>> load_vector_of_vectors_from_file(const std::string& filename);
void save_order(uint32_t W, uint32_t K, double err, const std::vector<uint64_t>& order, bool swapped);
bool does_order_exist(uint32_t W, uint32_t K, double err, bool swapped);
bool does_order_exists_specific(uint32_t W, uint32_t K, double err, bool swapped, const std::string& name);
std::vector<uint64_t> load_order(uint32_t W, uint32_t K, double err, bool swapped);
void save_order_specific(uint32_t W, uint32_t K, double err, const std::vector<uint64_t>& order, bool swapped, const std::string& name);
std::vector<uint64_t> load_order_specific(uint32_t W, uint32_t K, double err, bool swapped, std::string name);
void ensure_directories_exist();

double get_error_with_noise(double error, double noise);

void save_2d_to_csv(const std::map<uint32_t, std::map<uint32_t, double>>& data, const std::string& file_name, const std::string& row_name, const std::string& col_name, Config& config);

void save_1d_to_csv(const std::map<uint32_t, double>& data, const std::string file_name, Config& config);

void print_to_both(Config &config, const std::string& message);

std::tuple<std::vector<uint64_t>, std::vector<uint64_t>> load_odd_even_pair_from_file(const std::string& path, uint32_t w, uint32_t k);
#endif // TOOLS_H