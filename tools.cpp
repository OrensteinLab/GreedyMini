#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include "tools.h"
#include <sstream>  // Include this for istringstream
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <fstream>
#include <filesystem>
#include <set>
#include "tests.h"
#include "config.h"
#include <random>

void save_2d_to_csv(const std::map<uint32_t, std::map<uint32_t, double>>& data, const std::string& file_name, const std::string& row_name, const std::string& col_name, Config& config) {
    // Create folder for csvs
    std::filesystem::create_directories("output/csv");

    // open the file with trunc mode to ensure it's empty
    std::string dir_path = "output/v_" + config.version_id + "/csv";
    std::filesystem::create_directories(dir_path);

    // Construct the full file path
    std::string file_path = dir_path + "/" + file_name + ".csv";

    // Open the file with trunc mode to ensure it's empty
    std::ofstream file(file_path, std::ios::trunc);


    // Check if the file is open
    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file.");
    }

    // Collect all unique column keys
    std::set<uint32_t> column_keys;
    for (const auto& outer_pair : data) {
        for (const auto& inner_pair : outer_pair.second) {
            column_keys.insert(inner_pair.first);
        }
    }

    // Write the first row (column headers)
    file << row_name << "/" << col_name;
    for (const auto& col_key : column_keys) {
        file << "," << col_key;
    }
    file << std::endl;

    // For each row in data
    for (const auto& outer_pair : data) {
        uint32_t row_key = outer_pair.first;
        const auto& inner_map = outer_pair.second;

        // Write the row key
        file << row_key;

        // For each column key in the set of column keys
        for (const auto& col_key : column_keys) {
            // Check if the inner map has the column key
            auto it = inner_map.find(col_key);
            if (it != inner_map.end()) {
                // Write the value
                file << "," << it->second;
            }
            else {
                // Write placeholder (empty cell)
                file << ",";
            }
        }
        file << std::endl;
    }

    // Close the file
    file.close();
}

void save_1d_to_csv(const std::map<uint32_t, double>& data, const std::string file_name, Config& config) {
    // Create folder for csvs
    std::filesystem::create_directories("output/csv");

    // open the file with trunc mode to ensure it's empty
    std::string dir_path = "output/v_" + config.version_id + "/csv";
    std::filesystem::create_directories(dir_path);

    // Construct the full file path
    std::string file_path = dir_path + "/" + file_name + ".csv";

    // Open the file with trunc mode to ensure it's empty
    std::ofstream file(file_path, std::ios::trunc);


    // Check if the file is open
    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file.");
    }

    // write the first row as comma serperated values of the keys
    for (auto const& x : data) {
		file << x.first << ",";
	}
	file << std::endl;

	// write the second row as comma serperated values of the values
	for (auto const& x : data) {
		file << x.second << ",";
	}
	file << std::endl;

	// close the file
	file.close();
}

void print_to_both(Config& config, const std::string& message) {
    // Print to console
    std::cout << message << std::flush;

    // Check if the log file is open and print to it
    if (config.log_file.is_open()) {
        config.log_file << message << std::flush;
    }
    else {
        std::cerr << "Log file is not open or not initialized!" << std::endl;
    }

    return;
}

void save_vector_to_file(const std::vector<uint64_t>& vec, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        size_t size = vec.size();
        file.write(reinterpret_cast<const char*>(&size), sizeof(size));
        file.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(uint64_t));
        file.close();
    }
    else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
    }
}

std::vector<uint64_t> load_vector_from_file(const std::string& filename) {
    std::vector<uint64_t> vec;
    std::ifstream file(filename, std::ios::binary);
    if (file.is_open()) {
        size_t size;
        file.read(reinterpret_cast<char*>(&size), sizeof(size));
        vec.resize(size);
        file.read(reinterpret_cast<char*>(vec.data()), size * sizeof(uint64_t));
        file.close();
    }
    else {
        std::cerr << "Unable to open file for reading: " << filename << std::endl;
    }
    return vec;
}


std::tuple<std::vector<uint64_t>, std::vector<uint64_t>> load_odd_even_pair_from_file(const std::string& path, uint32_t w, uint32_t k) {
        // Copy the path to a temporary string
        std::string path_temp = path;
        path_temp += "/contexts_";
        path_temp += std::to_string(w + k);
        path_temp += "_odd_even.csv";

        //std::cout << "Loading odd-even pair from file: " << path_temp << std::endl;

        // open the file for reading
        std::ifstream file(path_temp);
        std::string line;


        std::vector<uint64_t> oddNumbers;
        std::vector<uint64_t> evenNumbers;

        // Check if file is open
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << path_temp << std::endl;
            return { oddNumbers, evenNumbers };
        }

        // Read and skip the header line
        std::getline(file, line);

        // Read each line of the file
        while (std::getline(file, line)) {
            // Split the line by tab ('\t')
            std::vector<std::string> tokens = split(line, ',');
            if (tokens.size() == 2) {  // Ensure there are two columns
                uint64_t oddNumber = std::stoull(tokens[0]);
                uint64_t evenNumber = std::stoull(tokens[1]);
                oddNumbers.push_back(oddNumber);
                evenNumbers.push_back(evenNumber);
            }
        }

        // Close the file
        file.close();

        return { oddNumbers, evenNumbers };
}

void save_vector_of_vectors_to_file(const std::vector<std::vector<uint64_t>>& vec_of_vecs, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        // Save the outer vector size
        size_t outer_size = vec_of_vecs.size();
        file.write(reinterpret_cast<const char*>(&outer_size), sizeof(outer_size));

        // Save each inner vector
        for (const auto& inner_vec : vec_of_vecs) {
            // Save the size of the inner vector
            size_t inner_size = inner_vec.size();
            file.write(reinterpret_cast<const char*>(&inner_size), sizeof(inner_size));

            // Save the elements of the inner vector
            file.write(reinterpret_cast<const char*>(inner_vec.data()), inner_size * sizeof(uint64_t));
        }

        file.close();
    }
    else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
    }
}

std::vector<std::vector<uint64_t>> load_vector_of_vectors_from_file(const std::string& filename) {
    std::vector<std::vector<uint64_t>> vec_of_vecs;
    std::ifstream file(filename, std::ios::binary);
    if (file.is_open()) {
        // Load the outer vector size
        size_t outer_size;
        file.read(reinterpret_cast<char*>(&outer_size), sizeof(outer_size));
        vec_of_vecs.resize(outer_size);

        // Load each inner vector
        for (auto& inner_vec : vec_of_vecs) {
            // Load the size of the inner vector
            size_t inner_size;
            file.read(reinterpret_cast<char*>(&inner_size), sizeof(inner_size));
            inner_vec.resize(inner_size);

            // Load the elements of the inner vector
            file.read(reinterpret_cast<char*>(inner_vec.data()), inner_size * sizeof(uint64_t));
        }

        file.close();
    }
    else {
        std::cerr << "Unable to open file for reading: " << filename << std::endl;
    }
    return vec_of_vecs;
}
//std::vector<std::vector<uint64_t>> load_vector_of_vectors_from_file(const std::string& filename) {
//    std::ifstream file(filename, std::ios::binary | std::ios::ate); // Open in binary mode and move to the end of the file
//    if (!file.is_open()) {
//        std::cerr << "Unable to open file for reading: " << filename << std::endl;
//        return {};
//    }
//
//    // Get the size of the file
//    std::streamsize file_size = file.tellg();
//    file.seekg(0, std::ios::beg);
//
//    // Read the entire file into a buffer
//    std::vector<char> buffer(file_size);
//    if (!file.read(buffer.data(), file_size)) {
//        std::cerr << "Error reading file: " << filename << std::endl;
//        return {};
//    }
//    file.close();
//
//    // Process the buffer to fill vec_of_vecs
//    const char* data_ptr = buffer.data();
//    size_t outer_size;
//    std::memcpy(&outer_size, data_ptr, sizeof(outer_size));
//    data_ptr += sizeof(outer_size);
//
//    std::vector<std::vector<uint64_t>> vec_of_vecs(outer_size);
//
//    for (auto& inner_vec : vec_of_vecs) {
//        size_t inner_size;
//        std::memcpy(&inner_size, data_ptr, sizeof(inner_size));
//        data_ptr += sizeof(inner_size);
//
//        inner_vec.resize(inner_size);
//        std::memcpy(inner_vec.data(), data_ptr, inner_size * sizeof(uint64_t));
//        data_ptr += inner_size * sizeof(uint64_t);
//    }
//
//    return vec_of_vecs;
//}




bool file_exists(const std::string& filename) {
    return std::filesystem::exists(filename);
}

void save_order(uint32_t W, uint32_t K, double min_alpha, double max_alpha, const std::vector<uint64_t>& order, bool swapped) {
	// construct the file path based on whether swapped is true or false
	std::string filePath = "output/minimizers/" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(min_alpha) + "_" + std::to_string(max_alpha);
	if (swapped) {
		filePath += "_swapped";
	}
	filePath += ".bin";

	// save the order to file
	save_vector_to_file(order, filePath);
}

bool does_order_exist(uint32_t W, uint32_t K, double min_alpha, double max_alpha, bool swapped) {
	// construct the file path based on whether swapped is true or false
	std::string filePath = "output/minimizers/" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(min_alpha) + "_" + std::to_string(max_alpha);
	if (swapped) {
		filePath += "_swapped";
	}
	filePath += ".bin";

	// check if the file exists
	return file_exists(filePath);
}

bool does_order_exists_specific(uint32_t W, uint32_t K, double min_alpha, double max_alpha, bool swapped, const std::string& name) {
    // make sure subfolder exists
    std::filesystem::create_directories("output/minimizers/" + name);

    // construct the file path based on whether swapped is true or false
    std::string filePath = "output/minimizers/" + name + "/" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(min_alpha) + "_" + std::to_string(max_alpha);
    if (swapped) {
        filePath += "_swapped";
    }
    filePath += ".bin";

    // check if the file exists
    return file_exists(filePath);
}

std::vector<uint64_t> load_order(uint32_t W, uint32_t K, double min_alpha, double max_alpha, bool swapped) {
	// construct the file path based on whether swapped is true or false
	std::string filePath = "output/minimizers/" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(min_alpha) + "_" + std::to_string(max_alpha);
	if (swapped) {
		filePath += "_swapped";
	}
	filePath += ".bin";

	// check if the file exists
	if (!file_exists(filePath)) {
		std::cerr << "File does not exist: " << filePath << std::endl;
        // stop the program
        exit(1);
	}

	// load original order from file
	return load_vector_from_file(filePath);
}

void save_order_specific(uint32_t W, uint32_t K, double min_alpha, double max_alpha, const std::vector<uint64_t>& order, bool swapped, const std::string& name) {
    // make sure subfolder exists
    std::filesystem::create_directories("output/minimizers/" + name);

    // construct the file path based on whether swapped is true or false
	std::string filePath = "output/minimizers/" + name + "/" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(min_alpha) + "_" + std::to_string(max_alpha);
	if (swapped) {
		filePath += "_swapped";
	}
	filePath += ".bin";

	// save the order to file
	save_vector_to_file(order, filePath);

}

std::vector<uint64_t> load_order_specific(uint32_t W, uint32_t K, double min_alpha, double max_alpha , bool swapped, std::string name) {
    // make sure subfolder exists
    std::filesystem::create_directories("output/minimizers/" + name);

    // construct the file path based on whether swapped is true or false
    std::string filePath = "output/minimizers/" + name + "/" + std::to_string(W) + "_" + std::to_string(K) + "_" + std::to_string(min_alpha) + "_" + std::to_string(max_alpha);
    if (swapped) {
        filePath += "_swapped";
    }
    filePath += ".bin";

    // check if the file exists
    if (!file_exists(filePath)) {
        std::cerr << "File does not exist: " << filePath << std::endl;
        // stop the program
        exit(1);
    }

    // load original order from file
    return load_vector_from_file(filePath);
}

void ensure_directories_exist() {
    std::filesystem::create_directories("temp");
    std::filesystem::create_directories("logs");
    std::filesystem::create_directories("output/minimizers");

}


double sample_alpha(Config& config) {
    // sample uniformly for log(1-err)
    double min_rescaled = std::log(1.0 - config.max_alpha);
    double max_rescaled = std::log(1.0 - config.min_alpha);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(min_rescaled, max_rescaled);
    double x = distribution(gen);
    double alpha = 1.0 - std::exp(x);
    return alpha;
}

double get_error_with_noise(double error, double noise) {
    // Compute log(1 - err)
    double log1m_err = std::log(1.0 - error);
    // Define the sampling range for x
    double x_min = log1m_err - noise;
    double x_max = std::min(0.0, log1m_err + noise);
    // Initialize random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(x_min, x_max);
    // Sample x and update err
    double x = distribution(gen);
    error = 1.0 - std::exp(x);

    return error;

}

std::vector<uint64_t> get_reversed_order(const std::vector<uint64_t>& order) {
	std::vector<uint64_t> reversed_order(order.size());
    // the length is the number of digits on binary of the size
    uint32_t length = static_cast<uint32_t>(std::log2(order.size()));
    //std::cout << "k: " << length << std::endl;

	for (uint32_t i = 0; i < order.size(); ++i) {
		// get the reversed binary of the first length bits of i
        uint64_t reversed = 0;
        for (uint32_t j = 0; j < length; ++j) {
			reversed |= ((i >> j) & 1) << (length - 1 - j);
		}
        reversed_order[reversed] = order[i];
	}
	return reversed_order;
}


// Function to split a string by a delimiter
std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<uint64_t> get_explicitly_extended_order(const std::vector<uint64_t>& order) {
    uint64_t original_size = static_cast<uint64_t>(order.size());
    uint64_t new_size = 2 * original_size;
    std::vector<uint64_t> extended_order(new_size);

    std::cout << "(DEBUG) New size: " << new_size << std::endl;

    for (uint64_t i = 0; i < new_size; ++i) {
		uint64_t first_k_bits = i >> 1;
        uint64_t last_bit = i & 1;
        uint64_t new_value = 2 * order[first_k_bits] + last_bit;
        // if above 2^k (can be 2^k+1) then set to 2^k
        new_value = std::min(new_value, new_size);
        extended_order[i] = new_value;
	}
    return extended_order;
}