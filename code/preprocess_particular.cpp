#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <boost/multiprecision/cpp_int.hpp> // Include Boost.Multiprecision
#include "config.h"
#include "tools.h"
#include "preprocess_particular.h"


using boost::multiprecision::cpp_int;


void process_sequence(Config& config, uint32_t w, uint32_t k) {
    namespace fs = std::filesystem;

    uint32_t window_size = w + k;

    // Construct the output directory
    std::string output_dir = "temp/preprocessed_contexts/" + config.name;
    if (!fs::exists(output_dir)) {
        fs::create_directories(output_dir);
    }

    std::string output_file = output_dir + "/contexts_" + std::to_string(window_size) + "_odd_even.csv";

    std::ifstream fasta_file(config.path);
    if (!fasta_file.is_open()) {
        print_to_both(config, "Error opening fasta file: " + config.path + "\n");
        return;
    }

    std::string line, sequence = "";
    bool first_sequence_found = false;

    while (std::getline(fasta_file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (first_sequence_found) break; // Stop after the first sequence
            first_sequence_found = true;
            continue; // Skip the header line
        }
        sequence += line;
    }

    fasta_file.close();

    for (auto& c : sequence) c = toupper(c); // Convert to uppercase

    std::unordered_map<char, std::string> binary_mapping = {
        {'A', "00"}, {'C', "01"}, {'G', "10"}, {'T', "11"}
    };

    std::ofstream result_file(output_file);
    if (!result_file.is_open()) {
        print_to_both(config, "Error creating results file: " + output_file);
        return;
    }

    result_file << "OddNumber,EvenNumber\n";

    for (size_t i = 0; i <= sequence.size() - window_size; ++i) {
        std::string window = sequence.substr(i, window_size);
        std::string binary_string = "";
        bool invalid_nucleotide = false;

        for (char nucleotide : window) {
            if (binary_mapping.find(nucleotide) != binary_mapping.end()) {
                binary_string += binary_mapping[nucleotide];
            }
            else {
                invalid_nucleotide = true;
                break;
            }
        }
        if (invalid_nucleotide) continue;

        std::string odd_positions = "", even_positions = "";
        for (size_t j = 0; j < binary_string.size(); ++j) {
            if (j % 2 == 0) {
                odd_positions += binary_string[j];
            }
            else {
                even_positions += binary_string[j];
            }
        }

        // Manually convert binary strings to cpp_int
        cpp_int odd_number = 0;
        cpp_int even_number = 0;

        for (char bit : odd_positions) {
            odd_number = (odd_number << 1) | (bit - '0'); // Shift and add bit
        }

        for (char bit : even_positions) {
            even_number = (even_number << 1) | (bit - '0'); // Shift and add bit
        }

        result_file << odd_number << "," << even_number << "\n";
    }

    result_file.close();
    //print_to_both(config, "Processing completed. Results saved to " + output_file);
    print_to_both(config, "Sequence processing completed\n");

}


bool is_sequence_processed(Config& config, uint32_t w, uint32_t k) {
    namespace fs = std::filesystem;

    // Calculate window size
    uint32_t window_size = w + k;

    // Construct output directory
    std::string output_dir = "temp/preprocessed_contexts/" + config.name;
    if (!fs::exists(output_dir)) {
        return false;
    }

    // Construct output file path
    std::string output_file = output_dir + "/contexts_" + std::to_string(window_size) + "_odd_even.csv";

    // Check if output file exists
    if (fs::exists(output_file)) {
        return true;
    }

    return false;
}

// TODO: this needs ot be void...
bool ensure_sequence_is_processed(Config& config, uint32_t w, uint32_t k) {
    uint32_t window_size = w + k;
	if (is_sequence_processed(config, w, k)) {
		print_to_both(config, "Sequence has already been processed for w+k = " + std::to_string(window_size) + "\n");
		return true;
	}
	else {
		print_to_both(config, "Processing sequence for w+k = " + std::to_string(window_size) + "\n");
		process_sequence(config, w, k);
		return true;
	}
}


