#include "results_manager.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm>

ResultsManager::ResultsManager(const std::string& directory, const std::string& filename)
    : directory(directory), filepath(directory + "/" + filename) {
    create_directories();
    load_file();
}

void ResultsManager::create_directories() {
    std::lock_guard<std::mutex> lock(mtx); // Lock the mutex to protect this section
    if (!std::filesystem::exists(directory)) {
        if (!std::filesystem::create_directories(directory)) {
            std::cerr << "Error creating directories: " << directory << std::endl;
        }
    }
}

void ResultsManager::load_file() {
    std::lock_guard<std::mutex> lock(mtx); // Lock the mutex to protect this section
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cout << "File does not exist, starting with an empty file: " << filepath << std::endl;
        return; // If the file doesn't exist, we'll start with an empty map
    }

    std::string line;
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::getline(ss, value, ','); // Skip the top-left corner cell
        while (std::getline(ss, value, ',')) {
            w_values.push_back(std::stoi(value));
        }
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        std::getline(ss, value, ',');
        uint32_t k = std::stoi(value); // Convert the first value in the line to uint32_t (k)
        k_values.push_back(k); // Add k to k_values

        for (size_t i = 0; i < w_values.size(); ++i) {
            if (std::getline(ss, value, ',')) {
                if (!value.empty()) {
                    data[k][w_values[i]] = std::stod(value); // Convert value to double and store in data[k][w_values[i]]
                }
                else {
                    data[k][w_values[i]] = 0.0; // Handle empty value, assign a default value (e.g., 0.0)
                }
            }
        }
    }


    //std::cout << "Loaded file successfully: " << filepath << std::endl;
}

void ResultsManager::ensure_column_exists(uint32_t w) {
    std::lock_guard<std::mutex> lock(mtx); // Lock the mutex to protect this section
    if (std::find(w_values.begin(), w_values.end(), w) == w_values.end()) {
        w_values.push_back(w);
        std::sort(w_values.begin(), w_values.end());
    }
}

void ResultsManager::ensure_row_exists(uint32_t k) {
    std::lock_guard<std::mutex> lock(mtx); // Lock the mutex to protect this section
    if (std::find(k_values.begin(), k_values.end(), k) == std::end(k_values)) {
        k_values.push_back(k);
        std::sort(k_values.begin(), k_values.end());
    }
}

void ResultsManager::update_cell(uint32_t w, uint32_t k, double value) {
    ensure_column_exists(w);
    ensure_row_exists(k);
    data[k][w] = value;
}

void ResultsManager::save_to_file() {
    std::lock_guard<std::mutex> lock(mtx); // Lock the mutex to protect this section
    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Unable to open file for writing: " << filepath << std::endl;
        return;
    }

    // Write the header
    file << ",";
    for (auto w : w_values) {
        file << w << ",";
    }
    file << "\n";

    // Write each row
    for (auto k : k_values) {
        file << k << ",";
        for (auto w : w_values) {
            file << (data[k].count(w) ? std::to_string(data[k][w]) : "") << ",";
        }
        file << "\n";
    }

    //std::cout << "File saved successfully: " << filepath << std::endl;
}
