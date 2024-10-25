#pragma once
#include <iostream>
#include <fstream>


#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>  // C++17 for directory creation

#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>

class Config {
public:
    std::string mode;
    std::string path;
    std::string name;
    uint32_t greedy_mini_runs;
    uint32_t max_mins_per_step;
    double min_alpha;
    double max_alpha;
    std::string version_id;
    std::ofstream log_file;
    uint32_t w;
    uint32_t k;
    uint32_t max_swapper_time_minutes;
    uint32_t n_cores;

    // Updated constructor to accept new parameters
    Config(const std::string& mode_,
        const std::string& path_,
        const std::string& name_,
        uint32_t greedy_mini_runs_,
        uint32_t max_mins_per_step_,
        double min_alpha_,
        double max_alpha_,
        const std::string& version_id_,
        uint32_t w_,
        uint32_t k_,
        uint32_t max_swapper_time_minutes_,
        uint32_t n_cores)
        : mode(mode_),
        path(path_),
        name(name_),
        greedy_mini_runs(greedy_mini_runs_),
        max_mins_per_step(max_mins_per_step_),
        min_alpha(min_alpha_),
        max_alpha(max_alpha_),
        version_id(version_id_),
        w(w_),
        k(k_),
        max_swapper_time_minutes(max_swapper_time_minutes_),
        n_cores(n_cores) {

        // Create the output directory
        std::string output_dir = "output/v_" + version_id;
        if (!std::filesystem::exists(output_dir)) {
            std::filesystem::create_directories(output_dir);
        }

        // Open the log file
        std::string log_file_path = output_dir + "/log_file.txt";
        log_file.open(log_file_path);
        if (!log_file.is_open()) {
            std::cerr << "Error: Could not open log file at " << log_file_path << std::endl;
        }
    }

    // Destructor to close the log file
    ~Config() {
        if (log_file.is_open()) {
            log_file.close();
        }
    }
};
