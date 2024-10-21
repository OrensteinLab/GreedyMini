#include <iostream>
#include <chrono>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <limits>
#include "functions.h" 
#include "runs.h"  
#include "tools.h"
#include "density_and_gc.h"
#include <algorithm> 
#include <cmath>
#include <random>     
#include <map>
#include <iomanip> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "tests.h"
#include "config.h"
#include <iostream>
#include <cstdlib>  
#include <ctime>    
#include <cstdlib>  // For std::stoul, std::stod
#include <stdexcept>
#include "Minimizers.h"




void update_config_particular(Config &config) {
    config.path = "sequences_1M";
    config.name = "chromosomeX_1m";
    config.version_id = "final-1m";
    config.max_mins_per_step = 60;

    // Re-initialize resources if necessary
    // For example, reopen the log file with the new version_id
    if (config.log_file.is_open()) {
        config.log_file.close();
    }
    std::string output_dir = "output/v_" + config.version_id;
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
    }
    std::string log_file_path = output_dir + "/log_file.txt";
    config.log_file.open(log_file_path);
    if (!config.log_file.is_open()) {
        std::cerr << "Error: Could not open log file at " << log_file_path << std::endl;
    }
}

void update_config_expected(Config &config) {
    config.path = "None";
    config.name = "None";
    config.version_id = "final";
    config.max_mins_per_step = 60;


    // Re-initialize resources if necessary
    if (config.log_file.is_open()) {
        config.log_file.close();
    }
    std::string output_dir = "output/v_" + config.version_id;
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
    }
    std::string log_file_path = output_dir + "/log_file.txt";
    config.log_file.open(log_file_path);
    if (!config.log_file.is_open()) {
        std::cerr << "Error: Could not open log file at " << log_file_path << std::endl;
    }
}



int main(int  argc, char* argv[])
{ 
    // Initialize variables with default values
    std::string mode;
    std::string path = "None";
    std::string name = "None";
    std::string version_id;
    uint32_t greedy_mini_runs = 4096;
    uint32_t max_swapper_time_minutes = std::numeric_limits<uint32_t>::max();
    std::string output_folder = "final";    
    double error = 0.995;
    double noise_for_error = 2.5;
    uint32_t w;
    uint32_t k;
    uint32_t n_cores = static_cast<uint32_t>(std::thread::hardware_concurrency() / 2);



    bool success = parse_arguments(argc, argv, mode, path, name, greedy_mini_runs, max_swapper_time_minutes, output_folder, error, noise_for_error, version_id, w,k, n_cores);
    if (!success) return -1;



    // Create Config object with all parameters
    Config config(
        mode,
        path,
        name,
        greedy_mini_runs,
        max_swapper_time_minutes,
        error,
        noise_for_error,
        version_id,
        w,
        k,
        max_swapper_time_minutes,
        n_cores
    );

    ensure_directories_exist();


    if (config.mode == "tests_e") {
        update_config_expected(config);
        print_to_both(config, "Starting tests\n");
        expected_density_tests(config, 15, 15);
        print_to_both(config, "Finished tests\n");
        config.log_file.close();
    }
    else if (config.mode == "tests_p") {
        update_config_particular(config);
        print_to_both(config, "Starting tests\n");
        particular_dna_tests(config);
        print_to_both(config, "Finished tests\n");
        config.log_file.close();
    } else if (config.mode == "tests") {
        update_config_expected(config);
        print_to_both(config, "Starting tests\n");
        expected_density_tests(config, 15, 15);
        print_to_both(config, "Finished tests\n");
        config.log_file.close();

        update_config_particular(config);
        print_to_both(config, "Starting tests\n");
        particular_dna_tests(config);
        print_to_both(config, "Finished tests\n");
        config.log_file.close();
	} else if (config.mode == "expected") {
        short_compute_and_store_densities(w, k, config);

	} else if (config.mode == "particular") {
        short_calculate_particular_density(w, k, config);
	} else {
		std::cerr << "Error: Unknown mode '" << config.mode << "'.\n";
		return 1;
	}

    return 0;
}

bool parse_arguments(int argc, char* argv[],
    std::string& mode,
    std::string& path,
    std::string& name,
    uint32_t& greedy_mini_runs,
    uint32_t& max_swapper_time_minutes,
    std::string& output_folder,
    double& error,
    double& noise_for_error,
    std::string& version_id,
    uint32_t& w,
    uint32_t& k,
    uint32_t& n_cores) {
    // Initialize flags to check if w and k were provided
    bool w_provided = false;
    bool k_provided = false;

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Check if arg starts with '-' or '--'
        if (arg.rfind("--", 0) == 0) {
            // Long option, remove '--'
            arg = arg.substr(2);
        }
        else if (arg.rfind("-", 0) == 0) {
            // Short option, remove '-'
            arg = arg.substr(1);
        }
        else {
            std::cerr << "Error: Unknown parameter '" << arg << "'.\n";
            return false;
        }

        // Now arg is the option name without the '-' or '--'
        // Check if the next argument exists
        if (i + 1 >= argc) {
            std::cerr << "Error: '--" << arg << "' requires a value.\n";
            return false;
        }

        // Parse the option
        if (arg == "mode") {
            mode = argv[++i];
        }
        else if (arg == "path") {
            path = argv[++i];
        }
        else if (arg == "name") {
            name = argv[++i];
        }
        else if (arg == "greedy_mini_runs") {
            greedy_mini_runs = static_cast<uint32_t>(std::stoul(argv[++i]));
        }
        else if (arg == "max_swapper_time_minutes") {
            std::string max_swapper_time_minutes_str = argv[++i];
            max_swapper_time_minutes = static_cast<uint32_t>(std::stoul(max_swapper_time_minutes_str));

        }
        else if (arg == "output_folder") {
            output_folder = argv[++i];
        }
        else if (arg == "n_cores") {
			n_cores = static_cast<uint32_t>(std::stoul(argv[++i]));
		}
        else if (arg == "error") {
            error = std::stod(argv[++i]);
        }
        else if (arg == "noise_for_error") {
            noise_for_error = std::stod(argv[++i]);
        }
        else if (arg == "w") {
            w = static_cast<uint32_t>(std::stoul(argv[++i]));
            w_provided = true;
        }
        else if (arg == "k") {
            k = static_cast<uint32_t>(std::stoul(argv[++i]));
            k_provided = true;
        }
        else {
            std::cerr << "Error: Unknown parameter '--" << arg << "'.\n";
            return false;
        }
    }

    // Validate required parameters
    if (mode.empty()) {
        std::cerr << "Error: '--mode' parameter is required.\n";
        return false;
    }

    if (mode != "expected" && mode != "particular" && mode != "tests_e" && mode != "tests_p" && mode != "tests") {
        std::cerr << "Error: '--mode' must be 'expected', 'particular', 'tests_e', 'tests_p', or 'tests'.\n";
        return false;
    }

    if ((mode == "expected" || mode == "particular") && (!w_provided || !k_provided)) {
        std::cerr << "Error: '--w' and '--k' are required when mode is 'expected' or 'particular'.\n";
        return false;
    }

    if (mode == "particular" && (path == "None" || name == "None")) {
        std::cerr << "Error: '--path' and '--name' are required when mode is 'particular'.\n";
        return false;
    }

    // check error is between 0 and 1
    if (error < 0 || error > 1) {
		std::cerr << "Error: '--error' must be between 0 and 1.\n";
		return false;
	}

    // check noise_for_error is positive
    if (noise_for_error < 0) {
        std::cerr << "Error: '--noise_for_error' must be positive.\n";
        return false;
    }


    // Use output_folder as version_id
    version_id = output_folder;

    return true; // Indicate success
}




