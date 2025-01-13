#include <iostream>
#include <chrono>
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
#include <fstream>
#include <sstream>
#include <string>
#include "tests.h"
#include "config.h"
#include <cstdlib>  
#include <ctime>    
#include <cstdlib>  
#include <stdexcept>
#include "Minimizers.h"



//#include "sampling_runtimes.h"




void update_config_particular(Config &config) {
    config.path = "chr_x_1m.fasta";
    config.name = "chr_x_1m";
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
    // USED FOR TESTING ONLY, NOT PART OF THE OFFICIAL IMPLEMENTATION
    //perform_all_sampling_tests();
    //return 0;



    // Initialize variables with default values
    std::string mode;
    std::string path = "None";
    std::string name = "None";
    std::string version_id;
    uint32_t greedy_mini_runs = 4096;
    uint32_t max_swapper_time_minutes = std::numeric_limits<uint32_t>::max();
    std::string output_folder = "final";    
    double min_alpha = 0.939088;
    double max_alpha = 0.999590;
    uint32_t w;
    uint32_t k;
    uint32_t n_cores = static_cast<uint32_t>(std::thread::hardware_concurrency() / 2);



    bool success = parse_arguments(argc, argv, mode, path, name, greedy_mini_runs, max_swapper_time_minutes, output_folder, min_alpha, max_alpha, version_id, w,k, n_cores);
    if (!success) return -1;



    // Create Config object with all parameters
    Config config(
        mode,
        path,
        name,
        greedy_mini_runs,
        max_swapper_time_minutes,
        min_alpha,
        max_alpha,
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
        clean_up_particular_temp_files(config, true);
        print_to_both(config, "Finished tests\n");
        config.log_file.close();

    } else if (config.mode == "tests") {
        update_config_expected(config);
        print_to_both(config, "Starting tests\n");
        expected_density_tests(config, 15, 15);
        print_to_both(config, "Finished expected density tests\n");
        config.log_file.close();

        update_config_particular(config);
        print_to_both(config, "Starting tests\n");
        particular_dna_tests(config);
        clean_up_particular_temp_files(config, true);
        print_to_both(config, "Finished particular density tests\n");
        config.log_file.close();
	} else if (config.mode == "expected") {
        short_compute_and_store_densities(w, k, config);

    }
    else if (config.mode == "particular") {
        short_calculate_particular_density(w, k, config);
        clean_up_particular_temp_files(config, true);
    } else if (config.mode == "swapper") {
        print_to_both(config, "Starting swaps\n");
        swaps_only(config);
        print_to_both(config, "Finished swaps\n");
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
    double& min_alpha,
    double& max_alpha,
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
        else if (arg == "min_alpha") {
            min_alpha = std::stod(argv[++i]);
        }
        else if (arg == "max_alpha") {
            max_alpha = std::stod(argv[++i]);
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

    if (mode != "expected" && mode != "particular" && mode != "tests_e" && mode != "tests_p" && mode != "tests" && mode != "swapper") {
        std::cerr << "Error: '--mode' must be 'expected', 'particular', 'swapper', 'tests_e', 'tests_p', or 'tests'.\n";
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

    if (mode == "swapper" && (path == "None" || !w_provided || !k_provided)) {
		std::cerr << "Error: '--path', '--w', and '--k' are required when mode is 'swapper'.\n";
		return false;
	}

    // check min_alpha is between 0 and 1
    if (min_alpha <= 0 || min_alpha >= 1) {
		std::cerr << "Error: '--min_alpha' must in the range of (0, 1).\n";
		return false;
	}

    // check max_alpha is between 0 and 1
    if (max_alpha <= 0 || max_alpha >= 1) {
        std::cerr << "Error: '--max_alpha' must in the range of (0, 1).\n";
        return false;
    }

    // check min_alpha is less or equal to max_alpha
    if (min_alpha > max_alpha) {
		std::cerr << "Error: '--min_alpha' must not be greater than '--max_alpha'.\n";
		return false;
	}

    // check that w + k < 64
    if (w + k >= 64) {
		std::cerr << "Error: '--w' + '--k' must be strictly less than 64.\n";
		return false;
	}

    // Use output_folder as version_id
    version_id = output_folder;

    return true; // Indicate success
}




