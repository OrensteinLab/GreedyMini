#pragma once

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
    uint32_t& n_cores);