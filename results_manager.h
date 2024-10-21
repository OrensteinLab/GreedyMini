#ifndef RESULTS_MANAGER_H
#define RESULTS_MANAGER_H

#include <string>
#include <unordered_map>
#include <vector>
#include <mutex>

class ResultsManager {
public:
    ResultsManager(const std::string& directory, const std::string& filename);
    void update_cell(uint32_t w, uint32_t k, double value);
    void save_to_file();

private:
    void load_file();
    void create_directories();

    std::string directory;
    std::string filepath;
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, double>> data; // Stores the CSV data
    std::vector<uint32_t> w_values;
    std::vector<uint32_t> k_values;

    void ensure_column_exists(uint32_t w);
    void ensure_row_exists(uint32_t k);

    std::mutex mtx; // Mutex to protect the critical sections
};

#endif // RESULTS_MANAGER_H
