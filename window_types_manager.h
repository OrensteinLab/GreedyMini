#ifndef WINDOW_TYPES_MANAGER_H
#define WINDOW_TYPES_MANAGER_H

#include <unordered_map>
#include <cstdint>

class WindowTypesManager {
public:
    // Method to add a binary vector to the dictionary or increment its count
    void add_vector(const uint64_t binary_vector);

    // Method to get the count of a specific binary vector
    uint64_t get_count(const uint64_t binary_vector) const;

    // Method to reset the counts
    void reset();

    uint64_t get_special_cases_count(uint32_t w, uint32_t k) const;

//private:
    // Dictionary to store counts of binary vectors of size (w+k)
    std::unordered_map<uint64_t, uint64_t> vector_counts;

};

#endif // WINDOW_TYPES_MANAGER_H
