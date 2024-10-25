#include "window_types_manager.h"

// Method to add a binary vector to the dictionary or increment its count
void WindowTypesManager::add_vector(const uint64_t binary_vector) {
    if (vector_counts.find(binary_vector) != vector_counts.end()) {
        vector_counts[binary_vector] += 1; // Increment count if it already exists
    }
    else {
        vector_counts[binary_vector] = 1; // Create new entry with count 1
    }
}

// Method to get the count of a specific binary vector
uint64_t WindowTypesManager::get_count(const uint64_t binary_vector) const {
    auto it = vector_counts.find(binary_vector);
    if (it != vector_counts.end()) {
        return it->second;
    }
    return 0; // Return 0 if the binary vector is not found
}

// will return the sum of the occurance of each special case times the probability of each special case to be a gamechanger
uint64_t WindowTypesManager::get_special_cases_count(uint32_t w, uint32_t k) const {
	uint64_t count = 0;
	for (auto const& x : vector_counts) {
		if (x.first == 0 || x.first == (1ULL << (w + k)) - 1) {
			count += x.second;
		}
	}
	return count;
}

// Method to reset the counts
void WindowTypesManager::reset() {
    vector_counts.clear();
}
