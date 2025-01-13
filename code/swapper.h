#pragma once
#ifndef SWAPPER_H
#define SWAPPER_H

#include <vector>
#include <cstdint>
#include <utility> 

uint64_t cost(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t u, uint64_t v, uint64_t y);

uint64_t cost_recursive(uint32_t w, uint32_t k, const std::vector<uint64_t>& order, uint64_t u, uint64_t v, uint64_t y);


std::pair<std::vector<uint64_t>, uint64_t> swapper_f(uint32_t w, uint32_t k, std::vector<uint64_t>& order, double max_time_seconds, bool verbose);

std::pair<std::vector<uint64_t>, uint64_t> swapper_f_v2(uint32_t w, uint32_t k, std::vector<uint64_t>& order, double max_time_seconds, bool verbose);

#endif // SWAPPER_H
