#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <boost/multiprecision/cpp_int.hpp> // Include Boost.Multiprecision
#include "config.h"
#include "tools.h"
#include "preprocess_particular.h"


bool ensure_sequence_is_processed(Config& config, uint32_t w, uint32_t k);
