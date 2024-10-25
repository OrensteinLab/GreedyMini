cmake_minimum_required(VERSION 3.15)
project(GreedyMini)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type to Release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Find required Boost libraries
find_package(Boost REQUIRED COMPONENTS
  system
  filesystem
  # Add other components as needed
)

# Include Boost headers
include_directories(${Boost_INCLUDE_DIRS})

# Collect all .cpp files
file(GLOB SOURCES "*.cpp")

# Create the executable
add_executable(GreedyMini ${SOURCES})

# Link against Boost libraries
target_link_libraries(GreedyMini ${Boost_LIBRARIES})
