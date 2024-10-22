# GreedyMini


## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation and Compilation](#installation-and-compilation)
  - [Step 1: Check for a Compiler](#step-1-check-for-a-compiler)
  - [Step 2: Install Boost Locally](#step-2-install-boost-locally)
  - [Step 3: Compile the Project](#step-3-compile-the-project)
- [Usage](#usage)
  - [Running Tests](#running-tests)
  - [Generating Minimizers](#generating-minimizers)
  - [Generating Minimizers For Particular Density](#generating-minimizers-for-particular-density)
  - [Additional Parameters](#additional-parameters)
- [Accesing the Minimizers](#accessing-the-minimizers)
  -[Locating the Minimizers](#locating-the-minimizers)
  -[Loading the Minimizers to Memory](#loading-the-minimizers-to-memory)


## Introduction

This document provides step-by-step instructions to compile and run the GreedyMini project. The project relies on the Boost Multiprecision library for handling large integers and precise arithmetic operations.

## Prerequisites

Before you begin, ensure you have the following:

- **C++ Compiler**: A compiler that supports C++20, such as `g++` version 10 or higher.

## Installation and Compilation

### Step 1: Check for a Compiler

First, verify that you have a suitable C++ compiler installed:

```bash
g++ --version
```

You should see output similar to:

```
g++ (GCC) 10.2.0
```

If you don't have `g++` or it's outdated, you may need to use an alternative compiler like `clang++` or install a newer version locally.

### Step 2: Install Boost Locally

We'll download the Boost headers to your home directory.

#### 2.1 Download Boost

Navigate to your home directory and download the Boost library:

```bash
cd ~
wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz
```

#### 2.2 Extract Boost

Extract the downloaded archive:

```bash
tar -xzf boost_1_82_0.tar.gz
```

This will create a directory `boost_1_82_0` in your home directory containing all the Boost headers.

### Step 3: Compile the Project

Navigate to the `GreedyMini` directory (where all your project files are located) and compile the project:

```bash
cd ~/GreedyMini
g++ -std=c++20 -O3 -march=native -I ~/boost_1_82_0 *.cpp -o GreedyMini
```

## Usage

After successful compilation, you can run the GreedyMini executable with various modes and parameters.

### Running Tests

Before running all the tests, please refer to [contexts setup](#setting-up-the-contexts) in order to perform tests on particular density.

To run all tests, execute:

```
./GreedyMini -mode tests
```

This will run the suite of tests which we used in our paper to generate expected and particular densities using GreedyMini+ and GreedyMiniParticular+.

### Generating Minimizers For Expected Density

To generate a minimizer with low expected density, use the following command:

```
./GreedyMini -mode expected -w {w} -k {k}
```

Replace `{w}` and `{k}` with your desired values:

- `{w}`: The window size.
- `{k}`: The k-mer size.

An example run would be:
```
./GreedyMini -mode expected -w 5 -k 4
```

### Generating Minimizers For Particular Density


#### Setting up the contexts
In order to run GreedyMiniParticular+, first we need to process our sequence to a specfic format where we first encode the entire DNA sequence as bits, then for each w+k long context we store the odd bits and even bits serperately. We included a python note book `preprocess.ipynb` that expects the file `GCA_009914755.4_T2T-CHM13v2.0_genomic.fna` from [Genome assembly T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) to be in the same folder, and from it creates all such contexts using the first 1M nucleotides of chromosome X.
Running the script generates a folder called `sequences_1M` which should be in the same folder as the executable. In order to use different sequences, only minor modifications are required of the python notebook.

#### Running 

To generate minimizers for a particular density, use the following command:

```
./GreedyMini -mode particular -w {w} -k {k} -path {path} -name {name}
```

Where:

- `{w}`: The window size.
- `{k}`: The k-mer size.
- `{path}`: The name of the folder (located in the same directory as the executable) containing the "contexts" for various `w` and `k` combinations.
- `{name}`: The name for the generated orders.

In our paper we used `sequences_1M` as the path.

An example run would be:
```
./GreedyMini -mode particular -w 5 -k 4 -path sequences_1M -name 1M
```

### Additional Parameters

You can customize the behavior of GreedyMini and GreedyMiniParticular using additional optional parameters:

- `-greedy_mini_runs`: **Number of runs of GreedyMini** (default: `4096`)
- `-n_cores`: **Number of cores** for GreedyMini and swapper (defaults to half the number of available cores)
- `-min_alpha`: **Minimum of the range** in which we sample alpha (as described in the paper, defeault: 0.939088)
- `-max_alpha`: **Maximum of the range** in which we sample alpha (default: 0.999590)
- `-max_swapper_time_minutes`: **Maximum time for the swapper in minutes** (defaults to the runtime of GreedyMini, is only relevant for expected density minimizers) 



## Accessing the Minimizers

### Locating the Minimizers
Generated minimizers will appear inside the `output/minimizers` folder. In case of particular density minimizers, they will appear inside a subfolder with the corresponding picked name.

As a rule of thumb, the most useful minimizer will be named 
```
{w}_{k}_{min_alpha}_{max_alpha}_swapped.bin
```
See [Additional Parameters](#additional-parameters) in order to see default values of `min_alpha` and `max_alpha`.

Other minimizers are saved as steps when building the best ones.

Note that running the tests will produce some orders with `min_alpha` > 1, which are used for internal naming.

### Loading the Minimizers to Memory

For convinience, we added a python notebook which showcases loading a minimizer to memory and prints the order of each k-mer in the `minimizer loading example` folder. Additionaly, for C++ we reccomend looking at the following functions: `load_order()` and `load_vector_from_file()` which are both located inside the file `tools.cpp`.

