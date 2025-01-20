# GreedyMini

A tool to create low-density minimizer orders using a greedy approach. The repository for the paper *Generating low-density minimizers*. We implemented two variants of GreedyMini: **GreedyMini+**, for generating low expected density minimizer orders, and **GreedyMiniParticular+**, for generating low particular density minimizer orders.

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation and Compilation](#installation-and-compilation)
  - [Precompiled Binaries](#precompiled-binaries)
  - [Compiling the Project Manually](#compiling-the-project-manually)
- [Usage](#usage)
  - [Generating Minimizers For Expected Density](#generating-minimizers-for-expected-density)
  - [Generating Minimizers For Particular Density](#generating-minimizers-for-particular-density)
  - [Additional Parameters](#additional-parameters)
  -  [Running Paper Tests](#running-paper-tests)
- [Accessing the Minimizers](#accessing-the-minimizers)
  - [Locating the Minimizers](#locating-the-minimizers)
  - [Loading the Minimizers to Memory](#loading-the-minimizers-to-memory)
- [Contact](#contact)

## Introduction

This document provides step-by-step instructions to use the **GreedyMini** variants.

## Prerequisites

Before you begin, ensure you have the following:

- **64-bit system**: Required to run the binaries or compile GreedyMini.
- **C++ Compiler**: Only required if you are compiling GreedyMini yourself (supports C++20, e.g., `g++` version 10 or higher).

## Installation and Compilation

### Precompiled Binaries

Precompiled binaries are available for **Ubuntu** and **macOS**. You can download them from the [GitHub release page](https://github.com/OrensteinLab/GreedyMini/releases).


#### Build Compatibility Note
To ensure that the compiled Linux binaries do not depend on newer GLIBC or libstdc++ symbols (e.g., `GLIBC_2.38`, `GLIBCXX_3.4.31`), we build on Ubuntu 20.04. This ensures that our released executables can run on older Linux distributions without encountering linker errors about missing GLIBC/GLIBCXX versions.

If your system is significantly older (e.g., it has GLIBC < 2.31), you may still need to build from source. Alternatively, you can run in a newer environment such as Docker or a more up-to-date Linux distribution.


### Compiling the Project Manually

If you prefer to compile the project from source, follow the steps below.

#### Step 1: Check for a Compiler

Verify that you have a suitable C++ compiler installed:

```bash
g++ --version
```

You should see output similar to:

```
g++ (GCC) 10.2.0
```

If you don't have `g++` or it's outdated, you may need to use `clang++` or install a newer version locally.

#### Step 2: Install Boost Locally

Download and extract the Boost headers:

```bash
cd ~
wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz
tar -xzf boost_1_82_0.tar.gz
```

This will create a `boost_1_82_0` directory in your home directory.

#### Step 3: Compile the Project

Navigate to the `GreedyMini` directory and compile the project:

```bash
cd ~/GreedyMini
g++ -std=c++20 -O3 -march=native -I ~/boost_1_82_0 *.cpp -o GreedyMini
```

## Usage

After compiling or downloading the binary, you can run **GreedyMini** using different modes and parameters.




### Generating Minimizers For Expected Density

To generate a minimizer with low expected density, run:

```bash
./GreedyMini -mode expected -w {w} -k {k}
```

Replace `{w}` and `{k}` with your desired values:

- `{w}`: The window size.
- `{k}`: The k-mer size.

> **Note:** Ensure that **`w + k < 64`** (due to 64 bit reliance).

An example run would be:
```
./GreedyMini -mode expected -w 5 -k 4
```

### Generating Minimizers For Particular Density

#### Sequence format

`GreedyMiniParticular+` expects a `.fasta` file containing exactly one sequence.

#### Running

Use the following command to generate minimizers with particular density:

```bash
./GreedyMini -mode particular -w {w} -k {k} -path {path} -name {name}
```

Where:

- `{w}`: The window size.
- `{k}`: The k-mer size.
- `{path}`: a path to the `fasta` file containing the sequence.
- `{name}`: The name for the generated orders.

> **Note:** Ensure that **`w + k < 64`** (due to 64 bit reliance).

In our paper we used `chr_x_1m.fasta` as the path.

An example run would be:
```
./GreedyMini -mode particular -w 5 -k 4 -path chr_x_1m.fasta -name 1M
```
### Additional Parameters

Customize the behavior of GreedyMini with the following options:

- `-greedy_mini_runs`: Number of runs of GreedyMini (default: `4096`)
- `-n_cores`: Number of CPU cores to use (default: half the total available threads due to hyper-threading)
- `-min_alpha`: Minimum alpha value (default: `0.939088`)
- `-max_alpha`: Maximum alpha value (default: `0.999590`)
- `-max_swapper_time_minutes`: Maximum swapper time in minutes

### Running Paper Tests

We ran our test on the first 1M nucleotides of chromosome X from [Genome assembly T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/). To do that we used the python notebook `shorten_fasta.ipynb` which is located in `various scripts/preprocessing chr x/`. We then put the resulting `.fasta` file in the same directory as the `GreedyMini` executable.

Execute the following command to run most of the tests from the paper:

```bash
./GreedyMini -mode tests
```


## Accessing the Minimizers

### Locating the Minimizers

Generated minimizers will appear inside the `output/minimizers` folder. For particular density minimizers, they will appear in a subfolder with the selected name.

Example filename:

```
{w}_{k}_{min_alpha}_{max_alpha}_swapped.bin
```

See [Additional Parameters](#additional-parameters) for the default values of `min_alpha` and `max_alpha`.

### Loading the Minimizers to Memory

We provide the Python notebook `load_order.ipynb` alongside the best minimizer orders from the paper, both located in the folder `minimizer loading example`. The notebook showcases how to load a minimizer to memory and print the order of each k-mer. For C++, we recommend looking at the functions `load_order()` and `load_vector_from_file()` located in `code/tools.cpp`.


## Contact
In case of issues with the tool, you may contact us at tziony.i@gmail.com.
