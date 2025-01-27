# GreedyMini

A toolkit to create low-density DNA minimizer orders using a greedy approach first proposed in the paper *GreedyMini: Generating low-density DNA minimizers*. 

## Table of Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Installation and Compilation](#installation-and-compilation)
- [Usage](#usage)
  - [GM-expected: Generating Minimizers For Expected Density](#gm-expected-generating-minimizers-for-expected-density)
  - [GM-particular: Generating Minimizers For Particular Density](#gm-particular-generating-minimizers-for-particular-density)
  - [GM-improve: Generating Minimizers for Large W](#gm-improve-generating-minimizers-for-large-w)
  - [GM-k: Generating Minimizers for Large K](#gm-k-generating-minimizers-for-large-k)
- [Accessing the Minimizers](#accessing-the-minimizers)
  - [Locating the Minimizers](#locating-the-minimizers)
  - [Exporting the Minimizers](#exporting-the-minimizers)
  - [Loading the Minimizers to Memory](#loading-the-minimizers-to-memory)
- [Miscellaneous](#miscellaneous)
  - [Compiling GreedyMini Manually](#compiling-greedymini-manually-1)
  - [Running Paper Tests](#running-paper-tests)
- [Contact](#contact)

## Introduction

All methods described here generate a *binary* minimizer order, in order to transform the binary minimizer to a DNA minimizer, we first encode the DNA sequence into binary, we then apply the GreedyMini order on the odd bits, apply a lexicographical order on the even bits (or any other order) and then concatenate the results. The upper bounds for DNA density shown in this toolkit are for the lowest density order between extending a binary minimizer to DNA using odd bits for GreedyMini and even for lexicographic and extending a binary minimizer to DNA using even bits for GreedyMini and odd bits for lexicographic. Note that in both cases the most significant bits come from GreedyMini.

<img src="github%20figures/example_binary_to_dna.png" alt="Extending to DNA" width="600">

To generate minimizers for various w and k values, we reccomend first running GM-expected to generate a minimizer for the set of required w and k values. In case GM-expected would take too much time (i.e., w+k is large), one can create an order for smaller set of w or k, and then extend them using GM-improve or GM-k, respectively. Usually it's better to run for smaller w and run GM-improve. 

Another option, is to generate an order for k' < k, and then rank the first k' bits using GreedyMini and rank the rest lexicographically (Before extending to DNA), very similarly to how we concatenate the ranks in the figure.  One can also generate an order for w' < w (as it would still rank all k-mers) and use that. If time allows it's always better to either generate an order from GM-expected, or GM-improve and use the aforementioned methods only if necessary.


## Prerequisites

Before you begin, ensure you have the following:

- **64-bit system**: Required to run the binaries or compile GreedyMini.
- **C++ Compiler**: Only required if you are compiling GreedyMini yourself (supports C++20, e.g., `g++` version 10 or higher).

## Installation and Compilation

Precompiled binaries are available for **Ubuntu** and **macOS**. You can download them from the [GitHub release page](https://github.com/OrensteinLab/GreedyMini/releases).
If you are using **Ubuntu**, ensure that the program has execution permissions by running the following command:  
```bash
chmod +x ./GreedyMini
```

To ensure that the compiled Linux binaries do not depend on newer GLIBC or libstdc++ symbols (e.g., `GLIBC_2.38`, `GLIBCXX_3.4.31`), we build on Ubuntu 20.04. This ensures that our released executables can run on older Linux distributions without encountering linker errors about missing GLIBC/GLIBCXX versions.

If your system is significantly older (e.g., it has GLIBC < 2.31), you may still need to build from source. Alternatively, you can run in a newer environment such as Docker or a more up-to-date Linux distribution.

Additionally, we provide instructions for manual compilation in the **[Compiling GreedyMini Manually](#compiling-greedymini-manually)** section.


## Usage

After compiling or downloading the binary, you can run **GreedyMini** using different modes and parameters.


### GM-expected: Generating Minimizers For Expected Density

To generate a minimizer with low expected density, run:

```bash
./GreedyMini -mode expected -w {w} -k {k}
```

Where:

- `{w}`: The window size.
- `{k}`: The k-mer size.

An example run would be:
```
./GreedyMini -mode expected -w 5 -k 4
```

#### Additional Parameters

- `-greedyE_runs`: Number of runs of GreedyE (default: `4096`)
- `-n_cores`: Number of CPU cores to use (default: half the total available threads due to hyper-threading)
- `-min_alpha`: Minimum alpha value (default: `0.939088`)
- `-max_alpha`: Maximum alpha value (default: `0.999590`)
- `-max_swapper_time_minutes`: Maximum SwapDFS time in minutes

### GM-particular: Generating Minimizers For Particular Density

#### Sequence format

`GreedyMini` expects a `.fasta` file containing exactly one sequence.

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

An example run would be:
```
./GreedyMini -mode particular -w 5 -k 4 -path chr_x_1m.fasta -name 1M
```
#### Additional Parameters

- `-greedyP_runs`: Number of runs of GreedyP (default: `4096`)
- `-n_cores`: Number of CPU cores to use (default: half the total available threads due to hyper-threading)
- `-min_alpha`: Minimum alpha value (default: `0.939088`)
- `-max_alpha`: Maximum alpha value (default: `0.999590`)

### GM-improve: Generating Minimizers for Large W
To generate a minimizer for large w (from a previous starting point), or to improve an existing minimizer, run:

```bash
./GreedyMini -mode improve -path {path} -w {w} -k {k} -max_swapper_time_minutes {max_swapper_time_minutes}
```

Where:
- `{path}`: A path to a GreedyMini order ('*.gm').
- `{w}`: The window size for the output minimizer (Not necessarily the w for which it was originally created for)
- `{k}`: The k-mer size.
- `-max_swapper_time_minutes`: Maximum SwapDP time in minutes

An example run would be:
```bash
./GreedyMini -mode improve -path output/minimizers/w5_k4.gm -w 30 -k 4 -max_swapper_time_minutes 5
```
#### Additional Parameters

- `-n_cores`: Number of CPU cores to use (default: half the total available threads due to hyper-threading)


### GM-k: Generating Minimizers for Large K

To generate a minimizer for large k (from a previous starting point), run:

```bash
./GreedyMini -mode extend_k -path {path} -w {w} -k {k} -k_extended {k_extended} -max_swapper_time_minutes {max_swapper_time_minutes}
```

Where:
- `{path}`: A path to a GreedyMini order ('*.gm').
- `{w}`: The window size.
- `{k}`: The k-mer size of the original order.
- `{k_extended}`: The k-mer size of the output order.
- `-max_swapper_time_minutes`: Maximum SwapDFS time in minutes - per each increase in k

An example run would be:
```bash
./GreedyMini -mode extend_k -path output/minimizers/w5_k4.gm -w 5 -k 4 -k_extended 8 -max_swapper_time_minutes 5
```

#### Additional Parameters

- `-n_cores`: Number of CPU cores to use (default: half the total available threads due to hyper-threading)


### Using the minimizers

## Accessing the Minimizers

### Locating the Minimizers

Generated minimizers will appear inside the `output/minimizers` folder. For particular density minimizers, they will appear in a subfolder with the selected name.

### Exporting the Minimizers

To export the minimizers to a `.csv` or `.txt` format, run:

```bash
./GreedyMini -mode export -path {path} -output_format {output_format}
```

Where:
- `{path}`: A path to a GreedyMini order ('*.gm').
- `-output_format`: Either 'csv' or 'txt'.

An example run would be:
```bash
./GreedyMini -mode export -path output/minimizers/w5_k4.gm -output_format txt
```


### Loading the Minimizers to Memory

We provide the Python notebook `load_order.ipynb` alongside the best minimizer orders from the paper, both located in the folder `minimizer loading example`. The notebook showcases how to load a minimizer to memory and print the order of each k-mer. For C++, we recommend looking at the functions `load_order()` and `load_vector_from_file()` located in `code/tools.cpp`.


## Misc.

### Compiling GreedyMini Manually

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

### Running Paper Tests

We ran our test on the first 1M nucleotides of chromosome X from [Genome assembly T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/). To do that we used the python notebook `shorten_fasta.ipynb` which is located in `various scripts/preprocessing chr x/`. We then put the resulting `.fasta` file in the same directory as the `GreedyMini` executable.

Execute the following command to run most of the tests from the paper:

```bash
./GreedyMini -mode tests
```



## Contact
In case of issues with GreedyMini, you may contact us at tziony.i@gmail.com.
