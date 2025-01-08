#!/bin/bash

# Define algorithms
#algorithms=("mod-minimizer" "lr-minimizer" "minimizer" "miniception" "open-closed-syncmer" "rot-minimizer-orig" "decycling" "double-decycling")
algorithms=("minimizer" "miniception" "open-closed-syncmer" "decycling" "double-decycling")

# Define specific (k, w) pairs as a list of tuples
k_w_pairs=(
    "15,17"
    "31,5"
    "7,22"
    "7,49"
    "9,16"
    "21,11"
    "15,10"
    "29,11"
    "19,30"
)

# Input file
file="random.10M.fa"

# Create a directory for the sequence file
dirname=$(basename "$file" .fa)
mkdir -p "$dirname"

# Loop through algorithms
for algo in "${algorithms[@]}"
do
    # Prepare the output CSV file for this algorithm
    output_file="$dirname/${algo}_specific_k_w_pairs.csv"

    # Start the CSV file with headers
    echo "k,w,density" > "$output_file"

    # Process each (k, w) pair
    for pair in "${k_w_pairs[@]}"
    do
        IFS="," read -r k w <<< "$pair"  # Split the pair into k and w

        # Run the density command and capture the output
        output=$(./density -i "$file" -k "$k" -w "$w" -a "$algo" --stream 2>&1)

        # Extract the density line that matches k and w
        density_line=$(echo "$output" | grep "^$w,$k," | head -n 1)

        # Extract the density value (assuming it's the fourth field)
        density=$(echo "$density_line" | awk -F',' '{print $4}')

        # If density is empty, set it to NA
        if [ -z "$density" ]; then
            density="NA"
        fi

        # Append k, w, and density to the CSV file
        echo "$k,$w,$density" >> "$output_file"

        # Output the information for debugging
        echo "$algo, k=$k, w=$w, density=$density"
    done

done
