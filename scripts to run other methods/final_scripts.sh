#!/bin/bash


algorithms=("mod-minimizer" "lr-minimizer" "minimizer" "miniception" "open-closed-syncmer" "rot-minimizer-orig" "decycling" "double-decycling")


# First script section (k x w matrix)
k_values=$(seq 3 30)
w_values=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)

for file in random.10M.fa  
do
    # Create a directory for each sequence file
    dirname=$(basename "$file" .fa)
    mkdir -p "$dirname"

    for algo in "${algorithms[@]}"
    do
        # Prepare the output CSV file for this algorithm
        output_file="$dirname/${algo}_densities_runs_over_k.csv"

        # Create an associative array to hold the densities
        declare -A density_matrix

        # Collect density values
        for w in "${w_values[@]}"
        #for w in $w_values
        do
            for k in $k_values
            do
                # Run your command and capture the output
                output=$(./density -i "$file" -k "$k" -w "$w" -a "$algo" --stream 2>&1)

                # Extract the line that contains the density value
                # Assuming the density line starts with 'w,k,'
                density_line=$(echo "$output" | grep "^$w,$k," | head -n 1)

                # Extract the density value (assuming it's the fourth field)
                density=$(echo "$density_line" | awk -F',' '{print $4}')

                # If density is empty, set it to NA
                if [ -z "$density" ]; then
                    density="NA"
                fi

                # Store the density in the matrix with keys as "w,k"
                density_matrix["$w,$k"]="$density"

                # Output the information for debugging
                echo "$algo, $w, $k, $density"
            done
        done

        # Output the CSV file
        {
            # Print the header row (k values)
            echo -n "w/k"
            for k in $k_values; do
                echo -n ",$k"
            done
            echo

            # Print the density values for each w
            for w in "${w_values[@]}"
            do
                echo -n "$w"
                for k in $k_values
                do
                    density="${density_matrix["$w,$k"]}"
                    echo -n ",$density"
                done
                echo
            done
        } > "$output_file"

        # Clean up the associative array
        unset density_matrix
    done
done

# Second script section (w=k values)
wk_values=$(seq 2 20)
for file in random.10M.fa 
do
    # Create a directory for each sequence file
    dirname=$(basename "$file" .fa)
    mkdir -p "$dirname"

    for algo in "${algorithms[@]}"
    do
        # Prepare the output CSV file for this algorithm
        output_file="$dirname/${algo}_wk.csv"

        # Start the CSV file with headers
        echo "w=k,density" > "$output_file"

        for wk in $wk_values
        do
            # Run the density command and capture the output
            output=$(./density -i "$file" -k "$wk" -w "$wk" -a "$algo" --stream 2>&1)

            # Extract the density line that matches w=k=wk
            density_line=$(echo "$output" | grep "^$wk,$wk," | head -n 1)

            # Extract the density value (assuming it's the fourth field)
            density=$(echo "$density_line" | awk -F',' '{print $4}')

            # If density is empty, set it to NA
            if [ -z "$density" ]; then
                density="NA"
            fi

            # Append the w=k and density to the CSV file
            echo "$wk,$density" >> "$output_file"

            # Output the information for debugging
            echo "$algo, $wk, $wk, $density"
        done
    done
done


# Third script section (w=k-1 values)
wk_values=$(seq 3 20)
for file in random.10M.fa 
do
    # Create a directory for each sequence file
    dirname=$(basename "$file" .fa)
    mkdir -p "$dirname"

    for algo in "${algorithms[@]}"
    do
        # Prepare the output CSV file for this algorithm
        output_file="$dirname/${algo}_w_k_minus_1.csv"

        # Start the CSV file with headers
        echo "k,w,density" > "$output_file"

        for wk in $wk_values
        do
            w=$((wk - 1))

            # Run the density command and capture the output
            output=$(./density -i "$file" -k "$wk" -w "$w" -a "$algo" --stream 2>&1)

            # Debug: print full output for w=k-1 cases
            #echo "Full density command output for algo=$algo, k=$wk, w=$w:"
            #echo "$output"

            # Extract the density line that matches k=wk, w=wk-1
            density_line=$(echo "$output" | grep -E "^ *$w *,$wk *," | head -n 1)

            # Debug: print extracted density line
            #echo "Extracted density line for k=$wk, w=$w: $density_line"

            # Extract the density value (assuming it's the fourth field)
            density=$(echo "$density_line" | awk -F',' '{print $4}')

            # If density is empty, set it to NA
            if [ -z "$density" ]; then
                density="NA"
            fi

            # Append k, w, and density to the CSV file
            echo "$wk,$w,$density" >> "$output_file"

            # Output the information for debugging
            echo "$algo, k=$wk, w=$w, density=$density"
        done
    done
done







# Fourth script section (wide w values)
w_values=$(seq 3 201)
k_values=(8 12)

for file in random.10M.fa 
do
    # Create a directory for each sequence file
    dirname=$(basename "$file" .fa)
    mkdir -p "$dirname"

    for algo in "${algorithms[@]}"
    do
        # Prepare the output CSV file for this algorithm
        output_file="$dirname/${algo}_wide_w.csv"

        # Create an associative array to hold the densities
        declare -A density_matrix

        # Collect density values
        for w in $w_values
        do
            for k in "${k_values[@]}"
            do
                # Run your command and capture the output
                output=$(./density -i "$file" -k "$k" -w "$w" -a "$algo" --stream 2>&1)

                # Extract the density line that matches w and k
                density_line=$(echo "$output" | grep "^$w,$k," | head -n 1)

                # Extract the density value (assuming it's the fourth field)
                density=$(echo "$density_line" | awk -F',' '{print $4}')

                # If density is empty, set it to NA
                if [ -z "$density" ]; then
                    density="NA"
                fi

                # Store the density in the matrix with keys as "w,k"
                density_matrix["$w,$k"]="$density"

                # Output the information for debugging
                echo "$algo, $w, $k, $density"
            done
        done

        # Output the CSV file
        {
            # Print the header row (k values)
            echo -n "w/k"
            for k in "${k_values[@]}"; do
                echo -n ",$k"
            done
            echo

            # Print the density values for each w
            for w in $w_values
            do
                echo -n "$w"
                for k in "${k_values[@]}"
                do
                    density="${density_matrix["$w,$k"]}"
                    echo -n ",$density"
                done
                echo
            done
        } > "$output_file"

        # Clean up the associative array
        unset density_matrix
    done
done


# Fifth script section (wide k values)
w_values=(8 12)
k_values=($(seq 2 1 40))

for file in random.10M.fa 
do
    # Create a directory for each sequence file
    dirname=$(basename "$file" .fa)
    mkdir -p "$dirname"

    for algo in "${algorithms[@]}"
    do
        # Prepare the output CSV file for this algorithm
        output_file="$dirname/${algo}_wide_k.csv"

        # Create an associative array to hold the densities
        declare -A density_matrix

        # Collect density values
        for w in "${w_values[@]}"
        do
            for k in "${k_values[@]}"
            do
                # Run your command and capture the output
                output=$(./density -i "$file" -k "$k" -w "$w" -a "$algo" --stream 2>&1)

                # Extract the density line that matches w and k
                density_line=$(echo "$output" | grep "^$w,$k," | head -n 1)

                # Extract the density value (assuming it's the fourth field)
                density=$(echo "$density_line" | awk -F',' '{print $4}')

                # If density is empty, set it to NA
                if [ -z "$density" ]; then
                    density="NA"
                fi

                # Store the density in the matrix with keys as "w,k"
                density_matrix["$w,$k"]="$density"

                # Output the information for debugging
                echo "$algo, $w, $k, $density"
            done
        done

        # Output the CSV file
        {
            # Print the header row (k values)
            echo -n "w/k"
            for k in "${k_values[@]}"; do
                echo -n ",$k"
            done
            echo

            # Print the density values for each w
            for w in "${w_values[@]}"
            do
                echo -n "$w"
                for k in "${k_values[@]}"
                do
                    density="${density_matrix["$w,$k"]}"
                    echo -n ",$density"
                done
                echo
            done
        } > "$output_file"

        # Clean up the associative array
        unset density_matrix
    done
done
