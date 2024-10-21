
algorithms=("mod-minimizer" "lr-minimizer" "minimizer" "miniception" "open-closed-syncmer" "rot-minimizer-orig" "decycling" "double-decycling")


k_values=$(seq 3 30)
w_values=(12)

# Adding fixed k and varying w case
fixed_k=8
varying_w_values=$(seq 3 30)

for file in sequence_100k.fasta sequence_1M.fasta 
do
    # Create a directory for each sequence file
    dirname=$(basename "$file" .fasta)
    mkdir -p "$dirname"

    for algo in "${algorithms[@]}"
    do
        # Prepare the output CSV file for this algorithm
        output_file="$dirname/${algo}_densities_runs_over_k_and_w.csv"

        # Create an associative array to hold the densities
        declare -A density_matrix_k_variable
        declare -A density_matrix_w_variable

        # First Case: Varying k, fixed w
        for w in "${w_values[@]}"
        do
            for k in $k_values
            do
                # Run your command and capture the output
                output=$(./density -i "$file" -k "$k" -w "$w" -a "$algo" --stream 2>&1)

                # Extract the line that contains the density value
                density_line=$(echo "$output" | grep "^$w,$k," | head -n 1)

                # Extract the density value (assuming it's the fourth field)
                density=$(echo "$density_line" | awk -F',' '{print $4}')

                # If density is empty, set it to NA
                if [ -z "$density" ]; then
                    density="NA"
                fi

                # Store the density in the matrix with keys as "w,k"
                density_matrix_k_variable["$w,$k"]="$density"

                # Output the information for debugging
                echo "$algo, $w, $k, $density"
            done
        done

        # Second Case: Varying w, fixed k
        for w in $varying_w_values
        do
            # Run your command and capture the output
            output=$(./density -i "$file" -k "$fixed_k" -w "$w" -a "$algo" --stream 2>&1)

            # Extract the line that contains the density value
            density_line=$(echo "$output" | grep "^$w,$fixed_k," | head -n 1)

            # Extract the density value (assuming it's the fourth field)
            density=$(echo "$density_line" | awk -F',' '{print $4}')

            # If density is empty, set it to NA
            if [ -z "$density" ]; then
                density="NA"
            fi

            # Store the density in the matrix with keys as "w"
            density_matrix_w_variable["$w"]="$density"

            # Output the information for debugging
            echo "$algo, $fixed_k, $w, $density"
        done

        # Output the CSV file for the first case (varying k, fixed w)
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
                    density="${density_matrix_k_variable["$w,$k"]}"
                    echo -n ",$density"
                done
                echo
            done
        } >> "$output_file"

        # Output the CSV file for the second case (varying w, fixed k)
        {
            # Print the header row (w values)
            echo -n "k/w"
            for w in $varying_w_values; do
                echo -n ",$w"
            done
            echo

            # Print the density values for fixed k
            echo -n "$fixed_k"
            for w in $varying_w_values
            do
                density="${density_matrix_w_variable["$w"]}"
                echo -n ",$density"
            done
            echo
        } >> "$output_file"

        # Clean up the associative arrays
        unset density_matrix_k_variable
        unset density_matrix_w_variable
    done
done
