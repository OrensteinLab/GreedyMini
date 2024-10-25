#!/bin/bash

# Log file to store all output
LOG_FILE="my_program_output_particular_v2.log"

# Clear the log file if it already exists
> "$LOG_FILE"

# Loop through all combinations of k and w
for w in $(seq 3 15); do
  for k in $(seq 3 15); do
    echo "Running: ./my_program -mode particular -w $w -k $k -path sequences_1M -name 1M" >> "$LOG_FILE"

    # Run the program with 'time' to capture resource usage
    /usr/bin/time -v ./my_program -mode particular -w "$w" -k "$k" -path sequences_1M -name 1M \
      >> "$LOG_FILE" 2>&1
      
    # Log that the combination has completed
    echo "Completed: w=$w, k=$k" >> "$LOG_FILE"
  done
done

echo "All processes completed. Check $LOG_FILE for output."