#!/bin/bash

# Log file to store all output
LOG_FILE="my_program_output_particular.log"

# Clear the log file if it already exists
> "$LOG_FILE"

# First: w = 12 and 3 <= k <= 16
w=12
for k in $(seq 3 16); do
  echo "Running: ./GreedyMini-ubuntu -mode particular -w $w -k $k -path chr_x_1m.fasta -name chrx_1M" >> "$LOG_FILE"

  # Run the program with 'time' to capture resource usage
  /usr/bin/time -v ./GreedyMini-ubuntu -mode particular -w "$w" -k "$k" -path chr_x_1m.fasta -name 1M \
    >> "$LOG_FILE" 2>&1

  # Log that the combination has completed
  echo "Completed: w=$w, k=$k" >> "$LOG_FILE"
done

# Second: k = 8 and 3 <= w <= 19
k=8
for w in $(seq 3 19); do
  echo "Running: ./GreedyMini-ubuntu -mode particular -w $w -k $k -path chr_x_1m.fasta -name 1M" >> "$LOG_FILE"

  # Run the program with 'time' to capture resource usage
  /usr/bin/time -v ./GreedyMini-ubuntu -mode particular -w "$w" -k "$k" -path chr_x_1m.fasta -name 1M \
    >> "$LOG_FILE" 2>&1

  # Log that the combination has completed
  echo "Completed: w=$w, k=$k" >> "$LOG_FILE"
done

echo "All processes completed. Check $LOG_FILE for output."
