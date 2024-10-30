#!/bin/bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=syv3ap@virginia.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=8000
#SBATCH --time=1-00:00:00
#SBATCH -A ratan

# Code starts here
# set -xe # Enable both debugging and exit-on-error

# Directory containing the chromosome info files
data_dir="/standard/vol169/cphg_ratan/ar7jq/TCGA/imputed_germline"

# Get the directory this script is in and set it as base_dir
base_dir="/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline"

# Loop through chromosomes 1 to 22 (adjust if you need more chromosomes)
for chr_num in {1..22}; do
    # Set the current chromosome file name
    chr_file="${data_dir}/chr${chr_num}.info.gz"

    # Get the base pair (bp) from the last row of the current chromosome file
    last_bp=$(zcat "$chr_file" | tail -n 1 | awk -F':' '{print $2}')

    # Calculate the total number of 10 Mb windows for this chromosome
    total_windows=$(( (last_bp / 10000000) + 1 ))

    # Generate options#.txt for the current chromosome (option1.txt, option2.txt, etc.)
    option_file="${base_dir}/option${chr_num}.txt"
    mkdir -p "$base_dir"  # Ensure the base directory exists

    # Save total_windows to a file for later use
    echo "$total_windows" > "${base_dir}/windows_chr${chr_num}.txt"

    # Generate the window ranges and save them to option#.txt
    for i in $(seq 0 10000000 $((total_windows * 10000000))); do
        echo "${i}-$((i + 10000000 - 1))"
    done > "$option_file"

    echo "Generated $option_file for chr${chr_num} with $total_windows windows."
done

