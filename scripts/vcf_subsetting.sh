#!/bin/bash
#
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=syv3ap@virginia.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=8000
#SBATCH --time=1-00:00:00
#SBATCH -A ratan

# Manually set the number of windows (manual) (Check relevant windows_chr#.txt)
#SBATCH --array=1-7  # Number of windows for current chromosome

# Load bcftools module
module load bcftools

# Set the chromosome number dynamically (manual)
chr_num=20  # Select chromosome

# Base directory for generated files
base_dir="/standard/cphg-RLscratch/syv3ap/rotation/TCGA_LGG_imputed_germline"
chr_dir="${base_dir}/chr${chr_num}"

# Create the directory structure if it doesn't exist
mkdir -p "$chr_dir"

# Get the window range from the appropriate option#.txt file
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ${base_dir}/option${chr_num}.txt)

# Extract the start and end of the window from the options
START=$(echo $OPTS | cut -d'-' -f1)
END=$(echo $OPTS | cut -d'-' -f2)

# Step 1: Directly pipe filtered SNPs from chr#.info.gz into bcftools
zcat /standard/vol169/cphg_ratan/ar7jq/TCGA/imputed_germline/chr${chr_num}.info.gz | \
awk -v start="$START" -v end="$END" -v rsq_cutoff=0.3 -F'\t' 'NR > 1{
    split($1, a, ":");
    pos = a[2];
    if (pos >= start && pos <= end && $7 > rsq_cutoff) {
        print a[1] "\t" a[2]   # Output in two-column tab-delimited format
    }
}' > "${chr_dir}/chr${chr_num}_${SLURM_ARRAY_TASK_ID}.txt"

# Check if the SNP list is empty
if [ ! -s "${chr_dir}/chr${chr_num}_${SLURM_ARRAY_TASK_ID}.txt" ]; then
    echo "No SNPs in window $START-$END, skipping further processing for this window."
    exit 0
fi

# Create a unique temporary VCF file for this array task
tmp_vcf="${chr_dir}/chr${chr_num}_LGG_${SLURM_ARRAY_TASK_ID}.vcf.gz"

# Step 2: Select sample based on LGG TSS and compress with bgzip
bcftools view -S LGG_samples.txt \
/standard/vol169/cphg_ratan/ar7jq/TCGA/imputed_germline/chr${chr_num}.dose.vcf.gz --force-samples -Oz -o "$tmp_vcf"

# Index the temporary VCF file
bcftools index "$tmp_vcf"

# Step 3: Pipe LGG samples SNP list into bcftools for filtering VCF
bcftools view -R "${chr_dir}/chr${chr_num}_${SLURM_ARRAY_TASK_ID}.txt" \
"$tmp_vcf" > "${chr_dir}/chr${chr_num}_LGG_filtered_${SLURM_ARRAY_TASK_ID}.vcf.gz"

# Clean up the temporary files
rm "$tmp_vcf" "${tmp_vcf}.csi"

