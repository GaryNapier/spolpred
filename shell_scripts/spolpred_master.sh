#! /bin/bash

# FILENAME
# spolpred_master.sh

set -e
set -u
set -o pipefail

# ------
# Setup 
# ------

cd ~/spolpred

# ----------
# Variables
# ----------

# ------------
# Directories
# ------------

# Existing directories
# tbp_results_dir=/mnt/storage7/jody/tb_ena/tbprofiler/freebayes/results/
# tbp_results_dir=/mnt/storage7/jody/tb_ena/tbprofiler/vcf_profile/results/
tbprofiler_results_dir=/mnt/storage7/jody/tb_ena/tbprofiler/gary/results/
metadata_remote_dir=../metadata/
metadata_local_dir=metadata/
database_dir=../pipeline/db/
tbdb_dir=../tbdb/
results_dir=results/
fasta_dir=fasta/
newick_dir=${results_dir}newick/
vcf_remote=/mnt/storage7/jody/tb_ena/per_sample/
vcf_db=~/vcf/
vcf_local_dir=vcf/


# ------
# Files
# ------

# Existing files

# Created files
# get_lineages.py
lineage_file=${results_dir}lineage_file.csv
# samples_by_lineage.R
all_lineages_samples_outfile=${results_dir}all_lineages_samples.txt
grouped_samples_prefix=${results_dir}grouped_samples_

# ------------
# Run scripts
# ------------

# Pull all lineage data from TB-profiler
if [ ! -f ${lineage_file} ]; then
    python python_scripts/get_lineages.py \
    --tbprofiler-results-dir ${tbprofiler_results_dir} \
    --lineage-file ${lineage_file}
fi

# Subset the lineages and make sample lists for trees
Rscript r_scripts/samples_by_lineage.R \
--lineage_file ${lineage_file} \
--all_lineages_samples_outfile ${all_lineages_samples_outfile}

# Run tree for all samples
tree_pipeline.sh all_lineages ${vcf_db} ${vcf_local_dir} ${all_lineages_samples_outfile} ${fasta_dir} ${newick_dir}

# Run trees for grouped lineages
for f in ${grouped_samples_prefix}*; do 
    code=$(basename ${f} .txt)
    tree_pipeline.sh ${code} ${vcf_db} ${vcf_local_dir} ${f} ${fasta_dir} ${newick_dir}
done

