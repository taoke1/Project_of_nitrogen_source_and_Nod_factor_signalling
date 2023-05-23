#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=5g
#SBATCH -c 20
#SBATCH --time=03:00:0
#SBATCH --account InRoot

fastspar --threshold 0.2 --otu_table Filtered_tables/root_starv.tsv --correlation Correlations/root_starv_cor.tsv --covariance Covariances/root_starv_cov.tsv

mkdir bootstrap_counts
fastspar_bootstrap --otu_table Filtered_tables/root_starv.tsv --number 20000 --prefix bootstrap_counts/root_starv --threads 20

mkdir bootstrap_correlation
parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

fastspar_pvalues --otu_table Filtered_tables/root_starv.tsv --correlation Correlations/root_starv_cor.tsv --prefix bootstrap_correlation/cor_root_starv_ --permutations 20000 --outfile p_vals/pvalues_root_starv.tsv --threads 20

rm -r bootstrap_correlation/
rm -r bootstrap_counts/
