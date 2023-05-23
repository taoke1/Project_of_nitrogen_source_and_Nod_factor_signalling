#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=5g
#SBATCH -c 20
#SBATCH --time=05:00:0
#SBATCH --account InRoot

fastspar --threshold 0.2 --otu_table Filtered_tables/root_inog.tsv --correlation Correlations/root_inog_cor.tsv --covariance Covariances/root_inog_cov.tsv

mkdir bootstrap_counts
fastspar_bootstrap --otu_table Filtered_tables/root_inog.tsv --number 50000 --prefix bootstrap_counts/root_inog --threads 20

mkdir bootstrap_correlation
parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

fastspar_pvalues --otu_table Filtered_tables/root_inog.tsv --correlation Correlations/root_inog_cor.tsv --prefix bootstrap_correlation/cor_root_inog_ --permutations 50000 --outfile p_vals/pvalues_root_inog.tsv --threads 20

rm -r bootstrap_correlation/
rm -r bootstrap_counts/
