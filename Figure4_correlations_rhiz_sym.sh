#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem=5g
#SBATCH -c 20
#SBATCH --time=02:30:0
#SBATCH --account InRoot

fastspar --threshold 0.2 --otu_table Filtered_tables/rhiz_sym.tsv --correlation Correlations/rhiz_sym_cor.tsv --covariance Covariances/rhiz_sym_cov.tsv

mkdir bootstrap_counts
fastspar_bootstrap --otu_table Filtered_tables/rhiz_sym.tsv --number 20000 --prefix bootstrap_counts/rhiz_sym --threads 20

mkdir bootstrap_correlation
parallel fastspar --otu_table {} --correlation bootstrap_correlation/cor_{/} --covariance bootstrap_correlation/cov_{/} -i 5 ::: bootstrap_counts/*

fastspar_pvalues --otu_table Filtered_tables/rhiz_sym.tsv --correlation Correlations/rhiz_sym_cor.tsv --prefix bootstrap_correlation/cor_rhiz_sym_ --permutations 20000 --outfile p_vals/pvalues_rhiz_sym.tsv --threads 20

rm -r bootstrap_correlation/
rm -r bootstrap_counts/
