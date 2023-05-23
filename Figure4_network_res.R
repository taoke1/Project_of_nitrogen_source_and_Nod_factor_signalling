# Code by Ib Thorsgaard Jensen
library(data.table)
library(magrittr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Changes working directory to the path of this file. Only works in Rstudio

## Rhizosphere ----
# p_val
P_rhiz_inog <- fread("p_vals/pvalues_rhiz_inog.tsv")
P_rhiz_sym <- fread("p_vals/pvalues_rhiz_sym.tsv")
P_rhiz_starv <- fread("p_vals/pvalues_rhiz_starv.tsv")

P_rhiz_inog <- as.matrix(P_rhiz_inog[,-1])
P_rhiz_sym <- as.matrix(P_rhiz_sym[,-1])
P_rhiz_starv <- as.matrix(P_rhiz_starv[,-1])

p_rhiz_inog_adj <- p.adjust(P_rhiz_inog[lower.tri(P_rhiz_inog)], method = "fdr", n = 571915)
p_rhiz_sym_adj <- p.adjust(P_rhiz_sym[lower.tri(P_rhiz_sym)], method = "fdr", n = 571915)
p_rhiz_starv_adj <- p.adjust(P_rhiz_starv[lower.tri(P_rhiz_starv)], method = "fdr", n = 571915)

sum(p_rhiz_inog_adj<0.1)
sum(p_rhiz_sym_adj<0.1)
sum(p_rhiz_starv_adj<0.1)
# correlations
cor_rhiz_inog <- fread("Correlations/rhiz_inog_cor.tsv")
cor_rhiz_sym <- fread("Correlations/rhiz_sym_cor.tsv")
cor_rhiz_starv <- fread("Correlations/rhiz_starv_cor.tsv")

C_rhiz_inog <- as.matrix(cor_rhiz_inog[,-1])
C_rhiz_sym <- as.matrix(cor_rhiz_sym[,-1])
C_rhiz_starv <- as.matrix(cor_rhiz_starv[,-1])

# Keep only significant correlations, inog
adj.mat <- matrix(0, nrow(P_rhiz_inog), ncol(P_rhiz_inog))
adj.mat[lower.tri(adj.mat)] <- p_rhiz_inog_adj
adj.mat <- adj.mat + t(adj.mat)
C_rhiz_inog[adj.mat>=0.1] <- 0
cor_rhiz_inog[,2:ncol(cor_rhiz_inog):= C_rhiz_inog %>% data.table()]

# Keep only significant correlations, sym
adj.mat <- matrix(0, nrow(P_rhiz_sym), ncol(P_rhiz_sym))
adj.mat[lower.tri(adj.mat)] <- p_rhiz_sym_adj
adj.mat <- adj.mat + t(adj.mat)
C_rhiz_sym[adj.mat>=0.1] <- 0
cor_rhiz_sym[,2:ncol(cor_rhiz_sym):= C_rhiz_sym %>% data.table()]

# Keep only significant correlations, starv
adj.mat <- matrix(0, nrow(P_rhiz_starv), ncol(P_rhiz_starv))
adj.mat[lower.tri(adj.mat)] <- p_rhiz_starv_adj
adj.mat <- adj.mat + t(adj.mat)
C_rhiz_starv[adj.mat>=0.1] <- 0
cor_rhiz_starv[,2:ncol(cor_rhiz_starv):= C_rhiz_starv %>% data.table()]

fwrite(cor_rhiz_inog, "Significant correlations/rhiz_inog.csv")
fwrite(cor_rhiz_sym, "Significant correlations/rhiz_sym.csv")
fwrite(cor_rhiz_starv, "Significant correlations/rhiz_starv.csv")

## Root ----
P_root_inog <- fread("p_vals/pvalues_root_inog.tsv")
P_root_sym <- fread("p_vals/pvalues_root_sym.tsv")
P_root_starv <- fread("p_vals/pvalues_root_starv.tsv")

P_root_inog <- as.matrix(P_root_inog[,-1])
P_root_sym <- as.matrix(P_root_sym[,-1])
P_root_starv <- as.matrix(P_root_starv[,-1])

p_root_inog_adj <- p.adjust(P_root_inog[lower.tri(P_root_inog)], method = "fdr", n = 590241)
p_root_sym_adj <- p.adjust(P_root_sym[lower.tri(P_root_sym)], method = "fdr", n = 590241)
p_root_starv_adj <- p.adjust(P_root_starv[lower.tri(P_root_starv)], method = "fdr", n = 590241)

sum(p_root_inog_adj<0.1)
sum(p_root_sym_adj<0.1)
sum(p_root_starv_adj<0.1)

# correlations
cor_root_inog <- fread("Correlations/root_inog_cor.tsv")
cor_root_sym <- fread("Correlations/root_sym_cor.tsv")
cor_root_starv <- fread("Correlations/root_starv_cor.tsv")

C_root_inog <- as.matrix(cor_root_inog[,-1])
C_root_sym <- as.matrix(cor_root_sym[,-1])
C_root_starv <- as.matrix(cor_root_starv[,-1])

# Keep only significant correlations, inog
adj.mat <- matrix(0, nrow(P_root_inog), ncol(P_root_inog))
adj.mat[lower.tri(adj.mat)] <- p_root_inog_adj
adj.mat <- adj.mat + t(adj.mat)
C_root_inog[adj.mat>=0.1] <- 0
cor_root_inog[,2:ncol(cor_root_inog):= C_root_inog %>% data.table()]

# Keep only significant correlations, sym
adj.mat <- matrix(0, nrow(P_root_sym), ncol(P_root_sym))
adj.mat[lower.tri(adj.mat)] <- p_root_sym_adj
adj.mat <- adj.mat + t(adj.mat)
C_root_sym[adj.mat>=0.1] <- 0
cor_root_sym[,2:ncol(cor_root_sym):= C_root_sym %>% data.table()]

# Keep only significant correlations, starv
adj.mat <- matrix(0, nrow(P_root_starv), ncol(P_root_starv))
adj.mat[lower.tri(adj.mat)] <- p_root_starv_adj
adj.mat <- adj.mat + t(adj.mat)
C_root_starv[adj.mat>=0.1] <- 0
cor_root_starv[,2:ncol(cor_root_starv):= C_root_starv %>% data.table()]

fwrite(cor_root_inog, "Significant correlations/root_inog.csv")
fwrite(cor_root_sym, "Significant correlations/root_sym.csv")
fwrite(cor_root_starv, "Significant correlations/root_starv.csv")

## Soil ----
P_soil_noN <- fread("p_vals/pvalues_soil_noN.tsv")
P_soil_N <- fread("p_vals/pvalues_soil_N.tsv")

P_soil_noN <- as.matrix(P_soil_noN[,-1])
P_soil_N <- as.matrix(P_soil_N[,-1])

P_soil_noN_adj <- p.adjust(P_soil_noN[lower.tri(P_soil_noN)], method = "fdr")
P_soil_N_adj <- p.adjust(P_soil_N[lower.tri(P_soil_N)], method = "fdr")

sum(P_soil_noN_adj<0.1) # Likely no statistically significant correlations
sum(P_soil_N_adj<0.1)
