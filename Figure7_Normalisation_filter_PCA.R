# Code by Ib Thorsgaard Jensen
library(data.table)
library(ggplot2)
library(lmerTest)
library(ggfortify)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Changes working directory to the path of this file. Only works in Rstudio

# Loading data
Metabolite_raw <- fread("../Data/FT_MS2_quant_remove_last_run.csv")
meta_data <- fread("../Data/KeTao_root_exudate_sample_list.csv")
colnames(Metabolite_raw)[1] <- "FeatureID"
Metabolite_raw[,FeatureID:=paste0("Feature", FeatureID)]
Metabolite_raw[,c("row m/z", "row retention time"):=NULL]

# Samples to omit due to experimental concerns
omit_samples <- c("Gifu_1_4", "Gifu_1_6", "Gifu_3_3", "C2_3_3")

# Reorganizing the meta data
meta_data <- meta_data[Priority %in% c("Yes", "yes")]
meta_data[,Sample_ID:=gsub("-", "_", Sample_ID)]
meta_data[,Sample_ID:=gsub(" ", "_", Sample_ID)]

# Constructing the state variable
meta_data[,State:=fcase(
  R7A == "no" & Nitrate == "no", "Starved",
  R7A == "yes" & Nitrate == "no", "Symbiotic",
  R7A == "no" & Nitrate == "yes", "Inorganic",
  R7A == "yes" & Nitrate == "yes", "Inorganic+R7A"
)]


# Including only the measurements of peak area
sample_names <- colnames(Metabolite_raw)[grepl("Peak area", colnames(Metabolite_raw))]

# Saving only the measurements of peak area in data table
Metabolite_peak_area <- Metabolite_raw[,..sample_names]

# Simplifying samples names
C2 <- substr(sample_names[!grepl("Gifu", sample_names)], 8, 13)
g <- substr(sample_names[grepl("Gifu", sample_names)], 8, 15)
sample_names[!grepl("Gifu", sample_names)] <- C2
sample_names[grepl("Gifu", sample_names)] <- g
colnames(Metabolite_peak_area) <- sample_names

# Adding the sampletype variable
meta_data[,Sample_type:=rep(c("Gifu","Control"),c(24,20))]

# Transposing metabolite data and adding feature names
Metabolite_data <- transpose(Metabolite_peak_area, keep.names = "Sample_ID")
colnames(Metabolite_data)[-1] <- Metabolite_raw$FeatureID

# Removing features with zeros in more than 90 % of samples
Absent <- apply(Metabolite_data[,-1], 2, function(x) mean(x==0)>0.9)
Present <- c("Sample_ID",names(Absent)[!Absent])
Metabolite_data <- Metabolite_data[,..Present]

# Combining metabolite data with meta data and removing poor samples
D <- merge(meta_data[,c(1:2,6:7)], Metabolite_data, by = "Sample_ID")
D <- D[!(Sample_ID %in% omit_samples)]

# Normalization - technique similar to the one used by default in DESeq2
Feature_table <- D[,-(1:4)]
non_zero <- apply(Feature_table, 2, function(x) all(x!=0))
non_zero <- names(non_zero)[non_zero]
D_mat <- as.matrix(Feature_table[,..non_zero])
geom_means <- apply(D_mat, 2, function(x) exp(mean(log(x))))
norm_fac <- apply(D_mat %*% diag(1/geom_means), 1, median)
D_norm <- as.matrix(D[,-(1:4)]*norm_fac)

# Log-transforming the data
D_norm_log <- D_norm
D_norm_log[D_norm_log == 0] <- 0.5 #Replacing zeros with small values prior to log transformation
D_norm_log <- log(D_norm_log)
D_norm_log <- data.table(D[,1:4], D_norm_log)

# Tables with colors
col.table <- data.table(gt = c("Inorganic", "Starved", "Symbiotic", "Inorganic+R7A"),
                        color = c("#CC79A7", "#0072B2", "#009E73", "#D95F02") )
# PCA - all samples
p <- prcomp(D_norm_log[,-(1:4)], center = T, scale. = T)
PCA_all <- autoplot(p, data = D_norm_log, colour = "State", shape = "Sample_type", size = 4) + 
  scale_color_manual(name = "State",
                     breaks = col.table$gt,
                     values = col.table$color)

# PCA - only gifu samples
D_gifu <- D_norm_log[Sample_type == "Gifu"]
all_zero_gifu <- apply(D_gifu[,-c(1:4)], 2, function(x) all(x==log(0.5)))
all_zero_gifu <- names(all_zero_gifu)[!all_zero_gifu]
col_names_keep <- c(colnames(D)[1:4], all_zero_gifu)
p2 <- prcomp(D_gifu[,..all_zero_gifu], center = T, scale. = T)
PCA_gifu <- autoplot(p2, data = D_gifu[,..col_names_keep], colour = "State", size = 4) + 
  scale_color_manual(name = "State",
                     breaks = col.table$gt,
                     values = col.table$color)

# Saving plots
ggsave(PCA_all, width = 230, height = 140, unit = "mm", file = "Figures/PCA/Metabolite_PCA_full.pdf")
ggsave(PCA_gifu, width = 230, height = 140, unit = "mm", file = "Figures/PCA/Metabolite_PCA_Gifu_full.pdf")

# Identifying all metabolites where the state-effect differ between the control and gifu samples
p_vals <- rep(NA, ncol(D_norm_log)-4)
for(i in 1:(ncol(D_norm_log)-4)){
  cat(i,"\r")
  j <- i+4
  peak <- D_norm_log[[j]]
  l <- lmer(peak ~ State*Sample_type + (1|Sample_ID), data = D_norm_log)
  al <- anova(l)
  p_vals[i] <- al["State:Sample_type","Pr(>F)"]
}
p_adj <- p.adjust(p_vals, method = "BY")
sig_feat <- colnames(D_norm_log)[-(1:4)][p_adj<0.05]
keep <- c(colnames(D)[1:4], sig_feat)
D_gifu_filt <- D_norm_log[Sample_type == "Gifu",..keep]
D_all_filt <- D_norm_log[,..keep]

# Removing all metabolite that remain, but are never detected in the gifu samples
all_zero <- apply(D_gifu_filt[,-(1:4)], 2, function(x) all(x == 0))
keep <- colnames(D_gifu_filt)[c(T,T,T,T,!all_zero)]
D_gifu_filt <- D_gifu_filt[,..keep]

# Constructing tablw after normalization and filtering
D_norm2 <- D_norm[,colnames(D_norm) %in% colnames(D_gifu_filt)[-(1:4)]]
D_norm2 <- data.table(D[,c(1,3,4)], D_norm2)
fwrite(D_norm2, "../Data/Metabolites_filtered_lmer.csv")

# PCA - only Gifu samples after filtering metabolites
p3 <- prcomp(D_gifu_filt[,-(1:4)], center = T, scale. = T)
PCA_gifu2 <- autoplot(p3, data = D_gifu_filt, colour = "State", size = 4) + 
  scale_color_manual(name = "State",
                     breaks = col.table$gt,
                     values = col.table$color)
ggsave(PCA_gifu2, width = 230, height = 140, unit = "mm", file = "Figures/PCA/Metabolite_PCA_Gifu_filter.pdf")

# PCA - all samples after filtering metabolites
p4 <- prcomp(D_all_filt[,-(1:4)], center = T, scale. = T)
PCA_all2 <- autoplot(p4, data = D_all_filt, colour = "State", shape = "Sample_type", size = 4) + 
  scale_color_manual(name = "State",
                     breaks = col.table$gt,
                     values = col.table$color)
ggsave(PCA_all2, width = 230, height = 140, unit = "mm", file = "Figures/PCA/Metabolite_PCA_all_filter.pdf")
