# Code by Ib Thorsgaard Jensen
library(data.table)
# library(magrittr)
library(edgeR)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading datasets
OTU_table <- fread("../Data/zotu_table_0.01_clean.csv")
meta_data <- fread("../Data/design_clean.csv")

# Subsetting and filtering
comp <- "root"
Nitrate <- "yes"
samples <- meta_data[compartments == comp & genotype != "epr3" & Nitrate_supplied == Nitrate, SampleID]
meta_data_subset <- meta_data[SampleID %in% samples]; samples <- meta_data_subset$SampleID
meta_data_subset[,genotype:=factor(genotype, levels = c("Gifu", "nfre", "chit5", "nfr5"))]

data_subset <- OTU_table[,c("OTUid", samples), with = FALSE]

absent_otus <- apply(data_subset[,-1], 1, function(x) mean(x == 0)==1)
data_subset <- data_subset[!absent_otus]

## Differential abundance analysis with edgeR
feature.table <- data.frame(data_subset[,-1], row.names = data_subset$OTUid)
data.DGE <- DGEList(counts=feature.table, group=meta_data_subset$genotype)
data.DGE <- calcNormFactors(data.DGE)

design.matrix <- model.matrix(~ genotype, data = meta_data_subset)
data.DGE.disp <- estimateDisp(data.DGE, design.matrix)

neg.bin.reg <- glmQLFit(data.DGE.disp, design.matrix)
logFCranks <- apply(abs(neg.bin.reg$coefficients[,-1]), 2, rank, ties.method = "max")
logFCsigns <- apply(neg.bin.reg$coefficients[,-1], 2, sign)
D <- data.table(OTUid = rownames(logFCranks), logFCranks, logFCsigns)

# filtered <- fread("../Data/OTU_table_norm_subset_NGroot0.3.csv")
gifu_samples <- meta_data_subset[genotype == "Gifu", SampleID]
gifu_data <- as.matrix(data_subset[,..gifu_samples])
gifu_data <- t(t(gifu_data)/colSums(gifu_data))
abd_metric <- apply(gifu_data, 1, max)
OTU_keep <- data_subset$OTUid[abd_metric > 0.003]

D[OTUid %in% OTU_keep,lapply(.SD, function(x) sum(x > nrow(D)-100)),.SDcols = paste0("genotype", c("nfre", "chit5", "nfr5"))]
D2 <- data.table(OTUid = rownames(logFCranks), logFCranks > nrow(D)-100 , logFCsigns,
                 logFC = neg.bin.reg$coefficients[,-1])
D2 <- D2[OTUid %in% OTU_keep]
colnames(D2)[2:7] <- c("nfre_rank", "chit5_rank", "nfr5_rank", "nfre", "chit5", "nfr5")
D2[nfre_rank == FALSE, nfre:=0]
D2[chit5_rank == FALSE, chit5:=0]
D2[nfr5_rank == FALSE, nfr5:=0]

D2 <- D2[,c(1,5:7)]

fwrite(D2, paste0("Results/top100FC_", comp, ".csv"))

