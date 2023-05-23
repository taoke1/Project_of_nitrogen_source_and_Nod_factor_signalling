# Code by Ib Thorsgaard Jensen
library(data.table)
library(magrittr)
library(edgeR)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Changes working directory to the path of this file. Only works in Rstudio

# Loading datasets
OTU_table <- fread("../Data/zotu_table_0.01_clean.csv")
meta_data <- fread("../Data/design_clean.csv")
taxonomy <- fread(file="../Data/Taxonomy.txt")
taxonomy[G_identity<0.7, Genus:="Unknown"]

# Choose compartment (root or rhizosphere)
comp <- "rhizosphere"

# subsetting meta data
samples <- meta_data[compartments == comp, SampleID]
meta_data <- meta_data[SampleID %in% samples]
data_subset <- OTU_table[,c("OTUid", samples), with = FALSE]

# Adding the nutritional state as a variable to the meta data
meta_data[,State:=fcase(
  Nitrate_supplied == "yes", "Inorganic",
  genotype %in% c("Gifu", "nfre") & Nitrate_supplied == "no", "Symbiotic",
  genotype %in% c("chit5", "nfr5") & Nitrate_supplied == "no", "Starved"
)]

# Collapsing OTUs to the genus level
taxonomy <- taxonomy[match(data_subset$OTUid, OTUid)]
feature.table <- rowsum(data_subset[,-1], taxonomy$Genus)

# Setting the significance level and effect-size cutoff
effect.cutoff <- 1
sig.lvl <- 0.05

# Setting up DGE object and fitting negative binomial model with edgeR
data.DGE <- DGEList(counts=feature.table, group=meta_data$SampleID)
data.DGE <- calcNormFactors(data.DGE)

design.matrix <- model.matrix(~ 0 + State, data = meta_data)
data.DGE.disp <- estimateDisp(data.DGE, design.matrix, robust = T)
neg.bin.reg <- glmQLFit(data.DGE.disp, design.matrix)

# Testing for differential abundance between states.
enrich_inog_sym_test <- glmTreat(neg.bin.reg, contrast = c(1, 0, -1), lfc = effect.cutoff) %>%
  decideTestsDGE(adjust.method = "fdr", p.value = sig.lvl)
enrich_inog_starv_test <- glmTreat(neg.bin.reg, contrast = c(1, -1, 0), lfc = effect.cutoff) %>%
  decideTestsDGE(adjust.method = "fdr", p.value = sig.lvl)
enrich_sym_starv_test <- glmTreat(neg.bin.reg, contrast = c(0, -1, 1), lfc = effect.cutoff) %>%
  decideTestsDGE(adjust.method = "fdr", p.value = sig.lvl)

# Genera enriched in the inorganic state compared to the other two
enrich_inog_sym <- rownames(feature.table)[enrich_inog_sym_test == 1]
enrich_inog_starv <- rownames(feature.table)[enrich_inog_starv_test == 1]
enrich_inog <- intersect(enrich_inog_sym, enrich_inog_starv)
enrich_inog <- enrich_inog[enrich_inog != "Unknown"]

# Genera enriched in the symbiotic state compared to the other two
enrich_sym_inog <- rownames(feature.table)[enrich_inog_sym_test == -1]
enrich_sym_starv <- rownames(feature.table)[enrich_sym_starv_test == 1]
enrich_sym <- intersect(enrich_sym_inog, enrich_sym_starv)
enrich_sym <- enrich_sym[enrich_sym != "Unknown"]

# Genera enriched in the starved state compared to the other two
enrich_starv_inog <- rownames(feature.table)[enrich_inog_starv_test == -1]
enrich_starv_sym <- rownames(feature.table)[enrich_sym_starv_test == -1]
enrich_starv <- intersect(enrich_starv_inog, enrich_starv_sym)
enrich_starv <- enrich_starv[enrich_starv != "Unknown"]

# Setting up vector denoting which genera are enriched in which state, if any
enrich <- rep("None", nrow(feature.table)); names(enrich) <- rownames(feature.table)
enrich[names(enrich) %in% enrich_inog] <- "Inorganic"
enrich[names(enrich) %in% enrich_sym] <- "Symbiotic"
enrich[names(enrich) %in% enrich_starv] <- "Starved"

# Genera depleted in the inorganic state compared to the other two
deplete_inog_sym <- rownames(feature.table)[enrich_inog_sym_test == -1]
deplete_inog_starv <- rownames(feature.table)[enrich_inog_starv_test == -1]
deplete_inog <- intersect(deplete_inog_sym, deplete_inog_starv)
deplete_inog <- deplete_inog[deplete_inog != "Unknown"]

# Genera depleted in the symbiotic state compared to the other two
deplete_sym_inog <- rownames(feature.table)[enrich_inog_sym_test == 1]
deplete_sym_starv <- rownames(feature.table)[enrich_sym_starv_test == -1]
deplete_sym <- intersect(deplete_sym_inog, deplete_sym_starv)
deplete_sym <- deplete_sym[deplete_sym != "Unknown"]

# Genera depleted in the starved state compared to the other two
deplete_starv_inog <- rownames(feature.table)[enrich_inog_starv_test == 1]
deplete_starv_sym <- rownames(feature.table)[enrich_sym_starv_test == 1]
deplete_starv <- intersect(deplete_starv_inog, deplete_starv_sym)
deplete_starv <- deplete_starv[deplete_starv != "Unknown"]

# Setting up vector denoting which genera are depleted in which state, if any
deplete <- rep("None", nrow(feature.table)); names(deplete) <- rownames(feature.table)
deplete[names(deplete) %in% deplete_inog] <- "Inorganic"
deplete[names(deplete) %in% deplete_sym] <- "Symbiotic"
deplete[names(deplete) %in% deplete_starv] <- "Starved"

DA_table <- data.table(Genus = rownames(feature.table), Enriched = enrich, Depleted = deplete)
DA_table <- DA_table[Enriched != "None" | Depleted != "None"]

# Saving table with enriched OTUs for later inspection
filename <- c("State_Associated_genera", comp) %>% paste(collapse = "_") %>% paste("csv", sep = ".")
dirname <- c("Results", "Nutrtion_associated_bacteria", filename) %>% paste(collapse = "/")
fwrite(DA_table, file = dirname)
