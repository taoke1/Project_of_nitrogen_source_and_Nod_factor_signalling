# Code by Ib Thorsgaard Jensen
library(data.table)
library(limma)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Changes working directory to the path of this file. Only works in Rstudio

filtering_otus <- function(otu_table, otus){
  otu_table2 <- otu_table
  otu_omit <- otu_table2[,lapply(.SD, function(x) mean(x == 0) >= 0.9), .SDcols = otus]
  otu_omit <- colnames(otu_omit)[unlist(otu_omit)]
  otu_table2[,c(otu_omit):=NULL]
  return(otu_table2)
}

filter_group_effect <- function(otu_table, group){
  otu_table_filtered_log <- removeBatchEffect(log(otu_table+1), batch = group)
  otu_table_filtered <- apply(otu_table_filtered_log, 2,
                                  function(x) sapply(x, function(y) max( exp(y)-1, 0 ))
                                  )
  return(otu_table_filtered)
}

# Importing data, initial subsetting and creation of variables ----


# Loading datasets
OTU_table <- fread("../Data/zotu_table_0.01_clean.csv")
meta_data <- fread("../Data/design_clean.csv")

# Subsetting data
samples <- meta_data[compartments != "pooled_nodules", SampleID]
meta_data <- meta_data[SampleID %in% samples]
data_subset <- OTU_table[,c("OTUid", samples), with = FALSE]

Microbiome_nitrate <- data.table(meta_data[,c("Nitrate_supplied", "genotype", "compartments", "BioRepID")],
                                 transpose(data_subset, make.names = "OTUid", keep.names = "SampleID"))
otus <- colnames(Microbiome_nitrate)[-(1:5)]

Microbiome_nitrate[,State:="Inorganic"]
Microbiome_nitrate[genotype %in% c("Gifu", "nfre") & Nitrate_supplied == "no", State:="Symbiotic"]
Microbiome_nitrate[genotype %in% c("chit5", "nfr5") & Nitrate_supplied == "no", State:="Starved"]

# Filtering OTU table ----
data_subset_rhiz <- Microbiome_nitrate[compartments == "rhizosphere"]
data_subset_root <- Microbiome_nitrate[compartments == "root"]

# Additional subsetting ----
# rhiz
rhiz_inog <- data_subset_rhiz[State == "Inorganic"]
rhiz_inog <- filtering_otus(rhiz_inog, otus = otus)
rhiz_sym <- data_subset_rhiz[State == "Symbiotic"]
rhiz_sym <- filtering_otus(rhiz_sym, otus = otus)
rhiz_starv <- data_subset_rhiz[State == "Starved"]
rhiz_starv <- filtering_otus(rhiz_starv, otus = otus)
# root
root_inog <- data_subset_root[State == "Inorganic"]
root_inog <- filtering_otus(root_inog, otus = otus)
root_sym <- data_subset_root[State == "Symbiotic"]
root_sym <- filtering_otus(root_sym, otus = otus)
root_starv <- data_subset_root[State == "Starved"]
root_starv <- filtering_otus(root_starv, otus = otus)

## Filtering out genotype effects ----
# rhiz-inog
otu_names_inog <- colnames(rhiz_inog)[-c(1:4,ncol(rhiz_inog))]
otu_table_rhiz_inog <- transpose(rhiz_inog[,..otu_names_inog], make.names = "SampleID", keep.names = "OTUid")
otu_table_rhiz_inog_filtered <- filter_group_effect(otu_table_rhiz_inog[,-1], rhiz_inog$genotype)
otu_table_rhiz_inog_filtered <- data.table(OTUid = otu_table_rhiz_inog$OTUid, otu_table_rhiz_inog_filtered)

# rhiz-sym
otu_names_sym <- colnames(rhiz_sym)[-c(1:4,ncol(rhiz_sym))]
otu_table_rhiz_sym <- transpose(rhiz_sym[,..otu_names_sym], make.names = "SampleID", keep.names = "OTUid")
otu_table_rhiz_sym_filtered <- filter_group_effect(otu_table_rhiz_sym[,-1], rhiz_sym$genotype)
otu_table_rhiz_sym_filtered <- data.table(OTUid = otu_table_rhiz_sym$OTUid, otu_table_rhiz_sym_filtered)

# rhiz-starv
otu_names_starv <- colnames(rhiz_starv)[-c(1:4,ncol(rhiz_starv))]
otu_table_rhiz_starv <- transpose(rhiz_starv[,..otu_names_starv], make.names = "SampleID", keep.names = "OTUid")
otu_table_rhiz_starv_filtered <- filter_group_effect(otu_table_rhiz_starv[,-1], rhiz_starv$genotype)
otu_table_rhiz_starv_filtered <- data.table(OTUid = otu_table_rhiz_starv$OTUid, otu_table_rhiz_starv_filtered)

# root-inog
otu_names_inog <- colnames(root_inog)[-c(1:4,ncol(root_inog))]
otu_table_root_inog <- transpose(root_inog[,..otu_names_inog], make.names = "SampleID", keep.names = "OTUid")
otu_table_root_inog_filtered <- filter_group_effect(otu_table_root_inog[,-1], root_inog$genotype)
otu_table_root_inog_filtered <- data.table(OTUid = otu_table_root_inog$OTUid, otu_table_root_inog_filtered)

# root-sym
otu_names_sym <- colnames(root_sym)[-c(1:4,ncol(root_sym))]
otu_table_root_sym <- transpose(root_sym[,..otu_names_sym], make.names = "SampleID", keep.names = "OTUid")
otu_table_root_sym_filtered <- filter_group_effect(otu_table_root_sym[,-1], root_sym$genotype)
otu_table_root_sym_filtered <- data.table(OTUid = otu_table_root_sym$OTUid, otu_table_root_sym_filtered)

# root-starv
otu_names_starv <- colnames(root_starv)[-c(1:4,ncol(root_starv))]
otu_table_root_starv <- transpose(root_starv[,..otu_names_starv], make.names = "SampleID", keep.names = "OTUid")
otu_table_root_starv_filtered <- filter_group_effect(otu_table_root_starv[,-1], root_starv$genotype)
otu_table_root_starv_filtered <- data.table(OTUid = otu_table_root_starv$OTUid, otu_table_root_starv_filtered)


# Saving filtered tables ----
data_list <- list(rhiz_inog = otu_table_rhiz_inog_filtered,
                  rhiz_sym = otu_table_rhiz_sym_filtered,
                  rhiz_starv = otu_table_rhiz_starv_filtered,
                  root_inog = otu_table_root_inog_filtered,
                  root_sym = otu_table_root_sym_filtered,
                  root_starv = otu_table_root_starv_filtered)
for(i in 1:6){
  A <- data_list[[i]]
  colnames(A)[1] <- "#OTU ID"
  fwrite(A, paste0("Filtered_tables/", names(data_list)[i], ".tsv"), sep = "\t")
}



