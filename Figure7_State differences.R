# Code by Ib Thorsgaard Jensen
library(data.table)
library(ggplot2)
library(ggtern)
library(lmerTest)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Changes working directory to the path of this file. Only works in Rstudio

# Loading data
metabolite_data <- fread("../Data/Metabolites_filtered_lmer.csv")
metabolite_data <- metabolite_data[Sample_type == "Gifu"]
meta_data <- metabolite_data[,1:3]
anno_file <- fread("../Data/canopus_compound_summary.tsv")

# Changes column names of annotation file
anno_file[,Feature_ID:=gsub(".*_","",Feature_ID)]
anno_file[,Feature_ID:=paste0("Feature", Feature_ID)]
colnames(anno_file)[1] <- "FeatureID"

meta_data[,State:=factor(State)]
states <- c("Starved", "Symbiotic", "Inorganic", "Inorganic+R7A")

# Feature level--------------------------------------------------------------------------------------------------------------------
# Testing the difference between states with a linear mixed model
A <- array(dim = c(4,4,ncol(metabolite_data)-3))
dimnames(A)[[1]] <- states
dimnames(A)[[2]] <- states
set.seed(1676559840)
for(i in 1:(ncol(metabolite_data)-3)){
  cat(i, "\r")
  # Metabolite variable
  peak <- metabolite_data[[i+3]]
  # Replacing zero-values with a small value prior to log transformation
  peak[peak == 0] <- 0.5
  peak <- log(peak)
  # Linear mixed models. A model is fitted with each state as a reference level to get all pairwise comparisons
  l_inog <- lmer(peak ~ State + (1|Sample_ID), data = meta_data)
  meta_data[,State:=relevel(State, ref = "Symbiotic")]
  l_sym <- lmer(peak ~ State + (1|Sample_ID), data = meta_data)
  meta_data[,State:=relevel(State, ref = "Starved")]
  l_starv <- lmer(peak ~ State + (1|Sample_ID), data = meta_data)
  meta_data[,State:=relevel(State, ref = "Inorganic+R7A")]
  l_inogR7A <- lmer(peak ~ State + (1|Sample_ID), data = meta_data)
  meta_data[,State:=relevel(State, ref = "Inorganic")]
  
  # Sets up matrix with p-values for all pairwise comparisons
  states_inog <- setdiff(states, "Inorganic")
  states_inogR7A <- setdiff(states, "Inorganic+R7A")
  states_sym <- setdiff(states, "Symbiotic")
  states_starv <- setdiff(states, "Starved")
  C1 <- coef(summary(l_inog)); rownames(C1) <- gsub("State", "", rownames(C1))
  C2 <- coef(summary(l_inogR7A)); rownames(C2) <- gsub("State", "", rownames(C2))
  C3 <- coef(summary(l_sym)); rownames(C3) <- gsub("State", "", rownames(C3))
  C4 <- coef(summary(l_starv)); rownames(C4) <- gsub("State", "", rownames(C4))
  
  # Inserting results from above into array
  A[,,i][states_inog, "Inorganic"] <- C1[states_inog, 5]
  A[,,i][states_inogR7A, "Inorganic+R7A"] <- C2[states_inogR7A, 5]
  A[,,i][states_sym, "Symbiotic"] <- C3[states_sym, 5]
  A[,,i][states_starv, "Starved"] <- C4[states_starv, 5]
  diag(A[,,i]) <- 1
}

# Correction for multiple testing
B <- array(0, c(4,4,ncol(metabolite_data)-3))
for(i in 1:4){
  for(j in 1:4){
    B[i,j,] <- p.adjust(A[i,j,], method = "fdr")
  }
}

dimnames(B)[[1]] <- states
dimnames(B)[[2]] <- states

# Table with results
res <- data.table(Feature = colnames(metabolite_data)[-(1:3)],
                  `Starved vs symbiotic` = B["Starved", "Symbiotic",]<0.05,
                  `Starved vs inorganic` = B["Starved", "Inorganic",]<0.05,
                  `Starved vs inorganic+R7A` = B["Starved", "Inorganic+R7A",]<0.05,
                  `Symbiotic vs inorganic` = B["Symbiotic", "Inorganic",]<0.05,
                  `Symbiotic vs inorganic+R7A` = B["Symbiotic", "Inorganic+R7A",]<0.05,
                  `Inorganic vs inorganic+R7A` = B["Inorganic", "Inorganic+R7A",]<0.05)
# Number of significantly different features for each comparison
res[,lapply(.SD, sum), .SDcols = colnames(res)[-1]]

# Setting up tables showing which metabolites are enriched in which conditions
metabolite_data[,lapply(.SD, mean),
                .SDcols = colnames(metabolite_data)[-(1:3)],
                State] -> Means
Means <- transpose(Means, make.names = "State", keep.names = "Feature")
Means[,":="(
  `Starved vs symbiotic` = ifelse(Starved > Symbiotic, 1, -1),
  `Starved vs inorganic` = ifelse(Starved > Inorganic, 1, -1),
  `Starved vs inorganic+R7A` = ifelse(Starved > `Inorganic+R7A`, 1, -1),
  `Symbiotic vs inorganic` = ifelse(Symbiotic > Inorganic, 1, -1),
  `Symbiotic vs inorganic+R7A` = ifelse(Symbiotic > `Inorganic+R7A`, 1, -1),
  `Inorganic vs inorganic+R7A` = ifelse(Inorganic > `Inorganic+R7A`, 1, -1)
)]

# Constructs matrix with 1, 0, -1 according to if and where a feature is enriched.
# A feature is enriched in state x when "x vs y" is 1 for that feature. If "x vs y" is
# -1, then the feature is enriched in y.
C <- as.matrix(Means[,6:11])
C2 <- as.matrix(res[,-1])
C[!C2] <- 0
Result <- data.table(Feature = colnames(metabolite_data)[-(1:3)], C)

fwrite(Result, "Results/Feature_differences.csv")

## Feature level ternary plot for symbiotic-starved-inorganic ----

# Identifying features enriched in each feature compared to the other two (excluding Inorganic+R7A)
Starv_enrich <- intersect(Result[`Starved vs symbiotic` == 1, Feature],
                          Result[`Starved vs inorganic` == 1, Feature])
Sym_enrich <- intersect(Result[`Starved vs symbiotic` == -1, Feature],
                        Result[`Symbiotic vs inorganic` == 1, Feature])
Inog_enrich <- intersect(Result[`Symbiotic vs inorganic` == -1, Feature],
                         Result[`Starved vs inorganic` == -1, Feature])

# Table of means by state
metabolite_data_log <- data.table(State = metabolite_data$State, metabolite_data[,-(1:3)])
Mean_met <- metabolite_data_log[,lapply(.SD, mean), by = State, .SDcols = colnames(metabolite_data_log)[-1]]
Mean_met <- transpose(Mean_met, make.names = "State", keep.names = "Feature")

# Vector indicating which of the states (if any) each feature is enriched in
Enrich <- fcase(
  !(Mean_met$Feature %in% c(Inog_enrich, Sym_enrich, Starv_enrich)), "None",
  Mean_met$Feature %in% Inog_enrich, "Inorganic",
  Mean_met$Feature %in% Sym_enrich, "Symbiotic",
  Mean_met$Feature %in% Starv_enrich, "Starved"
)

# Vector of overall means by feature
all_means <- metabolite_data_log[State != "Inorganic+R7A", lapply(.SD, function(x) mean(log(x+1))), .SDcols = colnames(metabolite_data_log)[-1]]
all_means <- transpose(all_means, keep.names = "Feature")

# Collects data necassary for ternary plot in a data table to pass to ggtern
plot_data <- data.table(Mean_met, Overall = all_means$V1, Enrich = Enrich)
plot_data[,Enrich:=factor(Enrich, levels = c("None", "Inorganic", "Starved", "Symbiotic"))]
plot_data <- plot_data[order(Enrich)]
col.table <- data.table(gt = c("Inorganic", "Starved", "Symbiotic", "Inorganic+R7A", "None"),
                        color = c("#CC79A7", "#0072B2", "#009E73", "#D95F02", "grey") )

# Ternary plot
ggtern(data = plot_data, aes(x = Symbiotic, y = Inorganic, z = Starved, size = Overall, col = Enrich))+
  geom_point(alpha = 0.5)+
  scale_color_manual(name = "Enriched Feature",
                     breaks = col.table$gt,
                     values = col.table$color)+
  scale_size_continuous(name = "Overall average")+
  theme(text = element_text(size = 20),
        legend.text=element_text(size=20))+
  theme(tern.axis.title.L = element_text(hjust = 0),
        tern.axis.title.R = element_text(hjust = 1)) -> tern_plot

# Saving ternary plot
ggsave(filename = paste0("Figures/Ternary/Ternary_plot_Feature_", as.character(Sys.Date()),".pdf"),
       tern_plot, width = 300, height = 180, units = "mm")

# Saving results shown on ternary plots
fwrite(plot_data, file = paste0("Results/Ternary_plot_Feature_", as.character(Sys.Date()), "_filter.csv"))

# Filtering for higher annotation levels---------------------------------------------------------------------------------------
# Removes features from the annotation files that are not in the filtered feature table
anno_file <- anno_file[FeatureID %in% colnames(metabolite_data)[-(1:3)]]
# Setting as unknown the annotations with low probability
anno_file[`NPC#pathway Probability`<0.6, `NPC#pathway`:="unknown"]
anno_file[`NPC#class Probability`<0.6, `NPC#class` := "unknown"]
anno_file[`ClassyFire#most specific class Probability`<0.6, `ClassyFire#most specific class`:="unknown"]
anno_file[`ClassyFire#superclass probability`<0.6, `ClassyFire#superclass`:="unknown"]
# Removes features from the feature table that are not annotated
filtered_features <- anno_file[FeatureID %in% colnames(metabolite_data), FeatureID]
non_annotated <- colnames(metabolite_data)[-(1:3)][ !(colnames(metabolite_data)[-(1:3)] %in% filtered_features) ]
metabolite_data[,c(non_annotated):=NULL]
anno_file <- anno_file[match(colnames(metabolite_data)[-(1:3)], FeatureID)]
# Sanity check - do the feature table and annotation file contain the same features, and
#  are they arranged in the same order?
all(colnames(metabolite_data)[-(1:3)] == anno_file$FeatureID)

# Most specific class ------------------------------------------------------------------------------------------------------------
# Setting up feature table collapsed by most specific class
M <- metabolite_data[,-(1:3)]
metabolite_trans <- transpose(M, keep.names = "FeatureID")
metabolite_MSC <- merge(anno_file[,c("FeatureID", "ClassyFire#most specific class")], metabolite_trans, by = "FeatureID")
metabolite_MSC <- metabolite_MSC[,-1]
colnames(metabolite_MSC)[1] <- "Most specific class"
metabolite_MSC <- metabolite_MSC[,lapply(.SD, sum), `Most specific class`]
metabolite_MSC <- metabolite_MSC[`Most specific class` != "unknown"]
colnames(metabolite_MSC)[-1] <- metabolite_data$Sample_ID
metabolite_MSC <- transpose(metabolite_MSC, make.names = "Most specific class", keep.names = "Sample_ID")
metabolite_MSC[,State:=factor(metabolite_data$State)]
setcolorder(metabolite_MSC, c("Sample_ID", "State"))

# Testing the difference between states with a linear mixed model
A <- array(dim = c(4,4,ncol(metabolite_MSC)-2))
dimnames(A)[[1]] <- states
dimnames(A)[[2]] <- states
set.seed(1676559840)
for(i in 1:(ncol(metabolite_MSC)-2)){
  cat(i, "\r")
  # Metabolite variable
  peak <- metabolite_MSC[[i+2]]
  # Replacing zero-values with a small value prior to log transformation
  peak[peak == 0] <- 0.5
  peak <- log(peak)
  # Linear mixed models. A model is fitted with each state as a reference level to get all pairwise comparisons
  l_inog <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_MSC)
  metabolite_MSC[,State:=relevel(State, ref = "Symbiotic")]
  l_sym <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_MSC)
  metabolite_MSC[,State:=relevel(State, ref = "Starved")]
  l_starv <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_MSC)
  metabolite_MSC[,State:=relevel(State, ref = "Inorganic+R7A")]
  l_inogR7A <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_MSC)
  metabolite_MSC[,State:=relevel(State, ref = "Inorganic")]
  
  # Sets up matrix with p-values for all pairwise comparisons
  states_inog <- setdiff(states, "Inorganic")
  states_inogR7A <- setdiff(states, "Inorganic+R7A")
  states_sym <- setdiff(states, "Symbiotic")
  states_starv <- setdiff(states, "Starved")
  C1 <- coef(summary(l_inog)); rownames(C1) <- gsub("State", "", rownames(C1))
  C2 <- coef(summary(l_inogR7A)); rownames(C2) <- gsub("State", "", rownames(C2))
  C3 <- coef(summary(l_sym)); rownames(C3) <- gsub("State", "", rownames(C3))
  C4 <- coef(summary(l_starv)); rownames(C4) <- gsub("State", "", rownames(C4))
  
  # Inserting results from above into array
  A[,,i][states_inog, "Inorganic"] <- C1[states_inog, 5]
  A[,,i][states_inogR7A, "Inorganic+R7A"] <- C2[states_inogR7A, 5]
  A[,,i][states_sym, "Symbiotic"] <- C3[states_sym, 5]
  A[,,i][states_starv, "Starved"] <- C4[states_starv, 5]
  diag(A[,,i]) <- 1
}
# Correction for multiple testing
B <- array(0, c(4,4,ncol(metabolite_MSC)-2))
for(i in 1:4){
  for(j in 1:4){
    B[i,j,] <- p.adjust(A[i,j,], method = "fdr")
  }
}

dimnames(B)[[1]] <- states
dimnames(B)[[2]] <- states

# Table with results
MSC_names <- colnames(metabolite_MSC)[-(1:2)]
res <- data.table(`Most specific class` = MSC_names,
                  `Starv vs symbiotic` = B["Starved", "Symbiotic",]<0.05,
                  `Starv vs inorganic` = B["Starved", "Inorganic",]<0.05,
                  `Starv vs inorganic+R7A` = B["Starved", "Inorganic+R7A",]<0.05,
                  `Symbiotic vs inorganic` = B["Symbiotic", "Inorganic",]<0.05,
                  `Symbiotic vs inorganic+R7A` = B["Symbiotic", "Inorganic+R7A",]<0.05,
                  `Inorganic vs inorganic+R7A` = B["Inorganic", "Inorganic+R7A",]<0.05)

# Number of significantly different features for each comparison
res[,lapply(.SD, sum), .SDcols = colnames(res)[-1]]

# Setting up tables showing which metabolites are enriched in which conditions
metabolite_MSC[,lapply(.SD, mean),
               .SDcols = MSC_names,
               State] -> Means
Means <- transpose(Means, make.names = "State", keep.names = "Most specific class")
Means[,":="(
  `Starved vs symbiotic` = ifelse(Starved > Symbiotic, 1, -1),
  `Starved vs inorganic` = ifelse(Starved > Inorganic, 1, -1),
  `Starved vs inorganic+R7A` = ifelse(Starved > `Inorganic+R7A`, 1, -1),
  `Symbiotic vs inorganic` = ifelse(Symbiotic > Inorganic, 1, -1),
  `Symbiotic vs inorganic+R7A` = ifelse(Symbiotic > `Inorganic+R7A`, 1, -1),
  `Inorganic vs inorganic+R7A` = ifelse(Inorganic > `Inorganic+R7A`, 1, -1)
)]

# Constructs matrix with 1, 0, -1 according to if and where a feature is enriched.
# A feature is enriched in state x when "x vs y" is 1 for that feature. If "x vs y" is
# -1, then the feature is enriched in y.
C <- as.matrix(Means[,6:11])
C2 <- as.matrix(res[,-1])
C[!C2] <- 0
Result <- data.table(`Most specific class` = MSC_names, C)

fwrite(Result, "Results/Most_specfic_class_differences.csv")

## Most specific class ternary plots - Inorganic ----

# Identifying features enriched in each feature compared to the other two (excluding Inorganic+R7A)
Starv_enrich <- intersect(Result[`Starved vs symbiotic` == 1, `Most specific class`],
                          Result[`Starved vs inorganic` == 1, `Most specific class`])
Sym_enrich <- intersect(Result[`Starved vs symbiotic` == -1, `Most specific class`],
                        Result[`Symbiotic vs inorganic` == 1, `Most specific class`])
Inog_enrich <- intersect(Result[`Symbiotic vs inorganic` == -1, `Most specific class`],
                         Result[`Starved vs inorganic` == -1, `Most specific class`])

# Table of means by state
metabolite_MSC_log <- data.table(State = metabolite_MSC$State, metabolite_MSC[,..MSC_names])
Mean_met <- metabolite_MSC_log[,lapply(.SD, mean), by = State, .SDcols = colnames(metabolite_MSC_log)[-1]]
Mean_met <- transpose(Mean_met, make.names = "State", keep.names = "Most specific class")

# Vector indicating which of the states (if any) each feature is enriched in
Enrich <- fcase(
  !(Mean_met$`Most specific class` %in% c(Inog_enrich, Sym_enrich, Starv_enrich)), "None",
  Mean_met$`Most specific class` %in% Inog_enrich, "Inorganic",
  Mean_met$`Most specific class` %in% Sym_enrich, "Symbiotic",
  Mean_met$`Most specific class` %in% Starv_enrich, "Starved"
)

# Vector of overall means by feature
all_means <- metabolite_MSC_log[State != "Inorganic+R7A", lapply(.SD, function(x) mean(log(x+1))), .SDcols = colnames(metabolite_MSC_log)[-1]]
all_means <- transpose(all_means, keep.names = "Feature")

# Collects data necassary for ternary plot in a data table to pass to ggtern
plot_data <- data.table(Mean_met, Overall = all_means$V1, Enrich = Enrich)
plot_data[,Enrich:=factor(Enrich, levels = c("None", "Inorganic", "Starved", "Symbiotic"))]
plot_data <- plot_data[order(Enrich)]
col.table <- data.table(gt = c("Inorganic", "Starved", "Symbiotic", "Inorganic+R7A", "None"),
                        color = c("#CC79A7", "#0072B2", "#009E73", "#D95F02", "grey") )

# Ternary plot
ggtern(data = plot_data, aes(x = Symbiotic, y = Inorganic, z = Starved, size = Overall, col = Enrich))+
  geom_point(alpha = 0.5)+
  scale_color_manual(name = "Enriched Feature",
                     breaks = col.table$gt,
                     values = col.table$color)+
  scale_size_continuous(name = "Overall average")+
  theme(text = element_text(size = 20),
        legend.text=element_text(size=20))+
  theme(tern.axis.title.L = element_text(hjust = 0),
        tern.axis.title.R = element_text(hjust = 1)) -> tern_plot

# Saving ternary plot
ggsave(filename = paste0("Figures/Ternary/Ternary_MSC_", as.character(Sys.Date()),".pdf"),
       tern_plot, width = 300, height = 180, units = "mm")

# Saving results shown on ternary plots
fwrite(plot_data, file = paste0("Results/Ternary_plot_table_", as.character(Sys.Date()), "_filter.csv"))

# Superclass --------------------------------------------------------------------------------------------------------------------
# Setting up feature table collapsed by most superclass
M <- metabolite_data[,-(1:3)]
metabolite_trans <- transpose(M, keep.names = "FeatureID")
metabolite_superclass <- merge(anno_file[,c("FeatureID", "ClassyFire#superclass")], metabolite_trans, by = "FeatureID")
metabolite_superclass <- metabolite_superclass[,-1]
colnames(metabolite_superclass)[1] <- "Superclass"
metabolite_superclass <- metabolite_superclass[,lapply(.SD, sum), Superclass]
metabolite_superclass <- metabolite_superclass[Superclass != "unknown"]
colnames(metabolite_superclass)[-1] <- metabolite_data$Sample_ID
metabolite_superclass <- transpose(metabolite_superclass, make.names = "Superclass", keep.names = "Sample_ID")
metabolite_superclass[,State:=metabolite_data$State]
metabolite_superclass[,State:=metabolite_data$State]
setcolorder(metabolite_superclass, c("Sample_ID", "State"))
metabolite_superclass[,State:=as.factor(State)]

# Testing the difference between states with a linear mixed model
A <- array(dim = c(4,4,ncol(metabolite_superclass)-2))
dimnames(A)[[1]] <- states
dimnames(A)[[2]] <- states
set.seed(1676559840)
for(i in 1:(ncol(metabolite_superclass)-2)){
  cat(i, "\r")
  # Metabolite variable
  peak <- metabolite_superclass[[i+2]]
  # Replacing zero-values with a small value prior to log transformation
  peak[peak == 0] <- 0.5
  peak <- log(peak)
  # Linear mixed models. A model is fitted with each state as a reference level to get all pairwise comparisons
  l_inog <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_superclass)
  metabolite_superclass[,State:=relevel(State, ref = "Symbiotic")]
  l_sym <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_superclass)
  metabolite_superclass[,State:=relevel(State, ref = "Starved")]
  l_starv <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_superclass)
  metabolite_superclass[,State:=relevel(State, ref = "Inorganic+R7A")]
  l_inogR7A <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_superclass)
  metabolite_superclass[,State:=relevel(State, ref = "Inorganic")]
  
  # Sets up matrix with p-values for all pairwise comparisons
  states_inog <- setdiff(states, "Inorganic")
  states_inogR7A <- setdiff(states, "Inorganic+R7A")
  states_sym <- setdiff(states, "Symbiotic")
  states_starv <- setdiff(states, "Starved")
  C1 <- coef(summary(l_inog)); rownames(C1) <- gsub("State", "", rownames(C1))
  C2 <- coef(summary(l_inogR7A)); rownames(C2) <- gsub("State", "", rownames(C2))
  C3 <- coef(summary(l_sym)); rownames(C3) <- gsub("State", "", rownames(C3))
  C4 <- coef(summary(l_starv)); rownames(C4) <- gsub("State", "", rownames(C4))
  
  # Inserting results from above into array
  A[,,i][states_inog, "Inorganic"] <- C1[states_inog, 5]
  A[,,i][states_inogR7A, "Inorganic+R7A"] <- C2[states_inogR7A, 5]
  A[,,i][states_sym, "Symbiotic"] <- C3[states_sym, 5]
  A[,,i][states_starv, "Starved"] <- C4[states_starv, 5]
  diag(A[,,i]) <- 1
}

# Correction for multiple testing
B <- array(0, c(4,4,ncol(metabolite_superclass)-2))
for(i in 1:4){
  for(j in 1:4){
    B[i,j,] <- p.adjust(A[i,j,], method = "fdr")
  }
}

dimnames(B)[[1]] <- states
dimnames(B)[[2]] <- states

# Table with results
res <- data.table(Superclass = colnames(metabolite_superclass)[-(1:2)],
                  `Starved vs symbiotic` = B["Starved", "Symbiotic",]<0.05,
                  `Starved vs inorganic` = B["Starved", "Inorganic",]<0.05,
                  `Starved vs inorganic+R7A` = B["Starved", "Inorganic+R7A",]<0.05,
                  `Symbiotic vs inorganic` = B["Symbiotic", "Inorganic",]<0.05,
                  `Symbiotic vs inorganic+R7A` = B["Symbiotic", "Inorganic+R7A",]<0.05,
                  `Inorganic vs inorganic+R7A` = B["Inorganic", "Inorganic+R7A",]<0.05)

# Number of significantly different features for each comparison
res[,lapply(.SD, sum), .SDcols = colnames(res)[-1]]

# Setting up tables showing which metabolites are enriched in which conditions
metabolite_superclass[,lapply(.SD, mean),
                      .SDcols = colnames(metabolite_superclass)[-(1:2)],
                      State] -> Means
Means <- transpose(Means, make.names = "State", keep.names = "Superclass")
Means[,":="(
  `Starved vs symbiotic` = ifelse(Starved > Symbiotic, 1, -1),
  `Starved vs inorganic` = ifelse(Starved > Inorganic, 1, -1),
  `Starved vs inorganic+R7A` = ifelse(Starved > `Inorganic+R7A`, 1, -1),
  `Symbiotic vs inorganic` = ifelse(Symbiotic > Inorganic, 1, -1),
  `Symbiotic vs inorganic+R7A` = ifelse(Symbiotic > `Inorganic+R7A`, 1, -1),
  `Inorganic vs inorganic+R7A` = ifelse(Inorganic > `Inorganic+R7A`, 1, -1)
)]

# Constructs matrix with 1, 0, -1 according to if and where a feature is enriched.
# A feature is enriched in state x when "x vs y" is 1 for that feature. If "x vs y" is
# -1, then the feature is enriched in y.
C <- as.matrix(Means[,6:11])
C2 <- as.matrix(res[,-1])
C[!C2] <- 0
Result <- data.table(Superclass = colnames(metabolite_superclass)[-(1:2)], C)

fwrite(Result, "Results/Superclass_differences.csv")

# NPC class ----------------------------------------------------------------------------------------------------------------------
# Setting up feature table collapsed by class
M <- metabolite_data[,-(1:3)]
metabolite_trans <- transpose(M, keep.names = "FeatureID")
metabolite_class <- merge(anno_file[,c("FeatureID", "NPC#class")], metabolite_trans, by = "FeatureID")
metabolite_class <- metabolite_class[,-1]
colnames(metabolite_class)[1] <- "Class"
metabolite_class <- metabolite_class[,lapply(.SD, sum), Class]
metabolite_class <- metabolite_class[Class != "unknown"]
colnames(metabolite_class)[-1] <- metabolite_data$Sample_ID
metabolite_class <- transpose(metabolite_class, make.names = "Class", keep.names = "Sample_ID")
metabolite_class[,State:=metabolite_data$State]
setcolorder(metabolite_class, c("Sample_ID", "State"))
metabolite_class[,State:=as.factor(State)]

# Testing the difference between states with a linear mixed model
A <- array(dim = c(4,4,ncol(metabolite_class)-2))
dimnames(A)[[1]] <- states
dimnames(A)[[2]] <- states
set.seed(1676559840)
for(i in 1:(ncol(metabolite_class)-2)){
  cat(i, "\r")
  # Metabolite variable
  peak <- metabolite_class[[i+2]]
  # Replacing zero-values with a small value prior to log transformation
  peak[peak == 0] <- 0.5
  peak <- log(peak)
  # Linear mixed models. A model is fitted with each state as a reference level to get all pairwise comparisons
  l_inog <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_class)
  metabolite_class[,State:=relevel(State, ref = "Symbiotic")]
  l_sym <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_class)
  metabolite_class[,State:=relevel(State, ref = "Starved")]
  l_starv <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_class)
  metabolite_class[,State:=relevel(State, ref = "Inorganic+R7A")]
  l_inogR7A <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_class)
  metabolite_class[,State:=relevel(State, ref = "Inorganic")]
  
  # Sets up matrix with p-values for all pairwise comparisons
  states_inog <- setdiff(states, "Inorganic")
  states_inogR7A <- setdiff(states, "Inorganic+R7A")
  states_sym <- setdiff(states, "Symbiotic")
  states_starv <- setdiff(states, "Starved")
  C1 <- coef(summary(l_inog)); rownames(C1) <- gsub("State", "", rownames(C1))
  C2 <- coef(summary(l_inogR7A)); rownames(C2) <- gsub("State", "", rownames(C2))
  C3 <- coef(summary(l_sym)); rownames(C3) <- gsub("State", "", rownames(C3))
  C4 <- coef(summary(l_starv)); rownames(C4) <- gsub("State", "", rownames(C4))
  
  # Inserting results from above into array
  A[,,i][states_inog, "Inorganic"] <- C1[states_inog, 5]
  A[,,i][states_inogR7A, "Inorganic+R7A"] <- C2[states_inogR7A, 5]
  A[,,i][states_sym, "Symbiotic"] <- C3[states_sym, 5]
  A[,,i][states_starv, "Starved"] <- C4[states_starv, 5]
  diag(A[,,i]) <- 1
}

# Correction for multiple testing
B <- array(0, c(4,4,ncol(metabolite_class)-2))
for(i in 1:4){
  for(j in 1:4){
    B[i,j,] <- p.adjust(A[i,j,], method = "fdr")
  }
}

dimnames(B)[[1]] <- states
dimnames(B)[[2]] <- states

class_names <- colnames(metabolite_class)[-(1:2)]

# Table with results
res <- data.table(Class = class_names,
                  `Starved vs symbiotic` = B["Starved", "Symbiotic",]<0.05,
                  `Starved vs inorganic` = B["Starved", "Inorganic",]<0.05,
                  `Starved vs inorganic+R7A` = B["Starved", "Inorganic+R7A",]<0.05,
                  `Symbiotic vs inorganic` = B["Symbiotic", "Inorganic",]<0.05,
                  `Symbiotic vs inorganic+R7A` = B["Symbiotic", "Inorganic+R7A",]<0.05,
                  `Inorganic vs inorganic+R7A` = B["Inorganic", "Inorganic+R7A",]<0.05)

# Number of significantly different features for each comparison
res[,lapply(.SD, sum), .SDcols = colnames(res)[-1]]

# Setting up tables showing which metabolites are enriched in which conditions
metabolite_class[,lapply(.SD, mean),
                 .SDcols = colnames(metabolite_class)[-(1:2)],
                 State] -> Means
Means <- transpose(Means, make.names = "State", keep.names = "Class")
Means[,":="(
  `Starved vs symbiotic` = ifelse(Starved > Symbiotic, 1, -1),
  `Starved vs inorganic` = ifelse(Starved > Inorganic, 1, -1),
  `Starved vs inorganic+R7A` = ifelse(Starved > `Inorganic+R7A`, 1, -1),
  `Symbiotic vs inorganic` = ifelse(Symbiotic > Inorganic, 1, -1),
  `Symbiotic vs inorganic+R7A` = ifelse(Symbiotic > `Inorganic+R7A`, 1, -1),
  `Inorganic vs inorganic+R7A` = ifelse(Inorganic > `Inorganic+R7A`, 1, -1)
)]

# Constructs matrix with 1, 0, -1 according to if and where a feature is enriched.
# A feature is enriched in state x when "x vs y" is 1 for that feature. If "x vs y" is
# -1, then the feature is enriched in y.
C <- as.matrix(Means[,6:11])
C2 <- as.matrix(res[,-1])
C[!C2] <- 0
Result <- data.table(Class = colnames(metabolite_class)[-(1:2)], C)

fwrite(Result, "Results/Class_differences.csv")

## NPC class ternary plot ----

# Identifying features enriched in each feature compared to the other two (excluding Inorganic+R7A)
Starv_enrich <- intersect(Result[`Starved vs symbiotic` == 1, Class],
                          Result[`Starved vs inorganic` == 1, Class])
Sym_enrich <- intersect(Result[`Starved vs symbiotic` == -1, Class],
                        Result[`Symbiotic vs inorganic` == 1, Class])
Inog_enrich <- intersect(Result[`Symbiotic vs inorganic` == -1, Class],
                         Result[`Starved vs inorganic` == -1, Class])

# Table of means by state
metabolite_class_log <- data.table(State = metabolite_class$State, metabolite_class[,..class_names])
Mean_met <- metabolite_class_log[,lapply(.SD, mean), by = State, .SDcols = colnames(metabolite_class_log)[-1]]
Mean_met <- transpose(Mean_met, make.names = "State", keep.names = "Class")

# Vector indicating which of the states (if any) each feature is enriched in
Enrich <- fcase(
  !(Mean_met$Class %in% c(Inog_enrich, Sym_enrich, Starv_enrich)), "None",
  Mean_met$Class %in% Inog_enrich, "Inorganic",
  Mean_met$Class %in% Sym_enrich, "Symbiotic",
  Mean_met$Class %in% Starv_enrich, "Starved"
)

# Vector of overall means by feature
all_means <- metabolite_class_log[State != "Inorganic+R7A", lapply(.SD, function(x) mean(log(x+1))), .SDcols = colnames(metabolite_class_log)[-1]]
all_means <- transpose(all_means, keep.names = "Feature")

# Collects data necassary for ternary plot in a data table to pass to ggtern
plot_data <- data.table(Mean_met, Overall = all_means$V1, Enrich = Enrich)
plot_data[,Enrich:=factor(Enrich, levels = c("None", "Inorganic", "Starved", "Symbiotic"))]
plot_data <- plot_data[order(Enrich)]
col.table <- data.table(gt = c("Inorganic", "Starved", "Symbiotic", "Inorganic+R7A", "None"),
                        color = c("#CC79A7", "#0072B2", "#009E73", "#D95F02", "grey") )

# Ternary plot
ggtern(data = plot_data, aes(x = Symbiotic, y = Inorganic, z = Starved, size = Overall, col = Enrich))+
  geom_point(alpha = 0.5)+
  scale_color_manual(name = "Enriched Feature",
                     breaks = col.table$gt,
                     values = col.table$color)+
  scale_size_continuous(name = "Overall average")+
  theme(text = element_text(size = 20),
        legend.text=element_text(size=20))+
  theme(tern.axis.title.L = element_text(hjust = 0),
        tern.axis.title.R = element_text(hjust = 1))+
  NULL -> tern_plot

# Saving ternary plot
ggsave(filename = paste0("Figures/Ternary/Ternary_class_", as.character(Sys.Date()),".pdf"),
       tern_plot, width = 300, height = 180, units = "mm")

# Saving results shown on ternary plots
fwrite(plot_data, file = paste0("Results/Ternary_plot_table_class_", as.character(Sys.Date()), ".csv"))

# NPC pathway --------------------------------------------------------------------------------------------------------------------
# Setting up feature table collapsed by most superclass
M <- metabolite_data[,-(1:3)]
metabolite_trans <- transpose(M, keep.names = "FeatureID")
metabolite_path <- merge(anno_file[,c("FeatureID", "NPC#pathway")], metabolite_trans, by = "FeatureID")
metabolite_path <- metabolite_path[,-1]
colnames(metabolite_path)[1] <- "Pathway"
metabolite_path <- metabolite_path[,lapply(.SD, sum), Pathway]
metabolite_path <- metabolite_path[Pathway != "unknown"]
colnames(metabolite_path)[-1] <- metabolite_data$Sample_ID
metabolite_path <- transpose(metabolite_path, make.names = "Pathway", keep.names = "Sample_ID")
metabolite_path[,State:=metabolite_data$State]
setcolorder(metabolite_path, c("Sample_ID", "State"))
metabolite_path[,State:=as.factor(State)]

# Testing the difference between states with a linear mixed model
A <- array(dim = c(4,4,ncol(metabolite_path)-2))
dimnames(A)[[1]] <- states
dimnames(A)[[2]] <- states
set.seed(1676559840)
for(i in 1:(ncol(metabolite_path)-2)){
  cat(i, "\r")
  # Metabolite variable
  peak <- metabolite_path[[i+2]]
  # Replacing zero-values with a small value prior to log transformation
  peak[peak == 0] <- 0.5
  peak <- log(peak)
  # Linear mixed models. A model is fitted with each state as a reference level to get all pairwise comparisons
  l_inog <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_path)
  metabolite_path[,State:=relevel(State, ref = "Symbiotic")]
  l_sym <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_path)
  metabolite_path[,State:=relevel(State, ref = "Starved")]
  l_starv <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_path)
  metabolite_path[,State:=relevel(State, ref = "Inorganic+R7A")]
  l_inogR7A <- lmer(peak ~ State + (1|Sample_ID), data = metabolite_path)
  metabolite_path[,State:=relevel(State, ref = "Inorganic")]
  
  # Sets up matrix with p-values for all pairwise comparisons
  states_inog <- setdiff(states, "Inorganic")
  states_inogR7A <- setdiff(states, "Inorganic+R7A")
  states_sym <- setdiff(states, "Symbiotic")
  states_starv <- setdiff(states, "Starved")
  C1 <- coef(summary(l_inog)); rownames(C1) <- gsub("State", "", rownames(C1))
  C2 <- coef(summary(l_inogR7A)); rownames(C2) <- gsub("State", "", rownames(C2))
  C3 <- coef(summary(l_sym)); rownames(C3) <- gsub("State", "", rownames(C3))
  C4 <- coef(summary(l_starv)); rownames(C4) <- gsub("State", "", rownames(C4))
  
  # Inserting results from above into array
  A[,,i][states_inog, "Inorganic"] <- C1[states_inog, 5]
  A[,,i][states_inogR7A, "Inorganic+R7A"] <- C2[states_inogR7A, 5]
  A[,,i][states_sym, "Symbiotic"] <- C3[states_sym, 5]
  A[,,i][states_starv, "Starved"] <- C4[states_starv, 5]
  diag(A[,,i]) <- 1
}

# Correction for multiple testing
B <- array(0, c(4,4,ncol(metabolite_path)-2))
for(i in 1:4){
  for(j in 1:4){
    B[i,j,] <- p.adjust(A[i,j,], method = "fdr")
  }
}

dimnames(B)[[1]] <- states
dimnames(B)[[2]] <- states

# Table with results
res <- data.table(Pathway = colnames(metabolite_path)[-(1:2)],
                  `Starved vs symbiotic` = B["Starved", "Symbiotic",]<0.05,
                  `Starved vs inorganic` = B["Starved", "Inorganic",]<0.05,
                  `Starved vs inorganic+R7A` = B["Starved", "Inorganic+R7A",]<0.05,
                  `Symbiotic vs inorganic` = B["Symbiotic", "Inorganic",]<0.05,
                  `Symbiotic vs inorganic+R7A` = B["Symbiotic", "Inorganic+R7A",]<0.05,
                  `Inorganic vs inorganic+R7A` = B["Inorganic", "Inorganic+R7A",]<0.05)

# Number of significantly different features for each comparison
res[,lapply(.SD, sum), .SDcols = colnames(res)[-1]]

# Setting up tables showing which metabolites are enriched in which conditions
metabolite_path[,lapply(.SD, mean),
                .SDcols = colnames(metabolite_path)[-(1:2)],
                State] -> Means
Means <- transpose(Means, make.names = "State", keep.names = "Pathway")
Means[,":="(
  `Starved vs symbiotic` = ifelse(Starved > Symbiotic, 1, -1),
  `Starved vs inorganic` = ifelse(Starved > Inorganic, 1, -1),
  `Starved vs inorganic+R7A` = ifelse(Starved > `Inorganic+R7A`, 1, -1),
  `Symbiotic vs inorganic` = ifelse(Symbiotic > Inorganic, 1, -1),
  `Symbiotic vs inorganic+R7A` = ifelse(Symbiotic > `Inorganic+R7A`, 1, -1),
  `Inorganic vs inorganic+R7A` = ifelse(Inorganic > `Inorganic+R7A`, 1, -1)
)]

# Constructs matrix with 1, 0, -1 according to if and where a feature is enriched.
# A feature is enriched in state x when "x vs y" is 1 for that feature. If "x vs y" is
# -1, then the feature is enriched in y.
C <- as.matrix(Means[,6:11])
C2 <- as.matrix(res[,-1])
C[!C2] <- 0
Result <- data.table(Pathway = colnames(metabolite_path)[-(1:2)], C)

fwrite(Result, "Results/Pathway_differences.csv")
