# Code by Ib Thorsgaard Jensen
library(data.table)
library(magrittr)
library(edgeR)
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Changes working directory to the path of this file. Only works in Rstudio

# Auxiliary helper function----
first.entries <- function(x){
  n <- length(x)
  k <- length(unique(x))
  fe <- rep(NA,k)
  u <- rep(NA,k)
  j <- 1
  for(i in 1:n){
    if(!(x[i] %in% u)){
      fe[j] <- i
      u[j] <- x[i]
      j <- j+1
    }
  }
  return(fe)
}

# Loading datasets----
OTU_table <- fread("../Data/zotu_table_0.01_clean.csv")
meta_data <- fread("../Data/design_clean.csv")
taxonomy <- fread(file="../Data/Taxonomy.txt")

orders_to_display <- c("Xanthomonadales", "Steroidobacterales", "Pseudomonadales",
                       "Burkholderiales", "Sphingomonadales", "Rhizobiales", "Dongiales",
                       "Azospirillales", "Streptomycetales", "Micrococcales", "AsCoM_o_570")

# Sorting by taxonomy
taxonomy <- taxonomy[order(Phylum, Class, Order, Family)]
rearrange <- match(taxonomy$OTUid, OTU_table$OTUid)
OTU_table <- OTU_table[rearrange]

# Family taxonomy and construction of order color-table
first_fam_idx <- first.entries(taxonomy$Family)
fam_tax <- taxonomy[first_fam_idx,c(4,6,8,10)]
sort_by_order <- unique(fam_tax$Order)
ord_colors <- fread(file = "../Data/Colors_order.txt")
ord_colors <- ord_colors[order %in% orders_to_display]
ord_colors <- ord_colors[match(orders_to_display, order)]

# Setting up vector of conditions with entries on the form "compartment_genotype"
comp <- c("root", "rhizosphere")
geno <- c("Gifu", "nfre", "chit5", "nfr5")
conditions <- expand.grid(comp, geno); colnames(conditions) <- c("compartment", "genotype")
conditions <- sapply(conditions, as.character)
conditions <- rbind(c("soil", "soil"), conditions)
cond.names <- apply(conditions, 1, paste0, collapse="_")
n <- nrow(conditions)

# Differential abundance analysis ----
diff.abn.otu.list <- list(); plot_list <- list(); diff.abn.otu.list2 <- list()
plot_data <- data.table(matrix(nrow=0,ncol=10))
colnames(plot_data) <- c("OTU", "logFC", "P_val_adjusted", "Genus", "Family", "Order", "RA", "cond", "s", "fill" )
DA_data_list <- list()
sig.lvl <- 0.05
effect.cutoff <- 1
for(i in 1:n){
  # Subsetting data
  cat(paste(cond.names[i], "          "), "\r")
  gen <- conditions[i, "genotype"]
  com <- conditions[i, "compartment"]
  samples <- meta_data[genotype == gen & compartments == com, SampleID]
  data_subset <- OTU_table[,c("OTUid",samples),with = FALSE]
  meta_data_subset <- meta_data[SampleID %in% samples]; samples <- meta_data_subset$SampleID

  # Differential abundance analysis
  data.DGE <- DGEList(counts=data_subset[,-1], group=meta_data_subset$Nitrate_supplied, genes = data_subset$OTUid)
  data.DGE <- calcNormFactors(data.DGE)
  
  design.matrix <- model.matrix(~ Nitrate_supplied, data = meta_data_subset)
  data.DGE.disp <- estimateDisp(data.DGE, design.matrix)
  
  neg.bin.reg <- glmQLFit(data.DGE.disp, design.matrix)
  neg.bin.reg.test <- glmTreat(neg.bin.reg, coef = "Nitrate_suppliedyes", lfc = effect.cutoff)
  p_adj <- p.adjust(neg.bin.reg.test$table$PValue, method = "BY")
  
  DA_res <- data.table(OTUid = neg.bin.reg.test$genes$genes, neg.bin.reg.test$table, PValueAdjusted = p_adj)
  
  diff.abn.otu <- DA_res[PValueAdjusted<sig.lvl, OTUid]
  
  diff.abn.otu.list[[i]] <- diff.abn.otu
  
  # Collecting information into a data table for plotting
  logFC <- DA_res[PValueAdjusted<sig.lvl, logFC]
  taxonomy_DA <- taxonomy[OTUid %in% diff.abn.otu]
  Genus_DA <- taxonomy_DA$Genus
  Family_DA <- taxonomy_DA$Family
  Order_DA <- taxonomy_DA$Order
  RA <- data.table(Nitrate_supplied = meta_data_subset$Nitrate_supplied, t(data_subset[,-1])/colSums(data_subset[,-1]))
  colnames(RA)[-1] <- data_subset$OTUid
  mean_RA <- RA[, lapply(.SD, mean, na.rm=TRUE), by=Nitrate_supplied, .SDcols=diff.abn.otu ]
  RA_vec_enrich <- apply(mean_RA[,-1], 2, max)
  DA_data <- data.table(OTU = taxonomy_DA$OTUid,
                        logFC,
                        P_val_adjusted = DA_res[PValueAdjusted < sig.lvl, PValueAdjusted],
                        Genus = taxonomy_DA$Genus,
                        Family=Family_DA,
                        Order=Order_DA,
                        RA = RA_vec_enrich,
                        compartment = conditions[i,1],
                        genotype = conditions[i,2],
                        DA_soil = ifelse(taxonomy_DA$OTUid %in% diff.abn.otu.list[[1]], "Same", "EE"))
  if(i != 1){
    # Identifying the OTUs that are enriched in the opposite direction in the soil than in the current condition
    soil_DA_OTU <- DA_data[DA_soil!="EE",OTU]
    soil_lfc <- DA_data_list[[1]][OTU %in% soil_DA_OTU, logFC]
    Soil_enriched_OTU <-  DA_data[DA_soil != "EE"][sign(logFC) != sign(soil_lfc), OTU]
    DA_data[OTU %in% Soil_enriched_OTU, DA_soil:= "Opposite"]
  }
  DA_data[,fill:=ifelse(DA_soil == "Same", "outline", "filled")]
  DA_data[,Order:=factor(Order, levels = unique(Order))]
  DA_data[,Family:=factor(Family, levels = unique(Family))]
  DA_data_list[[i]] <- DA_data
}
names(diff.abn.otu.list) <- cond.names
names(DA_data_list) <- cond.names

plot_data_list <- list()
for(i in 1:4){
  # Constructing a data table with the plotting information for soil, rhizosphere and root
  plot_data <- rbind(DA_data_list[[1]], DA_data_list[[2*i+1]], DA_data_list[[2*i]])
  
  # Setting the factor levels for the families and orders to ensure they are sorted by taxonomy
  plot_data[,Family:=factor(Family, levels = unique(taxonomy[Family %in% plot_data$Family, Family]))]
  plot_data[,Order:=factor(Order, levels = unique(taxonomy[Order %in% plot_data$Order, Order]))]
  
  # Setting the information used for filling the points such that all soil-points are filled
  plot_data[compartment == "soil", DA_soil:="EE"]
  plot_data[compartment == "soil", fill:="filled"]
  
  # Adding size information
  plot_data[,size:=log(RA*50000)]
  
  # Capitalising compartment names for use in facet_wrap
  plot_data[compartment == "soil",compartment:="Soil"]
  plot_data[compartment == "root",compartment:="Root"]
  plot_data[compartment == "rhizosphere",compartment:="Rhizosphere"]
  plot_data[,compartment:=factor(compartment, levels = c("Soil", "Rhizosphere", "Root"))]
  
  plot_data <- plot_data[Order %in% orders_to_display]
  plot_data <- plot_data[,Order:=droplevels(Order)]
  plot_data <- plot_data[,Family:=droplevels(Family)]
  
  plot_data_list[[i]] <- plot_data
}
names(plot_data_list) <- geno

# Defining functions to be used in stat_summary
nitrate_enrich <- function(x){
  return(c(y = -15, label = sum(x<0)))
}
nitrate_depleted <- function(x){
  return(c(y = 15, label = sum(x>0)))
}

# Constructing figures
Combined_plot <- list()
for(i in 1:4){
  ggplot(plot_data_list[[i]], aes(x=logFC, y=Family, color=Order, shape = fill, size = size)) + 
    geom_point()+
    geom_point(data = plot_data_list[[i]][DA_soil == "Same"], shape = 19) +
    geom_point(data = plot_data_list[[i]][DA_soil == "Same"], col = "black", shape = 1) +
    geom_point(data = plot_data_list[[i]][DA_soil == "Opposite"], col = "black", shape = 19) +
    scale_color_manual(breaks = ord_colors$order,
                       values = ord_colors$colors)+
    scale_shape_manual(name = 'In soil',
                       values =c(19,1), labels = c('Depleted', 'Enriched')) +
    scale_size(name = "Mean RA", breaks = log(c(50,500,5000)), labels = c("0.1%", "1%", "10%")) +
    theme(axis.text.x = element_text(vjust=0.5, size = 20),
          axis.text.y = element_text(size = 20),
          text = element_text(size=20),
          panel.grid.major = element_line(color = "grey"),
          panel.grid.minor = element_line(color = "grey"),
          panel.background = element_rect(color = "white", fill ="white"),
          strip.text = element_text(size = 20),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.text = element_text(size = 20)) +
    ylab("Family") +
    xlab("log2FC [Non-nitrate vs nitrate]") + 
    xlim(c(-15,15)) +
    stat_summary(aes(x=logFC, y=Family), fun.data = nitrate_enrich, geom = "text", col = "black", size = 5, inherit.aes = F) +
    stat_summary(aes(x=logFC, y=Family), fun.data = nitrate_depleted, geom = "text", col = "black", size = 5, inherit.aes = F) +
    geom_vline(xintercept = 0, size = 1) +
    guides(size=guide_legend(order = 1), color=guide_legend(order = 2, ncol = 6), shape=guide_legend(order = 3)) +
    facet_wrap(~compartment, ncol = 3) +
    NULL -> Combined_plot[[i]]
}
names(Combined_plot) <- geno

# Saving figures
for(i in 1:4){
  filename <- c("Nitrate", "DA", "plot", geno[i]) %>% paste(collapse = "_") %>% paste("pdf", sep = ".")
  dirname <- c("Figures", filename) %>% paste(collapse = "/")
  ggsave(filename = dirname, plot = Combined_plot[[i]], width = 600, height = 400, units = "mm", scale = 1, useDingbats=FALSE)
}

cond.names[1] <- "soil"
# Removing information from DA_data that is only used internally, and saving the rest for later inspection
for(i in 1:n){
  DA_data <- DA_data_list[[i]]
  DA_data[,c("genotype", "compartment", "DA_soil", "fill"):=NULL]
  setnames(DA_data, "RA", "Mean_RA_enriched_group")
  filename <- c("Nitrate", "DA", "table", cond.names[i]) %>% paste(collapse = "_") %>% paste("csv", sep = ".")
  dirname <- c("Tables", filename) %>% paste(collapse = "/")
  fwrite(DA_data, file = dirname)
}
