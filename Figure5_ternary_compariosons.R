# Code by Ib Thorsgaard Jensen
library(data.table)
library(magrittr)
library(edgeR)
library(ggplot2)
library(ggtern)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Changes working directory to the path of this file. Only works in Rstudio

# Loading datasets
OTU_table <- fread("../Data/zotu_table_0.01_clean.csv")
meta_data <- fread("../Data/design_clean.csv")
taxonomy <- fread(file="../Data/Taxonomy.txt")

# Choose compartment
comp <- "rhizosphere"

orders_to_display <- c("Burkholderiales", "Rhizobiales", "Streptomycetales", "Sphingomonadales",
                       "Micrococcales", "Pseudomonadales", "Xanthomonadales",
                       "Azospirillales", "Steroidobacterales", "Acholeplasmatales")

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

# Obtaining the mean RA of each OTU within each group
RA <- transpose(data_subset, make.names = "OTUid")
RA[, 1:ncol(RA) := RA %>% apply(1, function(x) x/sum(x) ) %>% t() %>% as.data.table()]
RA[,State:=meta_data$State]
mean_RAs <- RA[, lapply(.SD, mean), by=State]

# Obtaining the mean between groups
RA_grand_mean <- colMeans(mean_RAs[,-1])

# Setting the significance level and effect-size cutoff
effect.cutoff <- 1
sig.lvl <- 0.05

# Setting up DGE object and fitting negative binomial model with edgeR
data.DGE <- DGEList(counts=data_subset[,-1], genes = data_subset$OTUid)
data.DGE <- calcNormFactors(data.DGE)

design.matrix <- model.matrix(~ 0 + State, data = meta_data)
data.DGE.disp <- estimateDisp(data.DGE, design.matrix)
neg.bin.reg <- glmQLFit(data.DGE.disp, design.matrix)

# Testing for differential abundance between states. Using BY-correction rather 
# than BH to compensate for the lower between-group variance due to technical replicates
enrich_inog_sym_test <- glmTreat(neg.bin.reg, contrast = c(1, 0, -1), lfc = effect.cutoff) %>% 
  decideTestsDGE(adjust.method = "BY", p.value = sig.lvl)
enrich_inog_starv_test <- glmTreat(neg.bin.reg, contrast = c(1, -1, 0), lfc = effect.cutoff) %>% 
  decideTestsDGE(adjust.method = "BY", p.value = sig.lvl)
enrich_sym_starv_test <- glmTreat(neg.bin.reg, contrast = c(0, -1, 1), lfc = effect.cutoff) %>% 
  decideTestsDGE(adjust.method = "BY", p.value = sig.lvl)

# OTUs enriched in the inorganic state compared to the other two
enrich_inog_sym <- neg.bin.reg$genes$genes[enrich_inog_sym_test == 1]
enrich_inog_starv <- neg.bin.reg$genes$genes[enrich_inog_starv_test == 1]
enrich_inog <- intersect(enrich_inog_sym, enrich_inog_starv)

# OTUs enriched in the symbiotic state compared to the other two
enrich_sym_inog <- neg.bin.reg$genes$genes[enrich_inog_sym_test == -1]
enrich_sym_starv <- neg.bin.reg$genes$genes[enrich_sym_starv_test == 1]
enrich_sym <- intersect(enrich_sym_inog, enrich_sym_starv)

# OTUs enriched in the starved state compared to the other two
enrich_starv_inog <- neg.bin.reg$genes$genes[enrich_inog_starv_test == -1]
enrich_starv_sym <- neg.bin.reg$genes$genes[enrich_sym_starv_test == -1]
enrich_starv <- intersect(enrich_starv_inog, enrich_starv_sym)

# Setting up vector denoting which OTUs are enriched in which state, if any
enrich <- fcase(
  !(data_subset$OTUid %in% c(enrich_inog, enrich_sym, enrich_starv)), "None",
  data_subset$OTUid %in% enrich_inog, "Inorganic",
  data_subset$OTUid %in% enrich_sym, "Symbiotic",
  data_subset$OTUid %in% enrich_starv, "Starved"
)
names(enrich) <- data_subset$OTUid

# Constructing data set to be passed to ggplot2
mean_RAs <- transpose(mean_RAs, keep.names = "OTUid", make.names = "State")
plot_data <- cbind(mean_RAs, Enriched = enrich, size = RA_grand_mean)

# Reorder dataset such that the OTUs that are not enriched in any state are in the first rows. 
# This ensures that they will be plotted in the bottom-most layer by ggplot
plot_data <- plot_data[order(Enriched)]
plot_data <- rbind(plot_data[Enriched == "None"], plot_data[Enriched != "None"])

# Saving table with enriched OTUs for later inspection
filename <- c("State_Associated_OTUs", comp) %>% paste(collapse = "_") %>% paste("csv", sep = ".")
dirname <- c("Results", "Nutrtion_associated_bacteria", filename) %>% paste(collapse = "/")
fwrite(plot_data[Enriched != "None", c("OTUid", "Enriched")], file = dirname)

# Table with color information
col.table <- data.table(gt = c("Inorganic", "Starved", "Symbiotic", "None"),
                        color = c("#CC79A7", "#0072B2", "#009E73", "grey") )

# Ternary plots
ggtern(data = plot_data, aes(x=Symbiotic, y=Inorganic, z=Starved, colour = Enriched))+
  geom_point(size = log(40000*plot_data[,size]+1))+
  scale_color_manual(name = "Enriched OTUs",
                     breaks = col.table$gt,
                     values = col.table$color)+
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA))+
  theme(legend.text.align = 0,
        text = element_text(size=20),
        legend.text=element_text(size=20)) -> fig

# Saving ternary plot
filename <- c("ternary", "init", comp) %>% paste(collapse = "_") %>% paste("pdf",sep = ".")
dirname <- c("Figures","Ternary plots", filename) %>% paste(collapse = "/")
ggsave(filename = dirname, plot = fig, width = 40, height = 25, units = "cm", scale = 1, useDingbats=FALSE)

## Pie charts of order-level information on enriched OTUs ----
enrich_list <- list(Chem = enrich_inog, Sym = enrich_sym, Starv = enrich_starv)
tax_list <- list()
for(i in 1:3){
  a <- table(taxonomy[OTUid %in% enrich_list[[i]], Order])
  tax_list[[i]] <- data.table(Order=names(a), Percentage = as.numeric(a)/sum(a), abs = as.numeric(a))
}

# Color table
ord_colors <- fread("../Data/Colors_order.txt")
orders_to_display <- c(orders_to_display, "Other")
color.table <- ord_colors[order %in% orders_to_display]
color.table <- color.table[match(orders_to_display, order)]

# Collapses all orders that are not in orders_to_display. These make up
# at least 5% of the orders in at least one state in at least one compartment
collapse_order <- function(X){
  abs_other <- sum(X[!(Order %in% orders_to_display), abs])
  X <- X[Order %in% orders_to_display]
  Y <- rbind(X, data.table( Order = "Other",
                            Percentage = 1-sum(X$Percentage),
                            abs = abs_other ) )
  nms <- orders_to_display[!(orders_to_display %in% Y$Order)]
  if(length(nms)>0){
    Y <- rbind(Y, data.table(Order = nms, Percentage = 0, abs = 0))
  }
  perm <- color.table[order %in% Y$Order, order]
  idx <- match(perm, Y$Order)
  Y <- Y[idx,]
  Y[,Order:=factor(Order, levels = orders_to_display)]
  return(Y)
}

orders <- lapply(tax_list, collapse_order)

# Plotting
group.names <- c("Inorganic", "Symbiotic", "Starved")
pie_charts <- list()
for(i in 1:3){
  ggplot(data=orders[[i]], aes(x = "", y=Percentage, fill = Order)) +
    geom_bar(stat="identity") + 
    scale_fill_manual(breaks = color.table$order, values = color.table$color, name = "Order")+
    labs(title = paste(group.names[i], "enriched taxonomy", sum(orders[[i]]$abs))) +
    coord_polar("y", start=0) + 
    theme(plot.title=element_text(size=22, hjust = 0.5),
          panel.background = element_rect(color = "white", fill ="white"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          text = element_text(size = 20),
          legend.text.align = 0,
          legend.text=element_text(size=20))+
    NULL -> pie_charts[[i]]
}

gg <- grid.arrange(pie_charts[[1]], pie_charts[[2]], pie_charts[[3]], nrow = 2)

filename <- c("Pie", "chart", comp) %>% paste(collapse = "_") %>% paste("pdf",sep = ".")
dirname <- c("Figures","Pie Charts", filename) %>% paste(collapse = "/")
ggsave(filename = dirname, plot = gg, width = 40, height = 33, units = "cm", scale = 1, useDingbats=FALSE)

