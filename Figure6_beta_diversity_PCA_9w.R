
#
# 
#

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions


# directories

setwd("C:/Users/workfolder/")

results.dir <- "C:/Users/workfolder/"
figures.dir <- "C:/Users/workfolder/figures/"


# files

design.file <- paste(results.dir, "design_9w.txt", sep="")
otu_table.file <- paste(results.dir, "ASV_table_9w.txt", sep="")



# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names =1, check.names=F)


# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]


idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]


# conduct Beta diversity


library("ggplot2")
library("scales")
library("grid")
library("vegan")





# normalization and log-transform


otu_table_9w_norm <- apply(otu_table_9w, 2, function(x) x / sum(x))
otu_table_9w_norm_log <- log2(otu_table_9w_norm*1000+1)


### transpose the data frame

library(data.table)

otu_table_9w_norm_log <- as.data.frame(otu_table_9w_norm_log)

t_otu_table_9w_norm_log <- as.data.table(t(as.matrix(otu_table_9w_norm_log)))

t_otu_table_9w_norm_log <- as.data.frame(t_otu_table_9w_norm_log)


colnames(t_otu_table_9w_norm_log) <- rownames(otu_table_9w_norm_log)

rownames(t_otu_table_9w_norm_log) <- colnames(otu_table_9w_norm_log)


SampleID <- rownames(t_otu_table_9w_norm_log)

row.names(t_otu_table_9w_norm_log) <- NULL

df <- cbind(SampleID, t_otu_table_9w_norm_log)


### add compartment, nitrate, inoculum info on df

library(dplyr)

df2 <- data.frame(SampleID = design_9w$SampleID, Compartment = design_9w$Compartment,Inoculum=design_9w$Inoculum, Nitrate_supplement=design_9w$Nitrate_supplement)


df_PCA <- full_join(df, df2, by="SampleID")


row.names(df_PCA)<- df_PCA[,1]

df_PCA$SampleID <- NULL


df_PCA$Compartment <- factor(df_PCA$Compartment, levels = c("rhizosphere","root","control","inoculum"))

# PCA for all

library("FactoMineR")
library("factoextra")



pca <- PCA(df_PCA[,-35:-37], graph = FALSE)


pca

var <- get_pca_var(pca)

fviz_pca_var(pca, col.var = "black")


library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

fviz_cos2(pca, choice = "var", axes = 1:2)


p1 <- fviz_pca_ind(pca,
             geom.ind = "point",
             mean.point = FALSE, pointsize=5,
             col.ind = df_PCA$Compartment,
             palette = c("#d95f02","#1b9e77","#525252","#6a3d9a"),
             addEllipses = TRUE,
             legend.title="Compartment")

p10 <- ggpubr::ggpar(p1, title = "All",
                    xlab = "PC1", ylab = "PC2")


p10



### PCA for the none nitrate samples

idx <- df_PCA$Nitrate_supplement %in% c("none")

df_PCA_none <- df_PCA[idx,]

pca <- PCA(df_PCA_none[,-35:-37], graph = FALSE)


p2 <- fviz_pca_ind(pca,
                   geom.ind = "point", pointsize=5,
                   mean.point = FALSE,
                   col.ind = df_PCA_none$Compartment,
                   palette = c("#d95f02","#1b9e77","#525252","#6a3d9a"),
                   addEllipses = TRUE,
                   legend.title="Compartment")


p11 <- ggpubr::ggpar(p2, title = "None",
                     xlab = "PC1", ylab = "PC2")


p11




### PCA for the nitrate samples


idx <- df_PCA$Nitrate_supplement %in% c("3 mM KNO3")

df_PCA_N <- df_PCA[idx,]

pca <- PCA(df_PCA_N[,-35:-37], graph = FALSE)


p3 <- fviz_pca_ind(pca,
                   geom.ind = "point", pointsize=5,
                   mean.point = FALSE,
                   col.ind = df_PCA_N$Compartment,
                   palette = c("#d95f02","#1b9e77","#525252","#6a3d9a"),
                   addEllipses = TRUE,
                   legend.title="Compartment")


p12 <- ggpubr::ggpar(p3, title = "3 mM KNO3",
                     xlab = "PC1", ylab = "PC2")


p12




### PCA for the rhizosphere


idx <- df_PCA$Compartment %in% c("rhizosphere")

df_PCA_rhizo <- df_PCA[idx,]

pca <- PCA(df_PCA_rhizo[,-35:-37], graph = FALSE)


p4 <- fviz_pca_ind(pca,
                   geom.ind = "point",
                   mean.point = FALSE,
                  
                   fill.ind = df_PCA_rhizo$Nitrate_supplement, col.ind = "black",
                   pointshape =21, pointsize =5, 
                   palette = c("#e6ab02","#d95f02"),
                   addEllipses = TRUE, ellipse.level=0.95, 
                   legend.title="Nitrate_supplement")


p5 <- ggpubr::ggpar(p4, title = "Rhizosphere",
              xlab = "PC1", ylab = "PC2")
              

p5





### PCA for the root


idx <- df_PCA$Compartment %in% c("root")

df_PCA_root <- df_PCA[idx,]

pca <- PCA(df_PCA_root[,-35:-37], graph = FALSE)


p6 <- fviz_pca_ind(pca,
                   geom.ind = "point",
                   mean.point = FALSE,
                   
                   fill.ind = df_PCA_root$Nitrate_supplement, col.ind = "black",
                   pointshape =21, pointsize =5, 
                   palette = c("#66a61e","#1b9e77"),
                   addEllipses = TRUE, ellipse.level=0.95, 
                   legend.title="Nitrate_supplement")


p7 <- ggpubr::ggpar(p6, title = "Root",
                    xlab = "PC1", ylab = "PC2")

p7



# combine the figures into one


library(ggpubr)
library(cowplot)
library(patchwork)

Beta <- p10+p11+p12+p5+p7+plot_layout(guides = 'collect')

Beta

ggsave(paste(figures.dir, "PCA_9w_general.pdf", sep=""), Beta, width=25, height=15,useDingbats=F)

### rhizosphere none nitrate

idx <- df_PCA$Compartment%in% c("rhizosphere")&
        df_PCA$Nitrate_supplement%in% c("none")


df_PCA_rhizo_none <- df_PCA[idx,]

pca <- PCA(df_PCA_rhizo_none[,-35:-37], graph = FALSE)

p13 <- fviz_pca_biplot(pca,
                   mean.point = FALSE,
                   col.ind = df_PCA_rhizo_none$Inoculum,  pointsize =5, 
                   palette = c("#d95f02","#d95f02","#d95f02","#d95f02"),
                   addEllipses = TRUE, ellipse.type="convex",
                   label="var",
                   col.var="black", repel=TRUE)
                   

p14 <- ggpubr::ggpar(p13, title = "Rhizosphere_none", legend = "right")

p14


### rhizosphere nitrate

idx <- df_PCA$Compartment%in% c("rhizosphere")&
  df_PCA$Nitrate_supplement%in% c("3 mM KNO3")

df_PCA_rhizo_N <- df_PCA[idx,]


pca <- PCA(df_PCA_rhizo_N[,-35:-37], graph = FALSE)


p15 <- fviz_pca_biplot(pca,
                       mean.point = FALSE,
                       col.ind = df_PCA_rhizo_N$Inoculum,  pointsize =5, 
                       palette = c("#e6ab02","#e6ab02","#e6ab02","#e6ab02"),
                       addEllipses = TRUE, ellipse.type="convex",
                       label="var",
                       col.var="black", repel=TRUE)


p16 <- ggpubr::ggpar(p15, title = "Rhizosphere_3 mM KNO3", legend = "right")

p16



### root none


idx <- df_PCA$Compartment%in% c("root")&
  df_PCA$Nitrate_supplement%in% c("none")


df_PCA_root_none <- df_PCA[idx,]


pca <- PCA(df_PCA_root_none[,-35:-37], graph = FALSE)


p17 <- fviz_pca_biplot(pca,
                       mean.point = FALSE,
                       col.ind = df_PCA_root_none$Inoculum,  pointsize =5, 
                       palette = c("#1b9e77","#1b9e77","#1b9e77","#1b9e77"),
                       addEllipses = TRUE, ellipse.type="convex",
                       label="var",
                       col.var="black", repel=TRUE)



p18 <- ggpubr::ggpar(p17, title = "Root_none", legend = "right")

p18



### root nitrate

idx <- df_PCA$Compartment%in% c("root")&
  df_PCA$Nitrate_supplement%in% c("3 mM KNO3")


df_PCA_root_N <- df_PCA[idx,]

pca <- PCA(df_PCA_root_N[,-35:-37], graph = FALSE)

p19 <- fviz_pca_biplot(pca,
                       mean.point = FALSE,
                       col.ind = df_PCA_root_N$Inoculum,  pointsize =5, 
                       palette = c("#66a61e","#66a61e","#66a61e","#66a61e"),
                       addEllipses = TRUE, ellipse.type="convex",
                       label="var",
                       col.var="black", repel=TRUE)



p20 <- ggpubr::ggpar(p19, title = "Root_3 mM KNO3", legend = "right")

p20


PCA_root_rhizo <- p14+p16+p18+p20+plot_layout(guides = 'collect')

PCA_root_rhizo


ggsave(paste(figures.dir, "PCA_9w_root_rhizo.pdf", sep=""), PCA_root_rhizo, width=20, height=15,useDingbats=F)
