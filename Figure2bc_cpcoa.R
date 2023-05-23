
#
# modified from Ruben Garrido-Oter https://github.com/garridoo/lotus/blob/master/statistical_analyses/cpcoa.R
# garridoo@mpipz.mpg.de
#

options(warn=-1)

# cleanup

rm(list=ls())

# set working directory
setwd("C:/Users/workfolder/")

results.dir <- "C:/Users/workfolder/"
figures.dir <- "C:/Users/workfolder/figures/"

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("vegan")

# load plotting functions directory

source("cpcoa.func.R")

# files

design.file <- paste(results.dir, "design_clean.txt", sep="")
taxonomy.file <- paste(results.dir, "final_taxonomy_clean.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table_0.01_clean.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names =1, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

idx <- rownames(otu_table) %in% taxonomy[, 1]
otu_table <- otu_table[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# choose compartments

idx <- design$Compartment %in% c( "rhizosphere")


design <- design[idx, ]
otu_table <- otu_table[, idx]

# normalize otu table

otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))

# CPCoA

colors <- data.frame(group=c("Gifu","nfr5","nfre","chit5"), 
                     colors=c("#1b9e77","#666666","#d95f02","#7570b3")) ### color used here is from "Dark2"

shapes <- data.frame(group=c("none","10 mM KNO3"),
                     shape=c(18, 17))


colors <- colors[colors$group %in% design$Genotype, ]
shapes <- shapes[shapes$group %in% design$Nitrate_supplement,]


capscale.gen <- capscale(t(otu_table_norm) ~ Genotype*Nitrate_supplement + Condition(bio_replicate * tech_replicate), data=design, add=F, sqrt.dist=T, distance="bray")


# ANOVA-like permutation analysis

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)

# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

# extract the weighted average (sample) scores

points <- capscale.gen$CCA$wa[, 1:2]
points <- as.data.frame(points)

colnames(points) <- c("x", "y")

points <- cbind(points, design[match(rownames(points), design$SampleID), ])

points$Genotype <- factor(points$Genotype, levels=colors$group)
points$Nitrate_supplement <- factor(points$Nitrate_supplement, levels=shapes$group)

# plot CPCo 1 and 2



main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=15),
                    legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

p6 <- ggplot(points, aes(x=x, y=y, color=Genotype,shape=Nitrate_supplement)) +
  stat_ellipse(type = "norm", level = 0.8)+
  geom_point(size=8, alpha=0.7) +
  scale_colour_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape)+
  labs(x=paste("CPCo 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCo 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
  ggtitle(paste(format(100 * variance, digits=3), " % of variance; p=",
                format(p.val, digits=2),
                sep="")) +
  main_theme +
  theme(legend.position="bottom", 
        plot.title = element_text(size = 20, face="bold"), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

p6

ggsave(paste(figures.dir, "NF_cpcoa_Genotype_nitrate_root_ellipse.png", sep=""), p6, width=10, height=8)
ggsave(paste(figures.dir, "NF cpcoa_Genotype_nitrate_root_ellipse.pdf", sep=""), p6, width=7, height=8,useDingbats=F)








