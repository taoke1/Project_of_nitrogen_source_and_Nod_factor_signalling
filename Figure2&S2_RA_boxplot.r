## 2019.1.18
## Ke Tao
## stack bar plot for relative abundance


### 
options(warn=-1)

# cleanup

### 
rm(list=ls())


# directories

###
setwd("C:/Users/workfoler/")

results.dir <- "C:/Users/workfolder/"
figures.dir <- "C:/Users/workfolder/figures/"


# files

###
design.file <- paste(results.dir, "design_clean.txt", sep="")
###
taxonomy.file <- paste(results.dir, "final_taxonomy_clean.txt", sep="")
###
otu_table.file <- paste(results.dir, "otu_table_0.01_clean.txt", sep="")

# load data

###
design <- read.table(design.file, header=T, sep="\t")
###
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
###
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

# re-order data matrices

###
idx <- design$SampleID %in% colnames(otu_table)
###
design <- design[idx, ]

###
idx <- match(design$SampleID, colnames(otu_table))
###
otu_table <- otu_table[, idx]

###
idx <- rownames(otu_table) %in% taxonomy[,1]
###
otu_table <- otu_table[idx, ]

###
idx <- match(rownames(otu_table), taxonomy[,1])

###
taxonomy <- taxonomy[idx, ]



# otu_table normalization

###
otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))

df <- as.data.frame(otu_table_norm)


# add taxonomy info to otu_table_norm

tax_order <- taxonomy[ ,8]

df$order <- tax_order


OTUid <- row.names(df)

df <- cbind(OTUid, df)

row.names(df) <- NULL



### reshape data into long dataframe

library(reshape2)

df.long <- reshape(df, varying = 2:180, direction = "long", idvar = 'OTUid', timevar="sampleID", v.names=c("RA"), 
times=c(colnames(df[,2:180])), sep="")

row.names(df.long) <- NULL

df.long <- as.data.frame(df.long)



#  make a dataframe that sum RA by tax and sampleID

library(dplyr)

df.long_order <- df.long %>% select (2:4)


df.long_order <- df.long_order %>%
  group_by_(.dots = c("order", "sampleID")) %>%   ### this is trying to summarize counts by both tax and sampleID
  summarise_all(funs(sum))



# Add compartment and genotype information for df.long

library(dplyr)

df2 <- data.frame(sampleID = design$SampleID, Nitrate_supplement = design$Nitrate_supplement, Compartment = design$Compartment,Genotype=design$Genotype,Allele=design$Allele, Bio_replicate=design$bio_replicate)

df.long_order <- full_join(df.long_order, df2, by="sampleID")




# boxplot using ggplot2

library(ggplot2)
library(RColorBrewer)
library(ggh4x)




### replace the name of tax that are not included in the top

library(stringr)

do_not_replace <- c("Rhizobiales","Burkholderiales","Caulobacterales","Sphingomonadales","Streptomycetales","Pseudomonadales","Xanthomonadales","Azospirillales","Dongiales","Steroidobacterales","AsCoM_o_570","Micrococcales")

df.long_order$order <- ifelse(df.long_order$order%in%do_not_replace, df.long_order$order, "Other")



### rhizosphere samples

idx <- df.long_order$Compartment%in% c("rhizosphere")

df.long_rhizo <- df.long_order[idx,]



### plot

### set color

colors <- data.frame(group=c("Gifu","nfre","chit5","nfr5"), 
                     colors=c("#1b9e77","#d95f02","#7570b3","#666666"))
colors <- colors[colors$group %in% df.long_rhizo$Genotype, ]

### set level

df.long_rhizo$Genotype <- factor(df.long_rhizo$Genotype, levels =colors$group )
df.long_rhizo$Nitrate_supplement <- factor(df.long_rhizo$Nitrate_supplement, levels = c("none", "10 mM KNO3"))
df.long_rhizo$order <- factor(df.long_rhizo$order, levels = c("Streptomycetales","Micrococcales","Azospirillales","AsCoM_o_570","Dongiales","Sphingomonadales","Rhizobiales","Caulobacterales","Pseudomonadales","Steroidobacterales","Burkholderiales","Xanthomonadales","Other"))



library(ggplot2)

dodge <- position_dodge(width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))



p1 <- ggplot(df.long_rhizo, aes(x=order, y=RA, fill=Genotype)) +
  geom_boxplot( width=0.6,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.3), size=3, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~Nitrate_supplement,scale="free_x")+
  labs(x="", y="Aggregated relative abundance") +
  ggtitle("Rhizosphere")+
  main_theme +
  theme(legend.position= "bottom",
        plot.title = element_text(size = 20), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size= 20,angle=45,hjust=0.95))

p1


ggsave(paste(figures.dir, "RA_boxplot_topOrder_rhizo.pdf", sep=""), p1, width=20, height=8, useDingbats=F)




### root samples

idx <- df.long_order$Compartment%in% c("root")

df.long_root <- df.long_order[idx,]



### plot

### set color

colors <- data.frame(group=c("Gifu","nfre","chit5","nfr5"), 
                     colors=c("#1b9e77","#d95f02","#7570b3","#666666"))
colors <- colors[colors$group %in% df.long_root$Genotype, ]

### set level

df.long_root$Genotype <- factor(df.long_root$Genotype, levels =colors$group )
df.long_root$Nitrate_supplement <- factor(df.long_root$Nitrate_supplement, levels = c("none", "10 mM KNO3"))
df.long_root$order <- factor(df.long_root$order, levels = c("Streptomycetales","Micrococcales","Azospirillales","AsCoM_o_570","Dongiales","Sphingomonadales","Rhizobiales","Caulobacterales","Pseudomonadales","Steroidobacterales","Burkholderiales","Xanthomonadales","Other"))



library(ggplot2)

dodge <- position_dodge(width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))



p2 <- ggplot(df.long_root, aes(x=order, y=RA, fill=Genotype)) +
  geom_boxplot( width=0.6,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.3), size=3, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~Nitrate_supplement,scale="free_x")+
  labs(x="", y="Aggregated relative abundance") +
  ggtitle("Root")+
  main_theme +
  theme(legend.position= "bottom",
        plot.title = element_text(size = 20), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size= 20,angle=45,hjust=0.95))

p2


ggsave(paste(figures.dir, "RA_boxplot_topOrder_root.pdf", sep=""), p2, width=20, height=8, useDingbats=F)



### stack bar plot


### rhizosphere
### generate vector of colors

colors <- read.table("colors_order.txt", sep="\t", header=T)
taxon <- sort(unique(df.long_rhizo$order))
colors <- data.frame(taxon=taxon, colors=colors$color[match(taxon, colors$order)])

colors$colors <- as.character(colors$colors)


df.long_rhizo$order <- factor(df.long_rhizo$order, levels = colors$taxon)


p3 <- ggplot(df.long_rhizo, aes(x=sampleID, y = RA, fill = order)) +
  geom_bar(stat = "identity", width = .5) +
  facet_nested(~Nitrate_supplement+Genotype,scales ="free_x")+
  scale_fill_manual(values = colors$colors)+
  main_theme+
  ylab("Relative abundance")+
  ggtitle("Rhizosphere")+
  theme(legend.position = "bottom")+ guides(fill=guide_legend(nrow=2))+ 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=15),
        legend.title = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20))


p3


ggsave(paste(figures.dir, "RA_barplot_topOrder_rhizo.pdf", sep=""), p3, width=20, height=8, useDingbats=F)


### root


colors <- read.table("colors_order.txt", sep="\t", header=T)
taxon <- sort(unique(df.long_root$order))
colors <- data.frame(taxon=taxon, colors=colors$color[match(taxon, colors$order)])

colors$colors <- as.character(colors$colors)


df.long_root$order <- factor(df.long_root$order, levels = colors$taxon)


p4 <- ggplot(df.long_root, aes(x=sampleID, y = RA, fill = order)) +
  geom_bar(stat = "identity", width = .5) +
  facet_nested(~Nitrate_supplement+Genotype,scales ="free_x")+
  scale_fill_manual(values = colors$colors)+
  main_theme+
  ylab("Relative abundance")+
  ggtitle("Root")+
  theme(legend.position = "bottom")+ guides(fill=guide_legend(nrow=2))+ 
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=15),
        legend.title = element_blank(),
        plot.title = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20))


p4


ggsave(paste(figures.dir, "RA_barplot_topOrder_root.pdf", sep=""), p4, width=20, height=8, useDingbats=F)




