## 2019.1.18
## Ke Tao
## host preference index 


#

options(warn=-1)

# cleanup

rm(list=ls())


# directories

setwd("C:/Users/Ke/Desktop/host_preference/")

results.dir <- "C:/Users/Ke/Desktop/host_preference/"
figures.dir <- "C:/Users/Ke/Desktop/host_preference/figures/"

# files

design.file <- paste(results.dir, "design_clean.txt", sep="")
taxonomy.file <- paste(results.dir, "final_taxonomy_clean.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table_0.01_clean.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t",  fill=T)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

idx <- rownames(otu_table) %in% taxonomy[,1]
otu_table <- otu_table[idx, ]


idx <- match(rownames(otu_table), taxonomy[,1])
taxonomy <- taxonomy[idx, ]

rownames(taxonomy) <- NULL




#
# filter 0.3% RA ASVs in Gifu root in native soil condition for the hpi
#


# subset samples

idx <- design$Compartment %in% c("root")& 
  design$Genotype%in% c("Gifu")&
  design$Nitrate_supplement%in% c("10 mM KNO3")
design_subset <- design[idx, ]
otu_table_subset <- otu_table[, idx]



# normlise otu_table_subset

otu_table_subset_norm <- apply(otu_table_subset, 2, function(x) x / sum(x))

#### thresholding the normlized OTU table#####################################

threshold <- .3

idx <- rowSums(otu_table_subset_norm * 100 > threshold) >= 1

otu_table_subset_norm <- otu_table_subset_norm[idx, ]

idx <- row.names(otu_table)%in% row.names(otu_table_subset_norm)

otu_table_subset <- otu_table[idx,]

idx <- design$Compartment%in%c("root","soil")&
  design$Nitrate_supplement%in% c("10 mM KNO3")
otu_table_subset <- otu_table_subset[, idx]

idx <- match(rownames(otu_table_subset), taxonomy[,1])
taxonomy_subset <- taxonomy[idx, ]



# obtain normalized OTU_table for the 141 OTUs


otu_table_norm <- apply(otu_table, 2, function(x) x /sum(x))

idx <- row.names(otu_table_norm)%in% row.names(otu_table_subset)

otu_table_norm_subset <- otu_table_norm [idx, ]

otu_table_norm_subset <- as.data.frame(otu_table_norm_subset)


# subset only for soil and roots samples

idx <- design$Compartment%in%c("root","soil")&
  design$Nitrate_supplement%in% c("none")
otu_table_norm_subset <- otu_table_norm_subset[, idx]
design_subset <- design[idx,]


# add family info

otu_table_norm_subset$Family <- taxonomy_subset[,10]

# specify the Mesorhizobium

otu_table_norm_subset$Genus <- taxonomy_subset[,12]

otu_table_norm_subset$Family[otu_table_norm_subset$Genus=='Mesorhizobium'] <- 'Mesorhizobium'






OTUid <- row.names(otu_table_norm_subset)
row.names(otu_table_norm_subset)<- NULL

df <- cbind(OTUid, otu_table_norm_subset[,1:62])


# reshape data

df.long <- reshape(df, varying = 2:62, direction = "long", idvar = 'OTUid', timevar="sampleID", v.names=c("RA"), 
                   times=colnames(df[,2:62]), sep="")

row.names(df.long) <- NULL

df.long <- as.data.frame(df.long)




#  make a dataframe that sum RA by tax and sampleID

library(dplyr)

df.long_family <- df.long %>% select(2:4)

df.long_family <- df.long_family %>%
  group_by_(.dots = c("Family", "sampleID")) %>%   
  summarise_all(funs(sum))



### transfer long dataframe to wide dataframe


library(tidyr)

df_family <- df.long_family %>%
     spread(sampleID, RA)

family<- df_family$Family

df_family[1] <- NULL

row.names(df_family) <- family

### caculate mean RA for each genotype



idx <- design_subset$Genotype=="soil"
soil_means <- apply(df_family[, idx], 1, mean)

idx <- design_subset$Genotype=="Gifu"
Gifu_means <- apply(df_family[, idx], 1, mean)

idx <- design_subset$Genotype=="nfre"
nfre_means <- apply(df_family[, idx], 1, mean)

idx <- design_subset$Genotype=="chit5"
chit5_means <- apply(df_family[, idx], 1, mean)

idx <- design_subset$Genotype=="nfr5"
nfr5_means <- apply(df_family[, idx], 1, mean)


df_family_m <- data.frame(soil= soil_means, Gifu=Gifu_means, nfre=nfre_means, chit5 = chit5_means, nfr5 = nfr5_means)
family <- row.names(df_family)


df_family_m <- cbind(family, df_family_m)


df.long_family_m <- reshape(df_family_m, varying = 2:6, direction = "long", idvar = 'family', timevar="genotype", v.names=c("RA"), 
                       times=c("soil","Gifu","nfre","chit5","nfr5"), sep="")


row.names(df.long_family_m) <- NULL

df.long_family_m <- df.long_family_m[df.long_family_m$genotype!="soil",]



### strategy 2


Gifu_vs_nfre <- df_family_m$Gifu/df_family_m$nfre
Gifu_vs_chit5 <- df_family_m$Gifu/df_family_m$chit5
Gifu_vs_nfr5 <- df_family_m$Gifu/df_family_m$nfr5

df_hpi2 <- cbind(Gifu_vs_nfre, Gifu_vs_chit5, Gifu_vs_nfr5)

df_hpi2 <- as.data.frame(df_hpi2)

df_hpi2$family <- df_family_m$family

df.long_hpi2 <- reshape(df_hpi2, varying = 1:3, direction = "long", idvar = 'family', timevar="Pairs", v.names=c("hpi"), 
                                              times=c(colnames(df_hpi2[1:3])), sep="")

row.names(df.long_hpi2) <- NULL


### define color

colors <- data.frame(group=c("Gifu_vs_nfre","Gifu_vs_chit5","Gifu_vs_nfr5"), 
                     colors=c("#1b9e77","#1b9e77","#1b9e77"))

colors <- colors[colors$group %in% df.long_hpi2$Pairs, ]




# visualization
library(ggplot2)
library(RColorBrewer)
library(scales)

df.long_hpi2$Pairs <- factor(df.long_hpi2$Pairs, levels = c("Gifu_vs_nfre","Gifu_vs_chit5","Gifu_vs_nfr5"))

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))


p3 <- ggplot(df.long_hpi2, aes(x=hpi, y=Pairs, fill= Pairs))+
  geom_bar(stat="identity", position = "identity")+
  facet_wrap(~family, scales = "free_x", nrow = 4)+
  scale_fill_manual(values=as.character(colors$color)) +
  xlab("Ratio of RA between genotypes")+
  main_theme+
  theme(legend.position="none", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20, angle = 90, hjust=0.95, vjust=0.2),
        legend.key.size = unit(1,"cm"))

p3

ggsave(paste(figures.dir, "Geno_prefer_NGroot0.3_2023.pdf", sep=""), p3, width=15, height=10, useDingbats=FALSE)



### subset for main figure

idx <- df.long_hpi2$family%in% c("Mesorhizobium","Oxalobacteraceae","AsCoM_f_4827")

df.long_hpi2_sub <- df.long_hpi2[idx,]

df.long_hpi2_sub$Pairs <- factor(df.long_hpi2_sub$Pairs, levels = c("Gifu_vs_nfre","Gifu_vs_chit5","Gifu_vs_nfr5"))


p2 <- ggplot(df.long_hpi2_sub, aes(x=hpi, y=Pairs, fill= Pairs))+
  geom_bar(stat="identity",position = "identity")+
  facet_wrap(~family, scales = "free_x", nrow = 1)+
  scale_fill_manual(values=as.character(colors$color)) +
  xlim(-0.2,70)+
  ylab("Genotype preference index")+
  coord_flip() +
  main_theme+
  theme(legend.position="bottom", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, hjust=0.95, vjust=0.2),
        legend.key.size = unit(1,"cm"))

p2

ggsave(paste(figures.dir, "Geno_prefer_sub_Groot0.3_2.pdf", sep=""), p2, width=9, height=10, useDingbats=FALSE)





