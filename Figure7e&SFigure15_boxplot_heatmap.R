
#
# 
#

options(warn=-1)

# cleanup

rm(list=ls())




# directories

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

results.dir <- "C:/Users/exudate2023/"
figures.dir <- "C:/Users/exudate2023/Figures/"


# files
library(data.table)

design.file <- paste(results.dir, "exudate_design.txt", sep="")
exudate.file <- paste(results.dir, "filtered_feature_table_norm.txt", sep="")
annotation.file <- paste(results.dir, "anno_canopus.txt", sep = "")



# load data

design <- read.table(design.file, header=T, sep="\t")
exudate <- read.table(exudate.file, sep="\t", header=T, check.names=F)
anno <- read.table(annotation.file, sep = "\t", header = T, check.names = F)

enrichment <- fread("Ternary_plot_table_inog_Feature_2023-03-01_filter.csv")





# re-order data matrices


idx <- design$Sample_ID %in%  exudate$Sample_ID
design <- design[idx, ]

idx <- match(design$Sample_ID, exudate$Sample_ID)
exudate <- exudate[idx,]


idx <- match(colnames(exudate[4:1148]), anno[,1])
anno <- anno[idx, ]

rownames(anno) <- NULL

anno$id <- colnames(exudate[4:1148])

anno[is.na(anno)] <- "unknown"

anno$NPC_pathway[anno$NPC_pathway_Probability < 0.6] <- "low_score0.6"
anno$NPC_class[anno$NPC_class_Probability < 0.6] <- "low_score0.6"



### transpose exudate dataset



exudate_t <- as.data.frame (t(exudate))
colnames(exudate_t) <- exudate_t[1,]
exudate_t <- exudate_t[-(1:3),]
exudate_t[] <- lapply(exudate_t, as.numeric)

### subset Gifu samples
idx <- design$Run%in% c("Gifu")

design_subset <- design[idx,]
exudate_t_subset <- exudate_t[,idx]


exudate_t_subset$anno <- anno$NPC_class

FeatureID <- row.names(exudate_t_subset)
exudate_t_subset$FeatureID <- FeatureID
row.names(exudate_t_subset) <- NULL

### transform wide dataframe to long dataframe


df.long <- reshape(exudate_t_subset, varying = 1:40, direction = "long", idvar = 'FeatureID', timevar="sampleID", v.names=c("PA"), 
                   times=c(colnames(exudate_t_subset[,1:40])), sep="")


row.names(df.long) <- NULL



### summarize at the class level

library(dplyr)

df.long_anno <- df.long %>% select (1,3:4)


df.long_anno <- df.long_anno %>%
  group_by_(.dots = c("anno", "sampleID")) %>%   
  summarise_all(funs(sum))


log_PA <- log2(df.long_anno$PA+1)
df.long_anno$log_PA <- log_PA


### add state info

df2 <- data.frame(sampleID = design_subset$Sample_ID, State=design_subset$State)

df.long_anno <- full_join(df.long_anno, df2, by="sampleID")


### visualization heatmap

library(ggplot2)
library(RColorBrewer)

### order data


idx <- df.long_anno$anno != c("unknown")
df.long_anno <- df.long_anno[idx,]
idx <- df.long_anno$anno!= c("low_score0.6")
df.long_anno <- df.long_anno[idx,]


df.long_anno$State <- factor(df.long_anno$State, levels = c("Inorganic_N + R7A","Symbiotic_N","Inorganic_N","Starved"))

### subset df.long_anno

idx <- df.long_anno$anno %in% c("Phenylethylamines","Simple indole alkaloids","Aminoacids","Aminoglycosides","Monosaccharides","Dicarboxylic acids","Open-chain polyketides","Dipeptides","Fatty acyl homoserine lactones","Aminosugars","Glucosinolates","Pyrimidine nucleos(t)ides","Branched fatty acids","Hydroxy fatty acids","N-acyl amines","Other Octadecanoids","Oxygenated hydrocarbons","Unsaturated fatty acids","Depsidones","Cinnamic acids and derivatives","Gallotannins","Depsides","Flavandiols (Leucoanthocyanidins)","Monomeric stilbenes","Glycerophosphoserines","Jasmonic acids","Anthraquinones and anthrones","Plant xanthones","Flavones","Fatty acyl carnitines","Fatty acyl CoAs","Corynanthe type","Cyclic peptides","Mycosporine and Mycosporine-like amino acids","Thiodiketopiperazine alkaloids","Amino cyclitols","Disaccharides","Polysaccharides","Ascarosides","Fatty acyl glycosides of mono- and disaccharides","Glycerophosphates","Glycerophosphoethanolamines","Leukotrienes","Sophorolipids","Sphingoid bases","Erythromycins","Macrolide lactones","Iridoids monoterpenoids")


df.long_anno <- df.long_anno[idx,]

df.long_anno$anno <- factor(df.long_anno$anno, levels = c("Phenylethylamines","Simple indole alkaloids","Aminoacids","Aminoglycosides","Monosaccharides","Dicarboxylic acids","Open-chain polyketides","Dipeptides","Fatty acyl homoserine lactones","Aminosugars","Glucosinolates","Pyrimidine nucleos(t)ides","Branched fatty acids","Hydroxy fatty acids","N-acyl amines","Other Octadecanoids","Oxygenated hydrocarbons","Unsaturated fatty acids","Depsidones","Cinnamic acids and derivatives","Gallotannins","Depsides","Flavandiols (Leucoanthocyanidins)","Monomeric stilbenes","Glycerophosphoserines","Jasmonic acids","Anthraquinones and anthrones","Plant xanthones","Flavones","Fatty acyl carnitines","Fatty acyl CoAs","Corynanthe type","Cyclic peptides","Mycosporine and Mycosporine-like amino acids","Thiodiketopiperazine alkaloids","Amino cyclitols","Disaccharides","Polysaccharides","Ascarosides","Fatty acyl glycosides of mono- and disaccharides","Glycerophosphates","Glycerophosphoethanolamines","Leukotrienes","Sophorolipids","Sphingoid bases","Erythromycins","Macrolide lactones","Iridoids monoterpenoids"))


### heatmap

df.long_anno$log_PA <- cut(df.long_anno$log_PA, breaks = c(0,4,7,11,15,Inf), right = FALSE)

p1 <- ggplot(df.long_anno) +
  geom_tile(aes(x=sampleID, y= anno, fill = log_PA)) +
  scale_fill_manual(breaks=c("[0,4)","[4,7)","[7,11)", "[11,15)","[15,Inf)"),
                      values = c("#ffffff", "#d7acb4","#c4848f", "#b25c6b", "#8f4451")) +
  facet_wrap(~State,scale="free_x",nrow = 1)+
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold", size = 15),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 20))+
  labs(fill = "intensity")

p1



ggsave(paste(figures.dir, "heatmap.pdf", sep=""), p1, width=15, height=15)



### boxplot

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



p2 <- ggplot(df.long_anno, aes(x=anno, y=PA, fill=State)) +
  geom_boxplot( width=0.6,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=State), position= position_jitterdodge(jitter.width =0.3), size=2, alpha=0.3)+
  facet_wrap(~anno,scale="free_x",nrow = 2)+
  scale_y_log10()+
  labs(x="", y="Compound intensity") +
  main_theme +
  theme(legend.position= "bottom",
        plot.title = element_text(size = 20), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size= 20,angle=45,hjust=0.95))

p2


ggsave(paste(figures.dir, "boxplot_class.pdf", sep=""), p2, width=25, height=25)




### obtain the median value 

library(plyr)


meds <- ddply(df.long_anno,.(anno, State), summarize, med = median(PA), Mean=mean(PA), Max=max(PA), Min=min(PA),Std=sd(PA))

meds <- meds[,1:3]

df_meds <- meds %>%
  gather(key, value, -anno, -State) %>%  
  unite(new.col, c(key, State)) %>%   
  spread(new.col, value)

idx <- match(df_meds$anno, anno[,7])
anno <- anno[idx, ]
df_meds$supclass <- anno$NPC_superclass
df_meds$pathway <- anno$NPC_pathway
  
write.csv(df_meds, file = "median_compound_intensity_class.csv")
  

  
  