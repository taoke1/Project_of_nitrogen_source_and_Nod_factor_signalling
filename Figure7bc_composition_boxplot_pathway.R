
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


exudate_t_subset$anno <- anno$NPC_pathway

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

# visualization 

## boxplot

library(ggplot2)
library(RColorBrewer)

### order data


idx <- df.long_anno$anno != c("unknown")
df.long_anno <- df.long_anno[idx,]
idx <- df.long_anno$anno!= c("low_score0.6")
df.long_anno <- df.long_anno[idx,]

colors <- data.frame(group=c("Alkaloids","Amino acids and Peptides","Carbohydrates","Fatty acids","Polyketides","Shikimates and Phenylpropanoids","Terpenoids"), 
                     colors=c("#9c4a59","#b9808a","#d7b6bc","#679f62","#49633c","#5990b7", "#00538a"))
colors <- colors[colors$group %in% df.long_anno$anno, ]
df.long_anno$anno <- factor(df.long_anno$anno, levels = colors$group)


df.long_anno$State <- factor(df.long_anno$State, levels =c("Inorganic_N + R7A","Symbiotic_N","Inorganic_N","Starved"))


dodge <- position_dodge(width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="bottom",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))



p1 <- ggplot(df.long_anno, aes(x=State, y=PA, fill=anno)) +
  geom_boxplot( width=0.6,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=State), position= position_jitterdodge(jitter.width =0.3), size=2, alpha=0.3)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~anno,scale="free_x",nrow = 1)+
  scale_y_log10()+
  labs(x="", y="Metabolite intensity") +
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


ggsave(paste(figures.dir, "boxplot_pathway.pdf", sep=""), p1, width=20, height=12)


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



### conduct donut chart on average intensities


Inor_anno <- Inorganic %>% select(4,6)
Inor_anno <- Inor_anno %>%
  group_by_(.dots = c("Pathway")) %>%   
  summarise_all(funs(sum))


Inor_anno$fraction = Inor_anno$Inorganic/sum(Inor_anno$Inorganic)
Inor_anno$ymax = cumsum(Inor_anno$fraction)
Inor_anno$ymin = c(0, head(Inor_anno$ymax, n= -1))

Inor_anno$Pathway <- factor(Inor_anno$Pathway, levels = colors$group)

p2 <- ggplot(Inor_anno, aes(ymax=ymax, ymin=ymin, xmax=2, xmin=1, fill=Pathway)) +
      geom_rect() +
  coord_polar("y")+
  scale_fill_manual(values = colors$colors)+
  xlim(c(-1, 4)) +
  ggtitle("Inorganic")+
  theme_void() +
  theme(legend.position="none", 
        plot.title = element_text(size = 20, face="bold"), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(1,"cm"))
               

p2


ggsave(paste(figures.dir, "donut_inorganic_features.pdf", sep=""), p2, width=8, height=8)

### starved

Starved$Pathway <- star_anno$NPC_pathway


Star_anno <- Starved %>% select(2,6)
Star_anno <- Star_anno %>%
  group_by_(.dots = c("Pathway")) %>%   
  summarise_all(funs(sum))


Star_anno$fraction = Star_anno$Starved/sum(Star_anno$Starved)
Star_anno$ymax = cumsum(Star_anno$fraction)
Star_anno$ymin = c(0, head(Star_anno$ymax, n= -1))

Star_anno$Pathway <- factor(Star_anno$Pathway, levels = colors$group)

p3 <- ggplot(Star_anno, aes(ymax=ymax, ymin=ymin, xmax=2, xmin=1, fill=Pathway)) +
  geom_rect() +
  coord_polar("y")+
  scale_fill_manual(values = colors$colors)+
  xlim(c(-1, 4)) +
  ggtitle("Starved")+
  theme_void() +
  theme(legend.position="none", 
        plot.title = element_text(size = 20, face="bold"), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(1,"cm"))


p3

ggsave(paste(figures.dir, "donut_Starved_features.pdf", sep=""), p3, width=8, height=8)



### Symbiotic

Symbiotic$Pathway <- sym_anno$NPC_pathway


Sym_anno <- Symbiotic %>% select(3,6)
Sym_anno <- Sym_anno %>%
  group_by_(.dots = c("Pathway")) %>%   
  summarise_all(funs(sum))


Sym_anno$fraction = Sym_anno$Symbiotic/sum(Sym_anno$Symbiotic)
Sym_anno$ymax = cumsum(Sym_anno$fraction)
Sym_anno$ymin = c(0, head(Sym_anno$ymax, n= -1))

Sym_anno$Pathway <- factor(Sym_anno$Pathway, levels = colors$group)

p4 <- ggplot(Sym_anno, aes(ymax=ymax, ymin=ymin, xmax=2, xmin=1, fill=Pathway)) +
  geom_rect() +
  coord_polar("y")+
  scale_fill_manual(values = colors$colors)+
  xlim(c(-1, 4)) +
  ggtitle("Symbiotic")+
  theme_void() +
  theme(legend.position="none", 
        plot.title = element_text(size = 20, face="bold"), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(1,"cm"))


p4

ggsave(paste(figures.dir, "donut_Symbiotic_features.pdf", sep=""), p4, width=8, height=8)




  