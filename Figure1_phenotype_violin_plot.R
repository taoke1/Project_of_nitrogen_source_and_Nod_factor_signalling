options(warn=-1)

# cleanup

rm(list=ls())

# directory

setwd("C:/Users/workfoler")

results.dir <- "C:/Users/workfolder/"
figures.dir <- "C:/Users/workfolder/figures/"

# files

fresh_weight.file <- paste(results.dir, "lysm_soil_N_phenotype.txt", sep = "")

# load data

fresh_weight <- read.table(fresh_weight.file, header=T, sep="\t")

fresh_weight <- as.data.frame(fresh_weight)


### change comma to dot
fresh_weight$shoot_fresh_weight <- gsub("\\,",".",fresh_weight$shoot_fresh_weight) 
fresh_weight$shoot_fresh_weight <- as.numeric(fresh_weight$shoot_fresh_weight)


# subset samples for nodule number visulization


idx <- fresh_weight$Genotype%in%c("Gifu", "nfre","chit5","nfr5")
fresh_weight_subset <- fresh_weight[idx,]


### set levels for dataset

fresh_weight$Genotype <- factor(fresh_weight$Genotype, levels = c("Gifu", "nfre","chit5","nfr5") )

fresh_weight_subset$Genotype <- factor(fresh_weight_subset$Genotype, levels = c("Gifu", "nfre","chit5","nfr5"))


#Plot

## load plotting functions

library("ggplot2")
library("scales")
library("grid")
library(RColorBrewer)


main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major.y = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

## violin plot

dodge <- position_dodge(width = 0.9)

p1 <- ggplot(fresh_weight_subset, aes(x=Genotype, y=shoot_fresh_weight, fill=Genotype)) +
  geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_boxplot( width=0.2,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.7), size=3, alpha=0.3)+
  scale_fill_manual(values= c("#66c2a5","#66c2a5","#66c2a5","#66c2a5" )) +
  main_theme +
  ylab("Shoot fresh weight/plant (g)")+
  scale_y_continuous(breaks=c(0,1.0,2.0,3.0,4.0))+
  scale_x_discrete(labels=expression(Gifu, italic(nfre),italic(chit5),italic(nfr5)))+
  ggtitle(bquote("10 mM" ~KNO[3]))+
  theme(legend.position="none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size= 20),
        plot.title = element_text(size = 20),
        legend.key.size = unit(1,"cm"))


p1

ggsave(paste(figures.dir, "N_Lysm_fresh_weight_volin_Genotype.pdf", sep=""), p1, width=5, height=8)



idx<- fresh_weight_subset$Genotype%in% c("Gifu","nfre","chit5")
total_nods<- fresh_weight_subset[idx,]

p2 <- ggplot(total_nods, aes(x=Genotype, y=total, fill=Genotype)) +
  geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_boxplot( width=0.2,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.7), size=3, alpha=0.3)+
  scale_fill_manual(values= c("#b3b3b3","#b3b3b3","#b3b3b3" )) +
  main_theme +
  ylab("Number of total nodules/plant")+
  scale_y_continuous(breaks=c(0,8,16,24))+
  scale_x_discrete(labels=expression(Gifu, italic(nfre),italic(chit5)))+
  ggtitle(bquote("10 mM" ~KNO[3]))+
  theme(legend.position="none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size= 20),
        plot.title = element_text(size = 20),
        legend.key.size = unit(1,"cm"))


p2

ggsave(paste(figures.dir, "N_Lysm_nods_volin_Genotype.pdf", sep=""), p2, width=4, height=8)



# no nitrate ones

fresh_weight2.file <- paste(results.dir, "lysm_soil_phenotype.txt", sep = "")

### load data

fresh_weight2 <- read.table(fresh_weight2.file, header=T, sep="\t")
fresh_weight2 <- as.data.frame(fresh_weight2)

fresh_weight2$shoot_fresh_weight <- gsub("\\,",".",fresh_weight2$shoot_fresh_weight) 
fresh_weight2$shoot_fresh_weight <- as.numeric(fresh_weight2$shoot_fresh_weight)



### subset data for visualization

idx <- fresh_weight2$Genotype%in%c("Gifu", "nfre","chit5","nfr5")

fresh_weight_subset2 <- fresh_weight2[idx,]


### set levels for dataset

fresh_weight2$Genotype <- factor(fresh_weight2$Genotype, levels = c("Gifu", "nfre","chit5","nfr5") )
fresh_weight_subset2$Genotype <- factor(fresh_weight_subset2$Genotype, levels = c("Gifu", "nfre","chit5","nfr5"))

p3 <- ggplot(fresh_weight_subset2, aes(x=Genotype, y=shoot_fresh_weight, fill=Genotype)) +
  geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_boxplot( width=0.2,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.7), size=3, alpha=0.3)+
  scale_fill_manual(values= c("#66c2a5","#66c2a5","#66c2a5","#66c2a5" )) +
  main_theme +
  ylab("Shoot fresh weight/plant (g)")+
  scale_y_continuous(breaks=c(0,0.6,1.2,1.8,2.4))+
  scale_x_discrete(labels=expression(Gifu, italic(nfre),italic(chit5),italic(nfr5)))+
  ggtitle(paste("none"))+
  theme(legend.position="none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size= 20),
        plot.title = element_text(size = 20),
        legend.key.size = unit(1,"cm"))


p3

ggsave(paste(figures.dir, "Lysm_fresh_weight_volin_Genotype.png", sep=""), p3, width=5, height=8)
ggsave(paste(figures.dir, "Lysm_fresh_weight_volin_Genotype.pdf", sep=""), p3, width=5, height=8)






idx <- fresh_weight_subset2$Genotype%in% c("Gifu","nfre","chit5")
nods <- fresh_weight_subset2[idx,]

p4 <- ggplot(nods, aes(x=Genotype, y=pink_nodules, fill=Genotype)) +
  geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_boxplot( width=0.2,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.7), size=3, alpha=0.3)+
  scale_fill_manual(values= c("#fc8d62","#fc8d62","#fc8d62" )) +
  main_theme +
  ylab("Number of pink nodules/plant")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100))+
  scale_x_discrete(labels=expression(Gifu, italic(nfre),italic(chit5)))+
  ggtitle(paste("none"))+
  theme(legend.position="none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size= 20),
        plot.title = element_text(size = 20),
        legend.key.size = unit(1,"cm"))


p4

ggsave(paste(figures.dir, "N_Lysm_pinknods_volin_Genotype.png", sep=""), p4, width=4, height=8)
ggsave(paste(figures.dir, "N_Lysm_pinknods_volin_Genotype.pdf", sep=""), p4, width=4, height=8)




p5 <- ggplot(nods, aes(x=Genotype, y=total_nodules, fill=Genotype)) +
  geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_boxplot( width=0.2,position = dodge,outlier.color = NA)+
  geom_jitter( aes(group=Genotype), position= position_jitterdodge(jitter.width =0.7), size=3, alpha=0.3)+
  scale_fill_manual(values= c("#b3b3b3","#b3b3b3","#b3b3b3" )) +
  main_theme +
  ylab("Number of total nodules/plant")+
  scale_y_continuous(breaks=c(0,20,40,60,80,100,120))+
  scale_x_discrete(labels=expression(Gifu, italic(nfre),italic(chit5)))+
  ggtitle(paste("none"))+
  theme(legend.position="none", 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size= 20),
        plot.title = element_text(size = 20),
        legend.key.size = unit(1,"cm"))


p5

ggsave(paste(figures.dir, "Lysm_pinknods_volin_Genotype.png", sep=""), p5, width=4, height=8)
ggsave(paste(figures.dir, "Lysm_pinknods_volin_Genotype.pdf", sep=""), p5, width=4, height=8)




### group all the figures


library("cowplot")
library("patchwork")

FW <- plot_grid(p4,p3,p1, labels = c('a','b','c', label_size=20), nrow = 1) 

FW

NN <- plot_grid(p4,p5,p2,labels = c('d','e','f', label_size=20),nrow = 1 )

NN

FW/NN


ggsave(paste(figures.dir, "Fig1_pheno_v6.pdf", sep=""), FW/NN, width=16, height=16)


### statistic analysis for shoot fresh weight in sterile water condition

library(car)

ano <- aov(shoot_fresh_weight ~ Genotype, data=fresh_weight_subset2)
anova(ano) 

### Multiple pairwise-comparsions use pair-wise t test

pairwise <- TukeyHSD(ano)

### Generate lables for significance

library(multcompView)

generate_label_df <- function(pairwise, variable){
  
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}
LABELS=generate_label_df(pairwise , "Genotype")

## output Tukey results and label results

write.csv(LABELS, file =  "lysm_none_fresh_weight_TukeyHSD.csv")




# statistic summarize
library(dplyr)

### summarize sterile water samples
### remove NA values

fresh_weight_subset2 <- fresh_weight_subset2[complete.cases(fresh_weight_subset2[ , 2]),]

shoot_weight_summary <- fresh_weight_subset2%>% 
  group_by(Genotype)%>%
  summarise(Mean=mean(shoot_fresh_weight), Max=max(shoot_fresh_weight), Min=min(shoot_fresh_weight), Median=median(shoot_fresh_weight), Std=sd(shoot_fresh_weight))

shoot_weight_summary


nods <- nods[complete.cases(nods[,4]),]

red_nodule_summary <- nods%>% 
  group_by(Genotype)%>%
  summarise(Mean=mean(pink_nodules), Max=max(pink_nodules), Min=min(pink_nodules), Median=median(pink_nodules), Std=sd(pink_nodules))

total_nodule_summary <- nods%>%
  group_by(Genotype)%>%
  summarise(Mean=mean(total_nodules), Max=max(total_nodules), Min=min(total_nodules), Median=median(total_nodules), Std=sd(total_nodules))

write.csv(rbind(shoot_weight_summary,red_nodule_summary,total_nodule_summary), file="none_pheno_summarize.csv")


### summarize nitrate samples

shoot_weight_summary <- fresh_weight_subset%>%
  group_by(Genotype)%>%
  summarise(Mean=mean(shoot_fresh_weight), Max=max(shoot_fresh_weight), Min=min(shoot_fresh_weight), Median=median(shoot_fresh_weight), Std=sd(shoot_fresh_weight))


total_nods <- total_nods[complete.cases(total_nods[,6]),]

total_nodule_summary <- total_nods%>%
  group_by(Genotype)%>%
  summarise(Mean=mean(total), Max=max(total), Min=min(total), Median=median(total), Std=sd(total))

write.csv(rbind(shoot_weight_summary,total_nodule_summary), file="N_pheno_summarize.csv")


