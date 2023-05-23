options(warn=-1)

# cleanup

rm(list=ls())

# directory

setwd("C:/Users/workfolder/")

results.dir <- "C:/Users/workfolder/"
figures.dir <- "C:/Users/workfolder/figures/"

# files

length.file <- paste(results.dir, "shoot_length_9w.txt", sep = "")

# load data

length <- read.table(length.file, header=T, sep="\t")

length <- as.data.frame(length)

# set color for the Inoculum
# set order

length$Inoculum <- factor(length$Inoculum, levels = c("SC+R7A","SC+R7AnifH","SC+R7AnodC","SC","uninoculated"))

length$Nitrate_supplement <- factor(length$Nitrate_supplement, levels = c("3 mM KNO3","none"))

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library(RColorBrewer)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="top",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))

## violin plot

dodge <- position_dodge(width = 0.8)


# 9 weeks length

p1 <- ggplot(length, aes(x=Inoculum, y=Length, fill=Nitrate_supplement)) +
  geom_violin(trim=FALSE,position = dodge,scale = "width",alpha=0.3,color=NA) +
  geom_jitter( position=position_jitterdodge(jitter.width = 0.3), size=3, alpha=0.5)+
  geom_boxplot( width=0.1,position = dodge,outlier.color = NA)+
  facet_wrap(~Nitrate_supplement, scales = "free_x", nrow = 1)+
  scale_fill_manual(values=c("#66a61e","#1b9e77")) +
  main_theme +
  ggtitle("9wpi")+
  ylab("Shoot Length/plant(cm)")+
  theme(legend.position="none", 
        plot.title = element_text(size = 20, face="bold"), 
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90,hjust=0.95,vjust=0.2),
        legend.key.size = unit(1,"cm"))



p1

ggsave(paste(figures.dir, "shootlength_9w.png", sep=""), p1, width=12, height=8)
ggsave(paste(figures.dir, "shootlength_9w.pdf", sep=""), p1, width=12, height=8, useDingbats=F)




# Statistical analysis
# Order the dataframe
length_9week<-as.data.frame(length_9week)
length_9week$Inoculum <- factor(length_9week$Inoculum, levels = colors$group)

library(car)

ano <- aov(Length ~ Inoculum*Nitrate, data=length_9week)

anova(ano) ### shows significant difference

### Multiple pairwise-comparsions use pair-wise t test

pairwise <- TukeyHSD(ano)

### Generate lables for significance

library(multcompView)

generate_label_df <- function(pairwise, variable){
  
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  return(Tukey.labels)
}


pairwise[["Inoculum:Nitrate"]] <- na.omit(pairwise[["Inoculum:Nitrate"]])

# Generate the significant labels for each of the sample
LABELS=generate_label_df(pairwise , "Inoculum:Nitrate")

## output Tukey results and label results
write.csv(LABELS, file =  "sctest4_9wpi_TukeyHSD.csv")

