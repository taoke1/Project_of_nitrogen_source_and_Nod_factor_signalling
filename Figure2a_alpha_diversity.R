#
#2020.11.30 
#Ke Tao
#


alpha.file <- paste(results.dir, "alpha.txt", sep="")
design.file <- paste(results.dir, "design_clean.txt", sep="")

# load data


design <- read.table(design.file, header=T, sep="\t")
alpha <- read.table(alpha.file, sep="\t", header=T, check.names=F, row.names = 1)



### alpha diversity

shapes <- data.frame(group=c("soil","Gifu", "nfre","chit5","nfr5"), shape=c(8,19,0, 2, 15))

colors <- data.frame(group=c("soil","rhizosphere","root","nodules"), color=c("#666666","#d95f02","#1b9e77","#e7298a"))


# Chao1 index

index <- cbind(alpha[, 2], design[match(row.names(alpha), design$SampleID), ])
colnames(index)[1] <- "value"


# boxplots by Compartment

library(ggplot2)

index <- index[complete.cases(index),]

l1 <- c("soil", "rhizosphere", "root", "nodules")
index$Compartment <- factor(index$Compartment, levels=l1)
colors <- colors[match(l1, colors$group), ]
l2 <- c("soil","Gifu", "nfre","chit5","nfr5")
index$Genotype <- factor(index$Genotype, levels = l2)
shapes <- shapes[match(l2, shapes$group),]


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


p15 <- ggplot(index, aes(x=Nitrate_supplement, y=value, fill=Compartment)) +
  geom_boxplot(alpha=0.7, position= position_dodge(width = 0.7), outlier.color =NA, width=0.3) +
  geom_jitter(aes(shape=Genotype), position=position_jitter(0.17), size=5, alpha=1) +
  scale_fill_manual(values=as.character(colors$color)) +
  scale_shape_manual(values=shapes$shape) +
  facet_wrap(~index$Compartment,scale="free_x",nrow = 1)+
  ###scale_y_continuous(limits = c(0, 12))+
  labs(x="", y="Chao1 index") +
  main_theme +
  theme(legend.position="right", 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust=1),
        axis.title.y = element_text(size = 20),
        legend.key.size = unit(1,"cm"))





p15
ggsave(paste(figures.dir, "Chao1_lysm_Compartment_allASVs_v2.pdf", sep=""), p15, width=15, height=8)





# subset data

idx <- index$Compartments%in% c("soil","root", "rhizosphere")
index <- index[idx,]



l1 <- c("soil", "rhizosphere", "root", "nodules")
index$Compartments <- factor(index$Compartments, levels=l1)


l2 <- c("soil","Gifu", "nfre","chit5","nfr5")
index$Genotype <- factor(index$Genotype, levels = l2 )

l3 <- c("none","10 mM KNO3")
index$Nitrate_supplement <- factor(index$Nitrate_supplement, levels = l3)




# box plot

library(ggplot2)
library(ggh4x)
library(RColorBrewer)

p14 <- ggplot(index, aes(x=Genotype, y=value, fill=Compartments)) +
  geom_boxplot(alpha=0.7, position= position_dodge(width = 0.7), outlier.color =NA, width=0.3) +
  geom_jitter( position=position_jitter(0.17), size=5, alpha=0.3) +
  scale_fill_manual(values=c("#666666","#d95f02","#1b9e77")) +
  facet_nested(~Compartments+Nitrate_supplement,scales ="free_x")+
  scale_x_discrete(labels=expression(Gifu, italic(nfre),italic(chit5),italic(nfr5)))+
  labs(x="", y="Chao1 index") +
  main_theme +
  theme(legend.position="none", 
        plot.title = element_text(size = 20, face="bold"), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_text(size = 20),
        legend.key.size = unit(1,"cm"))


p14

ggsave(paste(figures.dir, "NF_chao1_root&rhizo.pdf", sep=""), p14, width=15, height=8)













# compute summary statistics

library(dplyr)



idx <- index$Compartment%in%c("rhizosphere")&
        index$Nitrate_supplement%in%c("none")
chao_rhizo_none <- index[idx,]

chao_rhizo_none_summarize <- chao_rhizo_none%>%
   group_by(Genotype)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))


idx <- index$Compartment%in%c("rhizosphere")&
  index$Nitrate_supplement%in%c("10 mM KNO3")
chao_rhizo_N <- index[idx,]

chao_rhizo_N_summarize <- chao_rhizo_N%>%
  group_by(Genotype)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))




idx <- index$Compartment%in%c("root")&
      index$Nitrate_supplement%in% c("none")
chao_root_none <- index[idx,]

chao_root_none_summarize <- chao_root_none%>%
  group_by(Genotype)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))


idx <- index$Compartment%in%c("root")&
  index$Nitrate_supplement%in% c("10 mM KNO3")
chao_root_N <- index[idx,]

chao_root_N_summarize <- chao_root_N%>%
  group_by(Genotype)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))



## SINK file to write all the summarize into one
sink("Chao1_summarize.csv")

cat('chao_rhizo_none_summarize')
write.csv(chao_rhizo_none_summarize)
cat('____________________________')

cat('\n')
cat('\n')

cat('chao_rhizo_N_summarize')
write.csv(chao_rhizo_N_summarize)
cat('____________________________')

cat('\n')
cat('\n')


cat('chao_root_none_summarize')
write.csv(chao_root_none_summarize)
cat('____________________________')

cat('\n')
cat('\n')



cat('chao_root_N_summarize')
write.csv(chao_root_N_summarize)
cat('____________________________')

cat('\n')
cat('\n')


sink()





# statistical analyses by Mann-Whitney U-test


### N vs none in compartment

idx <- index$Compartment%in% c("soil")

index_soil <- index[idx,]

wilcox.test(value ~ Nitrate_supplement, data = index_soil)



chao_rhizo <- rbind(chao_rhizo_N,chao_rhizo_none)

wilcox.test(value ~ Nitrate_supplement, data = chao_rhizo)


chao_root <- rbind(chao_root_N, chao_root_none)
wilcox.test(value ~ Nitrate_supplement, data = chao_root)




library(ggpubr)

p1 <- chao_rhizo %>%
  ggplot(aes(x=Genotype, y= value, fill=Nitrate_supplement))+
  geom_boxplot()+
  stat_compare_means(method = "wilcox.test")

p1

ggsave(paste(figures.dir, "Chao1_rhizo_wilcox_test_value.pdf", sep=""), p1, width=15, height=8)





p2 <- chao_root %>%
  ggplot(aes(x=Genotype, y= value, fill=Nitrate_supplement))+
  geom_boxplot()+
  stat_compare_means(method = "wilcox.test")
p2
ggsave(paste(figures.dir, "Chao1_root_wilcox_test_value.pdf", sep=""), p2, width=15, height=8)




### check the differences one by one

idx <- chao_rhizo_N$Genotype %in% c("Gifu","nfre")
chao_rhizo_N_nfre <- chao_rhizo_N[idx, ]

wilcox.test(value~Genotype, data = chao_rhizo_N_nfre)



idx <- chao_rhizo_N$Genotype %in% c("Gifu","chit5")
chao_rhizo_N_chit5 <- chao_rhizo_N[idx, ]

wilcox.test(value~Genotype, data = chao_rhizo_N_chit5)


idx <- chao_rhizo_N$Genotype %in% c("Gifu","nfr5")
chao_rhizo_N_nfr5 <- chao_rhizo_N[idx, ]

wilcox.test(value~Genotype, data = chao_rhizo_N_nfr5)



idx <- chao_rhizo_none$Genotype %in% c("Gifu","nfre")
chao_rhizo_none_nfre <- chao_rhizo_none[idx, ]

wilcox.test(value~Genotype, data = chao_rhizo_none_nfre)



idx <- chao_rhizo_none$Genotype %in% c("Gifu","chit5")
chao_rhizo_none_chit5 <- chao_rhizo_none[idx, ]

wilcox.test(value~Genotype, data = chao_rhizo_none_chit5)


idx <- chao_rhizo_none$Genotype %in% c("Gifu","nfr5")
chao_rhizo_none_nfr5 <- chao_rhizo_none[idx, ]

wilcox.test(value~Genotype, data = chao_rhizo_none_nfr5)






idx <- chao_root_none$Genotype %in% c("Gifu","nfre")
chao_root_none_nfre <- chao_root_none[idx, ]

wilcox.test(value~Genotype, data = chao_root_none_nfre)



idx <- chao_root_none$Genotype %in% c("Gifu","chit5")
chao_root_none_chit5 <- chao_root_none[idx, ]

wilcox.test(value~Genotype, data = chao_root_none_chit5)


idx <- chao_root_none$Genotype %in% c("Gifu","nfr5")
chao_root_none_nfr5 <- chao_root_none[idx, ]

wilcox.test(value~Genotype, data = chao_root_none_nfr5)



idx <- chao_root_N$Genotype %in% c("Gifu","nfre")
chao_root_N_nfre <- chao_root_N[idx, ]

wilcox.test(value~Genotype, data = chao_root_N_nfre)



idx <- chao_root_N$Genotype %in% c("Gifu","chit5")
chao_root_N_chit5 <- chao_root_N[idx, ]

wilcox.test(value~Genotype, data = chao_root_N_chit5)


idx <- chao_root_N$Genotype %in% c("Gifu","nfr5")
chao_root_N_nfr5 <- chao_root_N[idx, ]

wilcox.test(value~Genotype, data = chao_root_N_nfr5)


