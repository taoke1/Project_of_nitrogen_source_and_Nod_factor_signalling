#
#
#

library(edgeR)

### root compartment - Gifu vs nfre

# subset samples

idx <- design$Genotype %in% c("Gifu", "nfre")&
  design$Compartments %in% c("root")&
  design$Nitrate_supplement%in% c("none")
design_subset_nfre <- design[idx, ]
otu_table_subset_nfre <- otu_table[,idx]



### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset_nfre$Genotype

d <- DGEList(counts=otu_table_subset_nfre, group=groups)
d <- calcNormFactors(d)

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)

d2 <- estimateGLMCommonDisp(d, design.mat) ### estimate common dispersion for all the tags
d2 <- estimateGLMTagwiseDisp(d2, design.mat) ### estimate tagwise dispersion

fit <- glmQLFit(d2, design.mat)

lrt_Gifu_nfre <- glmTreat(fit, contrast=c(1, -1), lfc = 1) 



de_Gifu_nfre <- decideTestsDGE(lrt_Gifu_nfre, adjust.method="BY", p.value=0.05) 
Gifu_otus_root <- rownames(otu_table_subset_nfre)[de_Gifu_nfre==1]
nfe_otus_root <- rownames(otu_table_subset_nfre)[de_Gifu_nfre==-1]

nfre_Gifu_pvals <- lrt_Gifu_nfre$table

# write table of enriched OTUs

enriched_otus_root <- data.frame(Gifu_enriched=rownames(otu_table) %in% Gifu_otus_root,
                                 nfre_enriched=rownames(otu_table) %in% nfe_otus_root)
row.names(enriched_otus_root) <- rownames(otu_table_subset_nfre)

idx <- row.names(enriched_otus_root)%in% row.names(otu_table_subset)

enriched_otus_root <- enriched_otus_root[idx,]

write.csv(enriched_otus_root, file="enriched_nfre_root0.3_all.csv")


#

#
#
### root compartment - Gifu vs chit5

# subset samples

idx <- design$Genotype %in% c("Gifu", "chit5")&
  design$Compartments %in% c("root")&
  design$Nitrate_supplement%in% c("none")
design_subset_chit5 <- design[idx, ]
otu_table_subset_chit5 <- otu_table[, idx]




### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset_chit5$Genotype

d <- DGEList(counts=otu_table_subset_chit5, group=groups)
d <- calcNormFactors(d) 

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)

d2 <- estimateGLMCommonDisp(d, design.mat) ### estimate common dispersion for all the tags
d2 <- estimateGLMTagwiseDisp(d2, design.mat) ### estimate tagwise dispersion

fit <- glmQLFit(d2, design.mat)

lrt_Gifu_chit5 <- glmTreat(fit, contrast=c(1, -1), lfc = 1) 



###


de_Gifu_chit5 <- decideTestsDGE(lrt_Gifu_chit5, adjust.method="BY", p.value=0.05) 

Gifu_otus_root <- rownames(otu_table_subset_chit5)[de_Gifu_chit5==1]

Gifu_chit5_pvals <- lrt_Gifu_chit5$table

lrt_chit5_Gifu <- glmTreat(fit, contrast=c(-1, 1), lfc = 1)
de_chit5_Gifu <- decideTestsDGE(lrt_chit5_Gifu, adjust.method='BY', p.value=0.05)

nfe_otus_root <- rownames(otu_table_subset_chit5)[de_chit5_Gifu==1]

chit5_Gifu_pvals <- lrt_chit5_Gifu$table

# write table of enriched OTUs

enriched_otus_root <- data.frame(Gifu_enriched=rownames(otu_table) %in% Gifu_otus_root,
                                 chit5_enriched=rownames(otu_table) %in% nfe_otus_root)
row.names(enriched_otus_root) <- rownames(otu_table_subset_chit5)

idx <- row.names(enriched_otus_root)%in% row.names(otu_table_subset)

enriched_otus_root <- enriched_otus_root[idx,]


write.csv(enriched_otus_root, file="enriched_chit5_root0.3.csv")

#
#
#
### root compartment - Gifu vs nfr5

# subset samples

idx <- design$Genotype %in% c("Gifu", "nfr5")&
  design$Compartments %in% c("root")&
  design$Nitrate_supplement%in% c("none")
design_subset_nfr5 <- design[idx, ]
otu_table_subset_nfr5 <- otu_table[, idx]



### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset_nfr5$Genotype

d <- DGEList(counts=otu_table_subset_nfr5, group=groups)
d <- calcNormFactors(d) 

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)

d2 <- estimateGLMCommonDisp(d, design.mat) ### estimate common dispersion for all the tags
d2 <- estimateGLMTagwiseDisp(d2, design.mat) ### estimate tagwise dispersion

fit <- glmQLFit(d2, design.mat)

lrt_Gifu_nfr5 <- glmTreat(fit, contrast=c(1, -1), lfc = 1) 

### set parameters


de_Gifu_nfr5 <- decideTestsDGE(lrt_Gifu_nfr5, adjust.method="BY", p.value=0.05) 
Gifu_otus_root <- rownames(otu_table_subset_nfr5)[de_Gifu_nfr5==1]

Gifu_nfr5_pvals <- lrt_Gifu_nfr5$table

lrt_nfr5_Gifu <- glmTreat(fit, contrast=c(-1, 1), lfc = 1)
de_nfr5_Gifu <- decideTestsDGE(lrt_nfr5_Gifu, adjust.method='BY', p.value=0.05)

nfe_otus_root <- rownames(otu_table_subset_nfr5)[de_nfr5_Gifu==1]

nfr5_Gifu_pvals <- lrt_nfr5_Gifu$table

# write table of enriched OTUs

enriched_otus_root <- data.frame(Gifu_enriched=rownames(otu_table) %in% Gifu_otus_root,
                                 nfr5_enriched=rownames(otu_table) %in% nfe_otus_root)
row.names(enriched_otus_root) <- rownames(otu_table_subset_nfr5)

idx <- row.names(enriched_otus_root)%in% row.names(otu_table_subset)

enriched_otus_root <- enriched_otus_root[idx,]

write.csv(enriched_otus_root, file="enriched_nfr5_root0.3.csv")


### root compartment - Gifu vs soil

# subset samples

idx <- design$Genotype %in% c("Gifu", "soil")&
  design$Compartments %in% c("root","soil")&
  design$Nitrate_supplement%in% c("none")
design_subset_soil <- design[idx, ]
otu_table_subset_soil <- otu_table[, idx]



### Generalized Linear Model (GLM)

# create DGE list

groups <- design_subset_soil$Genotype

d <- DGEList(counts=otu_table_subset_soil, group=groups)
d <- calcNormFactors(d) 

# fit the GLM

design.mat <- model.matrix(~ 0 + d$samples$group)

d2 <- estimateGLMCommonDisp(d, design.mat) ### estimate common dispersion for all the tags
d2 <- estimateGLMTagwiseDisp(d2, design.mat) ### estimate tagwise dispersion

fit <- glmQLFit(d2, design.mat)

lrt_Gifu_soil <- glmTreat(fit, contrast=c(1, -1), lfc = 1) 


de_Gifu_soil <- decideTestsDGE(lrt_Gifu_soil, adjust.method="BY", p.value=0.05)

Gifu_otus_root <- rownames(otu_table_subset_soil)[de_Gifu_soil==1]

Gifu_soil_pvals <- lrt_Gifu_soil$table

lrt_soil_Gifu <- glmTreat(fit, contrast=c(-1, 1), lfc = 1)
de_soil_Gifu <- decideTestsDGE(lrt_soil_Gifu, adjust.method="BY", p.value=0.05)

nfe_otus_root <- rownames(otu_table_subset_soil)[de_soil_Gifu==1]

soil_Gifu_pvals <- lrt_soil_Gifu$table

# write table of enriched OTUs

enriched_otus_root <- data.frame(Gifu_enriched=rownames(otu_table) %in% Gifu_otus_root,
                                 soil_enriched=rownames(otu_table) %in% nfe_otus_root)
row.names(enriched_otus_root) <- rownames(otu_table_subset_soil)


idx <- row.names(enriched_otus_root)%in% row.names(otu_table_subset)

enriched_otus_root <- enriched_otus_root[idx,]


write.csv(enriched_otus_root, file="enriched_soil_root0.3_all.csv")

