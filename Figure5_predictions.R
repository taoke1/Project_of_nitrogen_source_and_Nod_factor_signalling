# Code by Ib Thorsgaard Jensen
library(data.table)
library(magrittr)
library(compositions)
library(randomForest)
library(e1071)
library(glmnet)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Changes working directory to the path of this file. Only works in Rstudio

comp <- "root" #compartment in which to perform the predictions

# Loading data from [Zgadzaj et. al 2016]
pnas_data <- fread("../../Data/otu_table_PNAS2016.txt")
pnas_design <- fread("../../Data/design_PNAS2016.txt")
# Subsetting data
pnas_design <- pnas_design[genotype %in% c("gifu", "nfr5_2", "nfr5_3") & compartment == comp]
# Constructing the state variable
pnas_design[,State:=fcase(
  KNO3_treatment == "treated", "Inorganic",
  genotype == "gifu" & KNO3_treatment == "untreated", "Symbiotic",
  genotype != "gifu" & KNO3_treatment == "untreated", "Starved"
)]

# Constructing a unique id for every biological replicate
b <- pnas_design$SampleID
if(comp == "root"){
  b <- gsub(".r1", "", b); b <- gsub(".r2", "", b); b <- gsub(".r3", "", b)
  b <- gsub(".R1", "", b); b <- gsub(".R2", "", b); b <- gsub(".R3", "", b)
}
if(comp == "rhizosphere"){
  b <- gsub(".rh1", "", b); b <- gsub(".rh2", "", b); b <- gsub(".rh3", "", b)
  b <- gsub(".Rh1", "", b); b <- gsub(".Rh2", "", b); b <- gsub(".Rh3", "", b)
}
pnas_design[, BioRep:=b]
BioReps <- unique(pnas_design$BioRep)

# Importing and formatting taxonomy-file
taxonomy.init <- RDPutils::import_sintax_file("../../Data/rep_seqs_ESV_db.sintax", confidence = 0.7)
taxonomy.init <- as.data.frame(taxonomy.init)
taxonomy <- sapply(taxonomy.init, function(x) substr(x, 3, 100))
taxonomy <- data.table(OTUid=rownames(taxonomy.init), taxonomy)
taxonomy <- taxonomy[match(pnas_data$OTUid,OTUid)]

# Collapsing data to the genus-level and constructing data.table with both metadata and abundances
genus_pnas <- rowsum(pnas_data[,-1], taxonomy$Genus)
genus_pnas <- data.table(SampleID = colnames(genus_pnas), t(genus_pnas))
genus_pnas <- genus_pnas[SampleID %in% pnas_design$SampleID]
pnas_design <- pnas_design[SampleID %in% genus_pnas$SampleID]
genus_pnas <- genus_pnas[match(pnas_design$SampleID, SampleID),]
full_data <- data.table(State = pnas_design$State, BioRep = pnas_design$BioRep, genus_pnas[,-1])
full_data[,State:=factor(State)]

set.seed(1635491528) # seed is the system-time as of the first time this code was run

# Biological replicates in each state
sym_rep <- pnas_design[State == "Symbiotic", BioRep] %>% unique()
inog_rep <- pnas_design[State == "Inorganic", BioRep] %>% unique()
starv_rep <- pnas_design[State == "Starved", BioRep] %>% unique()

# The number of replicates in each state to include in the test dataset
n_sym <- round(length(sym_rep)/3); n_inog <- round(length(inog_rep)/3); n_starv <- round(length(starv_rep)/3)

# Randomly selecting which biological replicates should be assigned to the test dataset
test_reps <- c( sample(sym_rep,n_sym), sample(inog_rep,n_inog), sample(starv_rep, n_starv) )

# Constructing training and test datasets
train_reps <- BioReps[!(BioReps %in% test_reps)]
test_data <- full_data[BioRep %in% test_reps]
train_data <- full_data[!(BioRep %in% test_reps)]

# Loading initial candidates of predictor genera
if(comp == "root"){
  genus_table <- fread(file = "Results/Nutrtion_associated_bacteria/State_Associated_genera_root.csv")
}
if(comp == "rhizosphere"){
  genus_table <- fread(file = "Results/Nutrtion_associated_bacteria/State_Associated_genera_rhizosphere.csv")
}
genus_pred <- genus_table$Genus

# Filtering and normalizing training data
absent <- apply(train_data[,-c(1,2)], 2, function(x) mean(x == 0)>=0.9)
absent.genus <- names(absent)[absent]
train_data[,(absent.genus):=NULL]
train_data[,3:ncol(train_data) := (as.matrix(train_data[,-c(1,2)])+1) %>% as.data.table()]
train_data[,3:ncol(train_data) := apply(train_data[,-c(1,2)], 1, clr) %>% t() %>% as.data.table()]

# Keeping the predictor genera that are present in the training dataset 
genus_pred <- intersect(genus_pred, colnames(train_data))

# Filtering and normalizing test data
absent <- apply(test_data[,-c(1,2)], 2, function(x) mean(x == 0)>=0.9)
absent.genus <- names(absent)[absent]
absent.genus <- absent.genus[!(absent.genus %in% genus_pred)]
test_data[,(absent.genus):=NULL]
test_data[,3:ncol(test_data) := (as.matrix(test_data[,-c(1,2)])+1) %>% as.data.table()]
test_data[,3:ncol(test_data) := apply(test_data[,-c(1,2)], 1, clr) %>% t() %>% as.data.table()]

# Initial random forest fit for the variable selection step
x <- train_data[,..genus_pred]
r <- randomForest(y = train_data$State,
                  x = x,
                  importance = T)
acc <- mean(r$predicted == train_data$State)

# Identifying the least important variable
imp <- importance(r, type = 1, scale = F)
min_imp_idx <- which.min(imp)
min_imp <- min(imp)
last_omitted <- colnames(x)[min_imp_idx]
acc_diff <- 0

# Iteratively leaving out the least important variable
while(acc_diff>=0){
  acc_prev <- acc
  x[,c(last_omitted):=NULL]
  r <- randomForest(y = train_data$State,
                    x = x,
                    importance = T)
  acc <- mean(r$predicted == train_data$State)
  acc_diff <- acc - acc_prev
  imp <- importance(r, type = 1, scale = F)
  min_imp_idx <- which.min(imp); min_imp <- min(imp)
  last_omitted <- colnames(x)[min_imp_idx]
}
genus_keep <- colnames(x) # Remaining covariates after variable selection

# Saving table of selected predictor variables for later inspection
filename <- c("Predictor","genera",comp) %>% paste(collapse = "_") %>% paste0(".csv")
dirname <- c("Results","Prediction analysis",filename) %>% paste(collapse = "/")
fwrite(genus_table[Genus %in% genus_keep], file = dirname)

# The training data to be used for fitting
x <- train_data[,..genus_keep]
y <- train_data$State

# Tuning and fitting a ridge regression on the training dataset with the selected variables
cc <- cv.glmnet(x = as.matrix(x), y = y, family = "multinomial", type.measure = "class", alpha = 0)
g <- glmnet(x = x, y = y, family = "multinomial", alpha = 0)
newx <- test_data[,..genus_keep]

# Predict on the test dataset
p_glm <- predict(g, as.matrix(newx), type = "class", s = cc$lambda.min)

# Tuning and fitting a support vector machine on the training dataset with the selected variables
ss <- tune.svm(x,y,
               data = train_data,
               cost = 2^(-10:10),
               kernel = "linear")
ss2 <- tune.svm(x, y,
                data = train_data,
                cost = 2^(-10:10),
                gamma = 2^(-10:2),
                kernel = "radial")
# The tuning CV-misclassification is lower for the radial kernel, so we pick that one
s <- svm(x, y,
         data = train_data,
         kernel = "radial",
         gamma = ss2$best.parameters[1],
         cost = ss2$best.parameters[2])

# Predict on the test dataset
p_svm <- predict(s, newx)

# Tuning and fitting a random forest on the training dataset with the selected variables
rr <- tune.randomForest(x, y,
                        data = train_data,
                        mtry = 1:length(genus_keep))
r <- randomForest(x, y,
                  data = train_data,
                  mtry = rr$best.parameters$mtry)

# Predict on the test dataset
p_rf <- predict(r, test_data)

# Evaluate prediction accuracy of the three methods
mean(p_glm == test_data$State)
mean(p_svm == test_data$State)
mean(p_rf == test_data$State)

C <- caret::confusionMatrix(p_svm, test_data$State)
C_svm <- C$table %*% diag(1/colSums(C$table))
colnames(C_svm) <- rownames(C_svm)

C <- caret::confusionMatrix(p_rf, test_data$State)
C_rf <- C$table %*% diag(1/colSums(C$table))
colnames(C_rf) <- rownames(C_rf)

C <- caret::confusionMatrix(factor(p_glm[,1]), test_data$State)
C_glm <- C$table %*% diag(1/colSums(C$table))
colnames(C_glm) <- rownames(C_glm)

C_sup <- cbind(C_rf, C_glm)

fwrite(C_svm, file = paste0("Results/Confusion_matrix_", comp, ".csv"))
fwrite(C_sup, file = paste0("Results/Confusion_matrix_supplementary_", comp, ".csv"))

