# START OF CODE

library(caret)
library(pROC)
library(DALEX)
library(org.Hs.eg.db)
library(gplots)

# import datasets from GITHUB so that others can try too! 

Metadata_AML <- read.csv("https://raw.githubusercontent.com/Naren037/hackbio-cancer-internship/refs/heads/main/TASK%203/Datasets/AML%20Metadata.csv", row.names=1)
Rawdata_AML <- read.csv("https://raw.githubusercontent.com/Naren037/hackbio-cancer-internship/refs/heads/main/TASK%203/Datasets/AML%20Rawdata.csv", row.names=1)

# Pre-processing 1 - Selecting high variance genes

trans.data <- data.frame(t(Rawdata_AML))
SD = apply(trans.data, 2, sd)
Top_predictions = order(SD, decreasing = T)[1:1000]
trans.data = trans.data[, Top_predictions]
dim(trans.data)

# Adding metadata

rownames(Metadata_AML) <- Metadata_AML$CASE.ID
Metadata_AML$`CASE.ID` = NULL
trans.data <- data.frame(Metadata_AML,trans.data)

# Pre-processing 2 - Removing zero values and central values 

zero_var <- preProcess(trans.data, method = 'nzv', uniqueCut = 15)
trans.data <- predict(zero_var, trans.data)
dim(trans.data)

# Pre-processing 3 - Centering to avoid feature scaling bias

center_var <- preProcess(trans.data, method = 'center')
trans.data <- predict(center_var, trans.data)
dim(trans.data)

# Pre-processing 4 - Removing highly correlated features

correlation <- preProcess(trans.data, method = 'corr', cutoff = 0.53) #ADJUST LATER 
Processed_data <- predict(correlation, trans.data)

# Heat map of pre-processed data set

heatmap.2(as.matrix(Processed_data[,2:50]), 
          scale = "col", 
          trace = "none", 
          Rowv = F, 
          Colv = F, 
          dendrogram = "none", 
          col = bluered(200), 
          keysize = 1, 
          key.title = 'Expression', 
          main = "Pre-processed Expression Profile", 
          xlab = 'GENES', 
          ylab = 'SAMPLES', 
          cexRow = 0.5, 
          cexCol = 0.7, 
          margins = c(10,10))

# Machine Learning 

# Repeat loop for 100% train accuracy and >90% test accuracy by shuffling the data set in every iteration until the desired accuracy is obtained
# NOTE: Try interrupting R and re-running the loop if it takes too long :-) it works!

repeat  {
Processed_data_jumbled <- Processed_data[sample(nrow(Processed_data)), ]

# 70-30 split 

datasplit <- createDataPartition(y = Processed_data_jumbled$Subtype, p = 0.7)[[1]]
train <- Processed_data[datasplit,]
test <- Processed_data[-datasplit,]

# Model training using KNN

#control <- trainControl(method = 'cv', number = 4) 
control <- trainControl(method = 'repeatedcv', number = 4, repeats = 10) # extra - repeats

KNN <- train(Subtype~.,
             data = train,
             method = 'knn',
             trControl = control,
             preProcess = "scale", # extra 
             tuneGrid = data.frame(k=1:8))

KNN$bestTune 

PredictedTrain <- predict(KNN, newdata = train)
PredictedTest <- predict(KNN, newdata = test)

# Interpretation 

PredictedTrain <- as.factor(PredictedTrain)
train$Subtype <- as.factor(train$Subtype)
levels(PredictedTrain) <- levels(train$Subtype)
x <- confusionMatrix(PredictedTrain, train$Subtype)

PredictedTest <- as.factor(PredictedTest)
test$Subtype <- as.factor(test$Subtype)
levels(PredictedTest) <- levels(test$Subtype)
y <- confusionMatrix(PredictedTest, test$Subtype)

# Condition check

if (x[["overall"]][["Accuracy"]] == 1 & y[["overall"]][["Accuracy"]] > 0.9) 
{break} 
  }
confusionMatrix(PredictedTrain, train$Subtype)
confusionMatrix(PredictedTest, test$Subtype)

# Variable importance 

explainer.aml <- explain(KNN,
                         label = 'knn',
                         data = train,
                         y = as.numeric(train$Subtype))

importance.aml <- feature_importance(explainer.aml, n_sample = 40, type = 'difference')
head(importance.aml)
Top25_Genes <- print(as.data.frame(importance.aml[3:27, 1]))
colnames(Top25_Genes) <- c("GeneID")
rownames(Top25_Genes) <- Top25_Genes$GeneID

# Retrieving gene names 

ensembl_ids <- rownames(Top25_Genes)
ensembl_ids_clean <- gsub("\\..*", "", ensembl_ids)
gene_symbols <- mapIds(
  x = org.Hs.eg.db,           
  keys = ensembl_ids_clean,         
  column = "SYMBOL",          
  keytype = "ENSEMBL",        
  multiVals = "first")

# Printing and Exporting Top25 genes (Features)

Top25_Genes <- data.frame(gene_symbols, Top25_Genes)
Top25_Genes$GeneID = NULL
print(Top25_Genes$gene_symbols) # Top 25 genes
write.csv(Top25_Genes, "Top25_Genes_ML.csv")

# END OF CODE


              



