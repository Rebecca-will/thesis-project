library(xgboost)
library(readr)
library(stringr)
library(caret)
library(car)
library(dplyr)
library(mclust)
library(tibble)
library(ggplot2)

?rownames_to_column

colnames(claraTEMTexpressionforRebecca)
#transpose first
#string, need to use fixed = TRUE. gsub is function 
#gsub '.', '-', fixed = TRUE 
#merge


Gene_sig <- t(claraTEMTexpressionforRebecca)
?gsub

Gene_sig <- rownames_to_column(Gene_sig, var = "RowNames")

Gene_sig<- as.data.frame(Gene_sig)
Gene_sig<- as.numeric(Gene_sig)

colnames(Gene_sig) <- Gene_sig[2,]
Gene_sig<- Gene_sig[-c(1,2), ]

cluster_df<- subset(combined_data, select = c(2,24))


names(Gene_sig)[names(Gene_sig) == "GENE"] <- "Sample"
Gene_sig$GENE

Gene_sig$Sample <- gsub('.', '-', Gene_sig$Sample, fixed = TRUE)

gene_sig_final <- merge(Gene_sig,cluster_df, by = 'Sample')

gene_sig_final$cluster

write.csv(gene_sig_final, file = 'gene_sig_final.csv', row.names = FALSE)

Gene_sig$RowNames <- as.character(rownames(Gene_sig))

col_to_convert <- names(gene_sig_final[2:1179])

gene_sig_final <- gene_sig_final %>%
  mutate_at(vars(col_to_convert), as.numeric)

gene_sig_final<- gene_sig_final[, -c(1)]



# load data
part_data_sig <- sample(nrow(gene_sig_final),nrow(gene_sig_final)*0.7)
df_train_sig <- gene_sig_final[part_data_sig,]
df_test_sig <- gene_sig_final[-part_data_sig,]

colnames(df_train)

# Loading labels of train data
labels_sig = df_train_sig['cluster']
df_train1_sig = subset(df_train_sig, select = -c(1179))#removes the cluster
# combine train and test data
#df_all = rbind(df_train,df_test)




param       = list("objective" = "multi:softmax", # multi class classification
                   "num_class"= 3 ,  		# Number of classes in the dependent variable.
                   "eval_metric" = "mlogloss",  	 # evaluation metric 
                   "nthread" = 8,   			 # number of threads to be used 
                   "max_depth" = 16,    		 # maximum depth of tree 
                   "eta" = 0.3,    			 # step size shrinkage 
                   "gamma" = 0,    			 # minimum loss reduction 
                   "subsample" = 0.7,    		 # part of data instances to grow tree 
                   "colsample_bytree" = 1, 		 # subsample ratio of columns when constructing each tree 
                   "min_child_weight" = 12  		 # minimum sum of instance weight needed in a child 
)


predictors_sig <-colnames(df_train1_sig)
predictors_sig
class(predictors_sig)
str(predictors_sig)

label_sig = as.numeric(df_train_sig$cluster)
print(table (label_sig))

label_sig =  as.numeric(df_train_sig$cluster)-1
print(table (label_sig))

set.seed(100)

cv.nround = 200;  # Number of rounds. This can be set to a lower or higher value, if you wish, example: 150 or 250 or 300  
bst.cv_sig = xgb.cv(
  param=param,
  data = as.matrix(df_train1_sig[,predictors_sig]),
  label = label_sig,
  nfold = 3,
  nrounds=cv.nround,
  prediction=TRUE)

print(bst.cv, verbose = TRUE)
cb.evaluation.log()


min.loss.idx_sig = which.min(bst.cv_sig$evaluation_log[, test_mlogloss_mean]) 
#min.loss.idx = which.min(bst.cv$dt[, test.mlogloss.mean]) 
cat ("Minimum logloss occurred in round : ", min.loss.idx_sig, "\n")

summary(bst.cv)


?xgboost
print(bst.cv$evaluation_log[min.loss.idx,])
#print(bst.cv$dt[min.loss.idx,])


bst_sig = xgboost(
  param=param,
  data = as.matrix(df_train1_sig), 
  label = label_sig,
  nrounds=min.loss.idx_sig)


setdiff(names(df_test_sig), names(df_train_sig))
common_predictors_sig <- intersect(predictors_sig, colnames(df_test_sig))

# Make prediction on the testing data.
df_test_sig$prediction = predict(bst_sig, as.matrix(df_test_sig[,predictors_sig]))
?predict

#Translate the prediction to the original class or Species.
df_test_sig$prediction = ifelse(df_test_sig$prediction==0,"1",ifelse(df_test_sig$prediction==1,"2",ifelse(df_test_sig$prediction==2,'3',NA)))

df_test_sig$cluster <- as.factor(df_test_sig$cluster)
df_test_sig$prediction <- as.factor(df_test_sig$prediction)

#Compute the accuracy of predictions.
confusionMatrix( df_test_sig$prediction,df_test_sig$cluster)


importance_matrix_sig = xgb.importance(colnames(bst_sig), model = bst_sig)
importance_matrix_sig

xgb.plot.importance(importance_matrix_sig[1:20,], col = my_colors, main = 'Hypothesis driven model')

summary(gene_sig_final)


set.seed(100)
train_control = trainControl(method="cv",number=10)
model.rf = train(cluster~., data=df_train, trControl=train_control, method="rf")

testing$prediction.rf = predict(model.rf,testing[,predictors])

summary(scaled_df$EMT_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)


intersect(Rokavec, EMT_transition)
intersect(Rokavec, Dry_functional)
intersect(EMT_transition, Dry_functional)


#################################################################################################################################
#################################################################################################################################
#xgboost non hypothesis driven 
#keeping only the data with the highest intesity 



rnaseq_prescreening_intensity <- function(rnaseqtrain, lowlimit=50, fractiontoremove =0.8, isalreadylog2=T, constantaddedwhenlog2performed=1.1){
  
  
  
  print(lowlimit)
  
  print(constantaddedwhenlog2performed)
  
  print(log2(lowlimit +constantaddedwhenlog2performed))
  
  intdat <- apply(rnaseqtrain, 2, function(x) length(which(as.numeric(x)<(log2(lowlimit +constantaddedwhenlog2performed))))/length(as.numeric(x))) ##########this is giving different answers for only some columns if columns such as "Samples" are included in the dataframe, so ensure you pass only genes in rnaseqtrain, you can cbind the samplename and/or end point you want to train against after
  
  
  
  keepintcols <- names(intdat[intdat <fractiontoremove])
  
  rnaseqtrain <- rnaseqtrain[, colnames(rnaseqtrain) %in% keepintcols]
  
  rnaseqtrain
  
}



rnaseq_prescreening_sd <- function(train, keepfeaturesvar=2000){
  
  
  
  ##keep features with the top 2000 standard deviation
  
  sddat <-as.data.frame(apply(train, 2, sd, na.rm=T))
  
  colnames(sddat) <-  "SD"
  
  sddat$FeatureName <- rownames(sddat)
  
  keepvarcols <- sddat[order(sddat$SD, decreasing=T),][1:keepfeaturesvar,]$FeatureName
  
  train <- train[, colnames(train) %in% keepvarcols]
  
  
  
  train
  
}




?rownames_to_column

colnames(claraTEMTexpressionforRebecca)
#transpose first
#string, need to use fixed = TRUE. gsub is function 
#gsub '.', '-', fixed = TRUE 
#merge


Gene_sig_nh <- t(allTCGAexpressionforRebecca)
?gsub

Gene_sig_nh <- rownames_to_column(Gene_sig, var = "RowNames")

Gene_sig_nh <- Gene_sig_nh[,-1]


cluster_df<- subset(combined_data, select = c(2,24))#the dataframe with the clusters 


write.csv(Gene_sig_nh, file = 'Gene_sig_nh.csv', row.names = FALSE)



gene_sig_final <- merge(Gene_sig,cluster_df, by = 'Sample')


col_to_convert_nh <- names(Gene_sig_nh[2:1179])

Gene_sig_nh <- Gene_sig_nh %>%
  mutate(across(all_of(col_to_convert_nh), as.numeric))

#want to remove the sample column to preform the sd/intensity then cbind back in then merge with cluster one.

selection_df<- Gene_sig_nh[,-1]

filtered_genes <- rnaseq_prescreening_intensity(selection_df)
sd_filtered_genes <- rnaseq_prescreening_sd(filtered_genes)

Sample <- Gene_sig_nh[,1]
sample_nh<-cbind(sd_filtered_genes,Sample)

final_nh<- merge(sample_nh, cluster_df, by = 'Sample')

write.csv(final_nh, file = 'Gene_sig_nh.csv', row.names = FALSE)



data_nh<- Gene_sig_nh[,-1]
#cluster=1132


# load data
part_data_nh <- sample(nrow(data_nh),nrow(data_nh)*0.7)
df_train_nh <- data_nh[part_data_nh,]
df_test_nh <- data_nh[-part_data_nh,]

df_train_nh$cluster<- as.factor(df_train_nh$cluster)

colnames(df_train)

# Loading labels of train data
labels_nh = df_train_nh['cluster']
df_train1_nh = subset(df_train_nh, select = -c(1132))#removes the cluster
# combine train and test data
#df_all = rbind(df_train,df_test)




param       = list("objective" = "multi:softmax", # multi class classification
                   "num_class"= 3 ,  		# Number of classes in the dependent variable.
                   "eval_metric" = "mlogloss",  	 # evaluation metric 
                   "nthread" = 8,   			 # number of threads to be used 
                   "max_depth" = 16,    		 # maximum depth of tree 
                   "eta" = 0.3,    			 # step size shrinkage 
                   "gamma" = 0,    			 # minimum loss reduction 
                   "subsample" = 0.7,    		 # part of data instances to grow tree 
                   "colsample_bytree" = 1, 		 # subsample ratio of columns when constructing each tree 
                   "min_child_weight" = 12  		 # minimum sum of instance weight needed in a child 
)


predictors_nh <-colnames(df_train1_nh)
predictors_nh
class(predictors_sig)
str(predictors_sig)

label_nh = as.numeric(df_train_nh$cluster)
print(table (label_nh))

label_nh =  as.numeric(df_train_nh$cluster)-1
print(table (label_nh))

set.seed(100)

cv.nround = 200;  # Number of rounds. This can be set to a lower or higher value, if you wish, example: 150 or 250 or 300  
bst.cv_nh = xgb.cv(
  param=param,
  data = as.matrix(df_train1_nh[,predictors_nh]),
  label = label_nh,
  nfold = 3,
  nrounds=cv.nround,
  prediction=TRUE)

min.loss.idx_nh = which.min(bst.cv_nh$evaluation_log[, test_mlogloss_mean]) 
#min.loss.idx = which.min(bst.cv$dt[, test.mlogloss.mean]) 
cat ("Minimum logloss occurred in round : ", min.loss.idx_nh, "\n")

summary(bst.cv)


?xgboost
print(bst.cv_nh$evaluation_log[min.loss.idx_nh,])
#print(bst.cv$dt[min.loss.idx,])


bst_nh = xgboost(
  param=param,
  data = as.matrix(df_train1_nh), 
  label = label_nh,
  nrounds=min.loss.idx_nh)


#setdiff(names(df_test_sig), names(df_train_sig))
#common_predictors_sig <- intersect(predictors_sig, colnames(df_test_sig))

# Make prediction on the testing data.
df_test_nh$prediction = predict(bst_nh, as.matrix(df_test_nh[,predictors_nh]))
?predict

#Translate the prediction to the original class or Species.
df_test_nh$prediction = ifelse(df_test_nh$prediction==0,"1",ifelse(df_test_nh$prediction==1,"2",ifelse(df_test_nh$prediction==2,'3',NA)))

df_test_nh$cluster <- as.factor(df_test_nh$cluster)
df_test_nh$prediction <- as.factor(df_test_nh$prediction)

#Compute the accuracy of predictions.
confusionMatrix( df_test_nh$prediction,df_test_nh$cluster)


importance_matrix_nh = xgb.importance(colnames(bst_nh), model = bst_nh)
importance_matrix_nh

xgb.plot.importance(importance_matrix_nh[1:20,], col = my_colors, main = 'Non-Hypothesis driven model')


#add the columns to see which it struggles with the most 
#testing df

df_test_nh$cluster_name <- NA

df_test_nh$cluster_name[df_test_nh$prediction==1]<- '1- mesenchymal'
df_test_nh$cluster_name[df_test_nh$prediction==2]<- '2- classical'
df_test_nh$cluster_name[df_test_nh$prediction==3]<- '3- no EMT'

df_test_nh$check <- df_test_nh$cluster==df_test_nh$prediction
summary(df_test_nh$check)


sample_id<- final_nh[1:2]


wrong_class<- df_test_nh[df_test_nh$check == 'FALSE',]

wrong_class_ID <- merge(wrong_class, sample_id, by = 'CFH')


disease_df<- Clara_T3[1:2]

wrong_class_disease <- merge(wrong_class_ID, disease_df, by = 'Sample')
false_disease<- table(wrong_class_disease$Disease, wrong_class_disease$cluster_name, wrong_class_ID$cluster)

wrong_class_disease %>%
  group_by(Disease, cluster_name) %>%
  summarise(n = n())%>%
  arrange(desc(n))%>%
  print(n=50)
  
wrong_class_disease %>%
  select(Disease, cluster_name, cluster)


table(wrong_class_disease)


#cross validate the data using kmeans 

# train_control_nh<- trainControl(method = 'cv',
#                                 number = 3,
#                                 verboseIter = TRUE,
#                                 allowParallel = TRUE)
# 
# grid_tune <- expand.grid(
#   nrounds = c(500,1000,1500),
#   mtry = 5,
#   max_depth = c(2,4,6),
#   eta = 0.3, 
#   gamma = 0, 
#   colsample_bytree = 1, 
#   min_child_weight = 1,
#   subsample = 1
# )
# 
# 
# 
# nh_tune <- train(x= df_train_nh,
#                  y= label_nh,
#                  trControl = train_control_nh,
#                  tuneGrid = grid_tune, 
#                  method = 'rf',
#                  verbose = TRUE)
# 
# 
# 
# 
# ctrl <- trainControl(method = "cv",   # Cross-validation method
#                      number = 10)     # Number of folds
# 
# # Define the training parameters for XGBoost
# params <- list(
#   objective = "multi:softmax",   # Objective function for multi-class classification
#   eval_metric = "mlogloss",      # Evaluation metric
#   num_class = 3,                 # Number of classes
#   nthread = 8                    # Number of threads to use
# )
# 
# 
# tuneGrid <- expand.grid(
#   nrounds = c(10, 50, 100),             # Number of boosting rounds
#   max_depth = c(3, 6, 9),               # Maximum tree depth
#   eta = c(0.01, 0.1, 0.3),              # Learning rate (step size shrinkage)
#   gamma = 0,                            # Minimum loss reduction required to make a further partition on a leaf node
#   colsample_bytree = 1,                 # Subsample ratio of columns when constructing each tree
#   min_child_weight = 1,                 # Minimum sum of instance weight (hessian) needed in a child
#   subsample = 1,                        # Subsample ratio of the training instances
#   num_class = 3                         # Number of classes (for multi-class classification)
# )



# Train the XGBoost classifier using cross-validation
model <- train(x = as.matrix(df_train_nh %>% select(-cluster)),              # Formula for the model, assuming all columns except the target are predictors
               y= factor(df_train_nh$cluster, labels = c("a", "b","c")),             # Data frame containing the data
               method = "xgbTree",      # XGBoost algorithm
               trControl = xgb_trcontrol_1)        # Cross-validation control parameters
               #tuneGrid = xgb_grid_1 )     # Additional XGBoost parameters

#remove the tuneGrid and see if the defaults in the caret package are good.


df_train_nh$cluster<- levels(factor(df_train_nh$cluster, labels = c("a", "b","c")))
#y= as.factor(df_train_nh$cluster)

#below y= df_train_nh[, colnames(df_train_nh) %in% "cluster"]


levels(df_train_nh$cluster) <- make.names(levels(df_train_nh$cluster))

RFE <- rfe(x = as.data.frame(df_train_nh[, !(colnames(df_train_nh) %in% "cluster")]),
           
           y= as.factor(df_train_nh$cluster) ,
           
           #rerank =TRUE,
           
           #saveDetails=TRUE,
           
           sizes = c(seq(5, 19, by=2), seq(20,50, 5), seq(60, 150, 10)),  
           
           rfeControl = rfeControl(functions = caretFuncs,
                                   
                                   method = "repeatedcv",
                                   
                                   number=3,
                                   
                                   repeats = 5,
                                   
                                   verbose = TRUE,
                                   allowParallel = TRUE ),
           
           method='xgbTree',
           
           #ntree=500,
           
           #tuneGrid = xgb_grid_1,
           
           trControl = trainControl(classProbs = TRUE)
           
)  
#try with the smaller tune grid before the xgboost, keep it consistent 

set.seed(100)

temp_df<- df_train_nh[,1:15]
temp_df$cluster <- df_train_nh$cluster
temp_df<- temp_df[sample(rownames(temp_df), 100),]


# Access the cross-validation results
print(model)

?expand.grid



xgb_trcontrol_1 = trainControl(
  method = "cv",
  number = 2,
  verboseIter = TRUE,
  returnData = FALSE,
  returnResamp = "all",
  classProbs = TRUE,
  allowParallel = TRUE
)


xgb_grid_1 = expand.grid(
  nrounds = 3,
  # scale_pos_weight = 0.32, # uncommenting this line leads to the error
  eta = c( 0.2, 0.3),
  max_depth = c(8,16),
  gamma = c(1, 2), 
  subsample = c(0.75),
  min_child_weight = c(3,5), 
  colsample_bytree = c(0.3,0.6,1)
)
#run again with different colsample_bytree
#which factors impact the performance the model 
#understand the model biology 
#leave out the xgb grid 
#0.005, 0.01, 0.05,


# computing the accuracy of the model 

cv_predictions <- predict(model, newdata = df_test_nh)

# Combine the predictions from all folds
combined_predictions <- unlist(cv_predictions)

new_model_df<- cbind(df_test_nh,combined_predictions)

new_model_df$combined_predictions = ifelse(new_model_df$combined_predictions=='a',"1",ifelse(new_model_df$combined_predictions=='b',"2",ifelse(new_model_df$combined_predictions=='c','3',NA)))


new_model_df$combined_predictions <- as.factor(new_model_df$combined_predictions)
new_model_df$cluster <- as.factor(new_model_df$cluster)

confusionMatrix( new_model_df$cluster ,new_model_df$combined_predictions)


importance_matrix_cv = varImp(model)
importance_matrix_cv

imp_mat<- as.data.frame(importance_matrix_cv)

top_features <- importance_matrix_cv$importance[1:20, ]
top_features$Features <- rownames(importance_matrix_cv)

top_features$Feature <- rownames(imp_mat)

# Create a bar plot
barplot(importance_matrix_cv, col = my_colors, main = 'Non-Hypothesis CV driven model',
        names.arg = rownames(top_features), las = 2)

ggplot(
  importance_matrix_cv,
  mapping = NULL,
  top = 20,
  ggtitle('Non-hypothesis CV model')
)

varimpPlo


plot(top_features)

ggplot(top_features, aes(x = rownames(importance), y = Overall)) +
  geom_bar(stat = "identity") +
  labs(x = "Features", y = "Importance", title = "Feature Importance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


##############

RFE <- rfe(x = train[, !(colnames(train) %in% "Endpoint")],
           
           y= train[, colnames(train) %in% "Endpoint"] ,
           
           rerank =TRUE,
           
           saveDetails=TRUE,
           
           sizes = c(seq(5, 19, by=2), seq(20,50, 5), seq(60, 150, 10)),  
           
           rfeControl = rfeControl(functions = caretFuncs,
                                   
                                   method = "repeatedcv",
                                   
                                   number=3,
                                   
                                   repeats = 5,
                                   
                                   verbose = TRUE),
           
           method='xgbTree',
           
           ntree=500,
           
           tuneGrid = expand.grid(mtry=floor(sqrt(ncol(train))-1)),
           
           trControl = trainControl(classProbs = TRUE)
           
)  