library(xgboost)
library(readr)
library(stringr)
library(caret)
library(car)
library(dplyr)
library(mclust)

set.seed(100)

colnames(Clara_T3)
colnames(sig_df)

sig_df<- subset(Clara_T3,select = c(3:24))


sig_df <- subset(sig_df, select = -c(21))

sig_df$cluster <- ifelse(sig_df$cluster == "#1F77B4", 1,
                                         ifelse(sig_df$cluster == "#2CA02C", 2,
                                                ifelse(sig_df$cluster == '#FF7F0E',3,
                                                NA)))



# class(sig_df$cluster)
# 
# sig_df$cluster <- as.numeric(sig_df$cluster)

scaled_df <- cbind(scaled_data,sig_df[21])



colnames(df_train)

# load data
part_data <- sample(nrow(scaled_df),nrow(scaled_df)*0.7)
df_train <- scaled_df[part_data,]
df_test <- scaled_df[-part_data,]

colnames(df_train)

# Loading labels of train data
labels = df_train['cluster']
df_train1 = subset(df_train, select = -c(21))
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


predictors = colnames(df_train1[-ncol(df_train)])
predictors
class(predictors)
str(predictors)

label = as.numeric(df_train$cluster)
print(table (label))

label =  as.numeric(df_train$cluster)-1
print(table (label))


cv.nround = 200;  # Number of rounds. This can be set to a lower or higher value, if you wish, example: 150 or 250 or 300  
bst.cv = xgb.cv(
  param=param,
  data = as.matrix(df_train1[,predictors]),
  label = label,
  nfold = 3,
  nrounds=cv.nround,
  prediction=TRUE)

print(bst.cv, verbose = TRUE)
cb.evaluation.log()


min.loss.idx = which.min(bst.cv$evaluation_log[, test_mlogloss_mean]) 
#min.loss.idx = which.min(bst.cv$dt[, test.mlogloss.mean]) 
cat ("Minimum logloss occurred in round : ", min.loss.idx, "\n")

summary(bst.cv)


?xgboost
print(bst.cv$evaluation_log[min.loss.idx,])
#print(bst.cv$dt[min.loss.idx,])


bst = xgboost(
  param=param,
  data = as.matrix(df_train1), 
  label = label,
  nrounds=min.loss.idx)


# Make prediction on the testing data.
df_test$prediction = predict(bst, as.matrix(df_test[,predictors]))
?predict

#Translate the prediction to the original class or Species.
df_test$prediction = ifelse(df_test$prediction==0,"1",ifelse(df_test$prediction==1,"2",ifelse(df_test$prediction==2,'3',NA)))

df_test$cluster <- as.factor(df_test$cluster)
df_test$prediction <- as.factor(df_test$prediction)

#Compute the accuracy of predictions.
confusionMatrix( df_test$prediction,df_test$cluster)


importance_matrix = xgb.importance(colnames(bst), model = bst)
importance_matrix

xgb.plot.importance(importance_matrix[1:20,])




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
#xgboost based on the two different clusters. 


colnames(sig_df)


blue_cluster <- subset(scaled_df, select = c(1,17,16,9,7,8,6,14,10,3,4,21))
blue_cluster <- subset(blue_cluster, select = -c(11))
colnames(blue_cluster)



part_data_blue <- sample(nrow(blue_cluster),nrow(blue_cluster)*0.7)
df_train_blue <- blue_cluster[part_data_blue,]
df_test_blue <- blue_cluster[-part_data_blue,]

colnames(df_train)

# Loading labels of train data
labels_blue = df_train_blue['cluster']
df_train1_blue = subset(df_train_blue, select = -c(11))
# combine train and test data
#df_all = rbind(df_train,df_test)




param       = list("objective" = "multi:softmax", # multi class classification
                   "num_class"= 3 ,  		# Number of classes in the dependent variable.
                   "eval_metric" = "mlogloss",  	 # evaluation metric 
                   "nthread" = 8,   			 # number of threads to be used 
                   "max_depth" = 16,    		 # maximum depth of tree 
                   "eta" = 0.5,    			 # step size shrinkage 
                   "gamma" = 0,    			 # minimum loss reduction 
                   "subsample" = 0.7,    		 # part of data instances to grow tree 
                   "colsample_bytree" = 1, 		 # subsample ratio of columns when constructing each tree 
                   "min_child_weight" = 30  		 # minimum sum of instance weight needed in a child 
)


predictors_blue = colnames(df_train1_blue[-ncol(df_train_blue)])
predictors_blue
class(predictors)
str(predictors)

label_blue = as.numeric(df_train_blue$cluster)
print(table (label_blue))

label_blue =  as.numeric(df_train_blue$cluster)-1
print(table (label_blue))


cv.nround = 500;  # Number of rounds. This can be set to a lower or higher value, if you wish, example: 150 or 250 or 300  
bst.cv_blue = xgb.cv(
  param=param,
  data = as.matrix(df_train1_blue[,predictors_blue]),
  label = label,
  nfold = 5,
  nrounds=cv.nround,
  prediction=TRUE)

print(bst.cv_blue, verbose = TRUE)
cb.evaluation.log()


min.loss.idx_blue = which.min(bst.cv_blue$evaluation_log[, test_mlogloss_mean]) 
#min.loss.idx = which.min(bst.cv$dt[, test.mlogloss.mean]) 
cat ("Minimum logloss occurred in round : ", min.loss.idx_blue, "\n")





summary(bst.cv)



print(bst.cv_blue$evaluation_log[min.loss.idx_blue,])
#print(bst.cv$dt[min.loss.idx,])


bst_blue = xgboost(
  param=param,
  data = as.matrix(df_train1_blue), 
  label = label_blue,
  nrounds=min.loss.idx_blue)



# Make prediction on the testing data.
df_test_blue$prediction = predict(bst_blue, as.matrix(df_test_blue[,predictors_blue]))

#Translate the prediction to the original class or Species.
df_test_blue$prediction = ifelse(df_test_blue$prediction==0,"1",ifelse(df_test_blue$prediction==1,"2",ifelse(df_test_blue$prediction==2,'3',NA)))

df_test_blue$cluster <- as.factor(df_test_blue$cluster)
df_test_blue$prediction <- as.factor(df_test_blue$prediction)

#Compute the accuracy of predictions.
confusionMatrix( df_test_blue$prediction,df_test_blue$cluster)


importance_matrix_blue = xgb.importance(colnames(bst_blue), model = bst_blue)
importance_matrix_blue

xgb.plot.importance(importance_matrix_blue[1:10,])



summary_stats_blue <- aggregate(blue_cluster$EMT_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, by = list(blue_cluster$cluster), FUN = function(x) c(mean = mean(x), median = median(x), min = min(x), max = max(x)))
summary_stats_blue





#################################################################################################################################
#################################################################################################################################

orange_cluster <- subset(scaled_df, select = c(13,12,11,2,18,19,15,4,20,5,21))

colnames(orange_cluster)



part_data_orange <- sample(nrow(orange_cluster),nrow(orange_cluster)*0.7)
df_train_orange <- orange_cluster[part_data_orange,]
df_test_orange <- orange_cluster[-part_data_orange,]

colnames(df_train_orange)

# Loading labels of train data
labels_orange = df_train_orange['cluster']
df_train1_orange = subset(df_train_orange, select = -c(11))
# combine train and test data
#df_all = rbind(df_train,df_test)




param       = list("objective" = "multi:softmax", # multi class classification
                   "num_class"= 3 ,  		# Number of classes in the dependent variable.
                   "eval_metric" = "mlogloss",  	 # evaluation metric 
                   "nthread" = 8,   			 # number of threads to be used 
                   "max_depth" = 16,    		 # maximum depth of tree 
                   "eta" = 0.5,    			 # step size shrinkage 
                   "gamma" = 0,    			 # minimum loss reduction 
                   "subsample" = 0.7,    		 # part of data instances to grow tree 
                   "colsample_bytree" = 1, 		 # subsample ratio of columns when constructing each tree 
                   "min_child_weight" = 30  		 # minimum sum of instance weight needed in a child 
)


predictors_orange= colnames(df_train1_orange[-ncol(df_train_orange)])
predictors_orange
class(predictors)
str(predictors)

label_orange = as.numeric(df_train_orange$cluster)
print(table (label_orange))

label_orange =  as.numeric(df_train_orange$cluster)-1
print(table (label_orange))


cv.nround = 500;  # Number of rounds. This can be set to a lower or higher value, if you wish, example: 150 or 250 or 300  
bst.cv_orange = xgb.cv(
  param=param,
  data = as.matrix(df_train1_orange[,predictors_orange]),
  label = label_orange,
  nfold = 5,
  nrounds=cv.nround,
  prediction=TRUE)

print(bst.cv_orange, verbose = TRUE)
cb.evaluation.log()


min.loss.idx_orange = which.min(bst.cv_orange$evaluation_log[, test_mlogloss_mean]) 
#min.loss.idx = which.min(bst.cv$dt[, test.mlogloss.mean]) 
cat ("Minimum logloss occurred in round : ", min.loss.idx_orange, "\n")




summary(bst.cv)



print(bst.cv_orange$evaluation_log[min.loss.idx_orange,])
#print(bst.cv$dt[min.loss.idx,])


bst_orange = xgboost(
  param=param,
  data = as.matrix(df_train1_orange), 
  label = label_orange,
  nrounds=min.loss.idx_orange)




# Make prediction on the testing data.
df_test_orange$prediction = predict(bst_orange, as.matrix(df_test_orange[,predictors_orange]))

#Translate the prediction to the original class or Species.
df_test_orange$prediction = ifelse(df_test_orange$prediction==0,"1",ifelse(df_test_orange$prediction==1,"2",ifelse(df_test_orange$prediction==2,'3',NA)))

df_test_orange$cluster <- as.factor(df_test_orange$cluster)
df_test_orange$prediction <- as.factor(df_test_orange$prediction)

#Compute the accuracy of predictions.
confusionMatrix( df_test_orange$prediction,df_test_orange$cluster)


importance_matrix_orange = xgb.importance(colnames(bst_orange), model = bst_orange)
importance_matrix_blue

xgb.plot.importance(importance_matrix_orange[1:10,])


summary_stats_orange <- aggregate(orange_cluster$EMT_Dry_et_al_MEK_Functional, by = list(orange_cluster$cluster), FUN = function(x) c(mean = mean(x), median = median(x), min = min(x), max = max(x)))
summary_stats_orange




##########################
##########################

importance_matrix <- xgb.importance(model = bst_orange)

# Print the feature importance scores
print(importance_matrix)
\

xgb.plot.importance(importance_matrix = importance_matrix)
xgb.plot.tree(model = bst_orange)


# Access specific feature importance scores
feature_names <- colnames(df_train_orange)
importance_scores <- importance_matrix$Feature
# Access importance scores by feature name
for (feature in feature_names) {
  importance_score <- importance_scores[feature]
  print(paste("Feature:", feature, "Importance:", importance_score))
}



#################################################################################################################################
#################################################################################################################################
#xgboost for t stages

TNM <- merge(Total_path, Clara_T3, by = 'PATIENT_ID')

colnames(TNM)

TNM <- subset(TNM, select = c(2,3,4,5,28))
unique(TNM)

TNM$cluster <- ifelse(TNM$cluster == "#1F77B4", 1,
                         ifelse(TNM$cluster == "#2CA02C", 2,
                                ifelse(TNM$cluster == '#FF7F0E',3,
                                       NA)))





TNM$PATH_N_STAGE <- ifelse(TNM$PATH_N_STAGE == 'N0', 1,
                              ifelse(TNM$PATH_N_STAGE == 'N0 (i-)', 1,
                                     ifelse(TNM$PATH_N_STAGE == 'N0 (i+)',1,
                                            ifelse(TNM$PATH_N_STAGE == 'N0 (mol+)',1, 
                                                   ifelse(TNM$PATH_N_STAGE == 'N1',2,
                                                          ifelse(TNM$PATH_N_STAGE == 'N1a',2,
                                                                 ifelse(TNM$PATH_N_STAGE == 'N1b',2,
                                                                        ifelse(TNM$PATH_N_STAGE == 'N1c',2,
                                                                               ifelse(TNM$PATH_N_STAGE == 'N1mi',2,
                                                                                      ifelse(TNM$PATH_N_STAGE == 'N2', 3,
                                                                                             ifelse(TNM$PATH_N_STAGE == 'N2a',3,
                                                                                                    ifelse(TNM$PATH_N_STAGE == 'N2b',3,
                                                                                                           ifelse(TNM$PATH_N_STAGE== 'N2c',3,
                                                                                                                  ifelse(TNM$PATH_N_STAGE == 'N3', 4,
                                                                                                                         ifelse(TNM$PATH_N_STAGE == 'N3a',4,
                                                                                                                                ifelse(TNM$PATH_N_STAGE == 'N3b',4,
                                                                                                                                       ifelse(TNM$PATH_N_STAGE == 'N3c',4,
                                                                                                                                              ifelse(TNM$PATH_N_STAGE == 'NX',5,
                                                                                                                                                     NA))))))))))))))))))

table(new_df$PATH_N_STAGE)

3724+154+28+1
1116+297+168+5+37


table(TNM$PATH_M_STAGE)

TNM$PATH_M_STAGE <- ifelse(TNM$PATH_M_STAGE == 'cM0 (i+)',1,
                              ifelse(TNM$PATH_M_STAGE == 'M0',1,
                                     ifelse(TNM$PATH_M_STAGE == 'M1',2,
                                            ifelse(TNM$PATH_M_STAGE == 'M1a',2,
                                                   ifelse(TNM$PATH_M_STAGE == 'M1b',2,
                                                          ifelse(TNM$PATH_M_STAGE == 'M1c',2,
 
                                                                 
                                                                                                                                        NA))))))

table(TNM$PATH_T_STAGE)
TNM$PATH_T_STAGE <- ifelse(TNM$PATH_T_STAGE == 'T0',1,
                           ifelse(TNM$PATH_T_STAGE == 'T1',2,
                                  ifelse(TNM$PATH_T_STAGE == 'T1a',2,
                                         ifelse(TNM$PATH_T_STAGE == 'T1a1',2,
                                                ifelse(TNM$PATH_T_STAGE == 'T1b',2,
                                                       ifelse(TNM$PATH_T_STAGE == 'T1b1',2,
                                                              ifelse(TNM$PATH_T_STAGE == 'T1b2',2,
                                                                     ifelse(TNM$PATH_T_STAGE == 'T1c',2,
                                                                            ifelse(TNM$PATH_T_STAGE == 'T2',3,
                                                                                   ifelse(TNM$PATH_T_STAGE == 'T2a',3,
                                                                                          ifelse(TNM$PATH_T_STAGE == 'T2a1',3,
                                                                                                 ifelse(TNM$PATH_T_STAGE == 'T2a2',3,
                                                                                                        ifelse(TNM$PATH_T_STAGE == 'T2b',3,
                                                                                                               ifelse(TNM$PATH_T_STAGE == 'T2c',3,
                                                                                                                      ifelse(TNM$PATH_T_STAGE == 'T3',4,
                                                                                                                             ifelse(TNM$PATH_T_STAGE == 'T3a',4,
                                                                                                                                    ifelse(TNM$PATH_T_STAGE == 'T3a',4,
                                                                                                                                           ifelse(TNM$PATH_T_STAGE == 'T3b',4,
                                                                                                                                                  ifelse(TNM$PATH_T_STAGE == 'T3c',4,
                                                                                                                                                         ifelse(TNM$PATH_T_STAGE == 'T4',5,
                                                                                                                                                                ifelse(TNM$PATH_T_STAGE == 'T4a',5,
                                                                                                                                                                       ifelse(TNM$PATH_T_STAGE == 'T4b',5,
                                                                                                                                                                              ifelse(TNM$PATH_T_STAGE == 'T4d',5,
                                                                                                                                                                                     ifelse(TNM$PATH_T_STAGE == 'Tis',6,
                                                                                                                                                                                            NA))))))))))))))))))))))))
table(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE)


TNM$AJCC_PATHOLOGIC_TUMOR_STAGE <- ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'I/II NOS',1,
                                          ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'IS',1,
                                                 ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage 0',2,
                                                        ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage I',2,
                                                               ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IA',2,
                                                                      ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IA1',2,
                                                                             ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IA2',2,
                                                                                    ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IB',3,
                                                                                           ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IB1',3,
                                                                                                  ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IB2',3,
                                                                                                         ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IC',4,
                                                                                                                ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage II',5,
                                                                                                                       ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE II',5,
                                                                                                                              ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIA',5,
                                                                                                                                     ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIA1',5,
                                                                                                                                            ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIA2',5,
                                                                                                                                                   ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIB',6,
                                                                                                                                                          ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIC',6,
                                                                                                                                                                 ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage III',7,
                                                                                                                                                                        ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE III',7,
                                                                                                                                                                               ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIIA',7,
                                                                                                                                                                                      ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIIB',8,
                                                                                                                                                                                             ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIIC',9,
                                                                                                                                                                                                    ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIIC1',9,
                                                                                                                                                                                                           ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IIIC2',9,
                                                                                                                                                                                                                  ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IV',10,
                                                                                                                                                                                                                         ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE IV',10,
                                                                                                                                                                                                                                ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IVA',10,
                                                                                                                                                                                                                                       ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IVB',11,
                                                                                                                                                                                                                                              ifelse(TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IVC',12,
                                                                                                                                                                                                                                                     NA))))))))))))))))))))))))))))))
                                                                                                         

part_data_TNM <- sample(nrow(TNM),nrow(TNM)*0.7)
df_train_TNM <- TNM[part_data_TNM,]
df_test_TNM <- TNM[-part_data_TNM,]

colnames(df_train_TNM)

# Loading labels of train data
labels_TNM = df_train_TNM['cluster']
df_train1_TNM = subset(df_train_TNM, select = -c(5))
# combine train and test data
#df_all = rbind(df_train,df_test)




param       = list("objective" = "multi:softmax", # multi class classification
                   "num_class"= 3 ,  		# Number of classes in the dependent variable.
                   "eval_metric" = "mlogloss",  	 # evaluation metric 
                   "nthread" = 8,   			 # number of threads to be used 
                   "max_depth" = 16,    		 # maximum depth of tree 
                   "eta" = 0.5,    			 # step size shrinkage 
                   "gamma" = 0,    			 # minimum loss reduction 
                   "subsample" = 0.7,    		 # part of data instances to grow tree 
                   "colsample_bytree" = 1, 		 # subsample ratio of columns when constructing each tree 
                   "min_child_weight" = 30  		 # minimum sum of instance weight needed in a child 
)


predictors_TNM= colnames(df_train1_TNM[-ncol(df_train_TNM)])
predictors_TNM
class(predictors)
str(predictors)

label_TNM = as.numeric(df_train_TNM$cluster)
print(table (label_TNM))

label_TNM =  as.numeric(df_train_TNM$cluster)-1
print(table (label_TNM))


cv.nround = 500;  # Number of rounds. This can be set to a lower or higher value, if you wish, example: 150 or 250 or 300  
bst.cv_TNM = xgb.cv(
  param=param,
  data = as.matrix(df_train1_TNM[,predictors_TNM]),
  label = label_TNM,
  nfold = 5,
  nrounds=cv.nround,
  prediction=TRUE)


print(bst.cv_orange, verbose = TRUE)
cb.evaluation.log()


min.loss.idx_TNM = which.min(bst.cv_TNM$evaluation_log[, test_mlogloss_mean]) 
#min.loss.idx = which.min(bst.cv$dt[, test.mlogloss.mean]) 
cat ("Minimum logloss occurred in round : ", min.loss.idx_TNM, "\n")





summary(bst.cv)#if it fits well on the test data then it has not overfit to the training data



print(bst.cv_orange$evaluation_log[min.loss.idx_orange,])
#print(bst.cv$dt[min.loss.idx,])


bst_TNM = xgboost(
  param=param,
  data = as.matrix(df_train1_TNM), 
  label = label_TNM,
  nrounds=min.loss.idx_TNM)




# Make prediction on the testing data.
df_test_TNM$prediction = predict(bst_TNM, as.matrix(df_test_TNM[,predictors_TNM]))

#Translate the prediction to the original class or Species.
df_test_TNM$prediction = ifelse(df_test_TNM$prediction==0,"1",ifelse(df_test_TNM$prediction==1,"2",ifelse(df_test_TNM$prediction==2,'3',NA)))

df_test_TNM$cluster <- as.factor(df_test_TNM$cluster)
df_test_TNM$prediction <- as.factor(df_test_TNM$prediction)

#Compute the accuracy of predictions.
confusionMatrix( df_test_TNM$prediction,df_test_TNM$cluster)


importance_matrix_TNM = xgb.importance(colnames(bst_TNM), model = bst_TNM)
importance_matrix_TNM

xgb.plot.importance(importance_matrix_TNM[1:4,])



# mc = Mclust(sig_df)
# dr = MclustDR(object = mc, lambda = 1, normalized = T)
# plot(dr, main = "MclustDR Plot for Data", legend = TRUE) 
# plot(dr)
# 
# mylabels<-rownames(dr$x)
# boxed.labels(x,y,mylabels)
# 
# ?plot
# 
# library(plotrix)
#################################################################################################################################
#################################################################################################################################
#enrichment of TNM.

str(TNM_2)

TNM_2 <- merge( Total_path, Clara_T3, by = 'PATIENT_ID')

colnames(TNM_2)

TNM_2 <- subset(TNM_2, select = c(1,28,7,2,3,4,5))
TNM_2$PATH_N_STAGE<- as.factor(TNM_2$PATH_N_STAGE)



sum(is.na(TNM$PATH_N_STAGE))



unique(TNM$PATH_N_STAGE)

TNM_2$cluster <- ifelse(TNM_2$cluster == "#1F77B4", '1-gbm/lgg',
                      ifelse(TNM_2$cluster == "#2CA02C", '2- classical',
                             ifelse(TNM_2$cluster == '#FF7F0E','3- no emt',
                                    NA)))

TNM_2$cluster <- as.factor(TNM_2$cluster)
TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE <-


TNM %>%
  count(PATH_T_STAGE, cluster) %>%
  mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = PATH_T_STAGE, y = n, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#1F77B4", "#2CA02C", '#FF7F0E')) +
  labs(x = "PATH_T_STAGE", y = "Frequency", fill = "Cluster") +
  ggtitle("Frequency of PATH_T_STAGE by Cluster") +
  theme_bw()
  
TNM %>%
  count(PATH_N_STAGE, cluster) %>%
  mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = PATH_N_STAGE, y = n, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#1F77B4", "#2CA02C", '#FF7F0E')) +
  labs(x = "PATH_N_STAGE", y = "Frequency", fill = "Cluster") +
  ggtitle("Frequency of PATH_N_STAGE by Cluster") +
  theme_bw()

TNM %>%
  count(PATH_M_STAGE, cluster) %>%
  mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = PATH_M_STAGE, y = n, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#1F77B4", "#2CA02C", '#FF7F0E')) +
  labs(x = "PATH_M_STAGE", y = "Frequency", fill = "Cluster") +
  ggtitle("Frequency of PATH_M_STAGE by Cluster") +
  theme_bw()

TNM %>%
  count(AJCC_PATHOLOGIC_TUMOR_STAGE, cluster) %>%
  mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = AJCC_PATHOLOGIC_TUMOR_STAGE, y = n, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#1F77B4", "#2CA02C", '#FF7F0E')) +
  labs(x = "Stage", y = "Frequency", fill = "Cluster") +
  ggtitle("Frequency of Stage by Cluster") +
  theme_bw()

#########
TNM_2$PATH_N_STAGE[TNM_2$PATH_N_STAGE == 'NX'] <- NA
TNM_2$PATH_N_STAGE[TNM_2$PATH_N_STAGE == '[Not Available]'] <- NA
TNM_2$PATH_N_STAGE[TNM_2$PATH_N_STAGE == '[Discrepancy]'] <- NA


frequency_data <- TNM_2 %>%
  count(PATH_N_STAGE, cluster) %>%
  group_by(cluster) %>%
  mutate(proportion = n / sum(n))

TNM_2$cluster <- as.factor(TNM_2$cluster)
TNM_2$PATH_N_STAGE <- as.factor(TNM_2$PATH_N_STAGE)


252/(252+1906+1319)
dim(TNM_2)




# Plot the proportions
ggplot(frequency_data, aes(x = PATH_N_STAGE, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#1F77B4", "#2CA02C", "#FF7F0E")) +
  labs(x = "PATH_N_STAGE", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH_N_STAGE by Cluster") +
  theme_bw()


TNM_2$PATH_T_STAGE <- as.factor(TNM_2$PATH_T_STAGE)

table(unique(TNM_2$PATH_T_STAGE))

TNM_2$PATH_T_STAGE[TNM_2$PATH_T_STAGE == 'TX'] <- NA
TNM_2$PATH_T_STAGE[TNM_2$PATH_T_STAGE == '[Not Available]'] <- NA
TNM_2$PATH_T_STAGE[TNM_2$PATH_T_STAGE == '[Discrepancy]'] <- NA


frequency_data_T <- TNM_2 %>%
  count(PATH_T_STAGE, cluster) %>%
  group_by(cluster) %>%
  mutate(proportion = n / sum(n))

ggplot(frequency_data_T, aes(x = PATH_T_STAGE, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#1F77B4", "#2CA02C", "#FF7F0E")) +
  labs(x = "PATH_T_STAGE", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH_T_STAGE by Cluster") +
  theme_bw()


TNM_2$PATH_M_STAGE[TNM_2$PATH_M_STAGE == 'MX'] <- NA
TNM_2$PATH_M_STAGE[TNM_2$PATH_M_STAGE == '[Not Available]'] <- NA
TNM_2$PATH_M_STAGE[TNM_2$PATH_M_STAGE == '[Discrepancy]'] <- NA


TNM_2$PATH_M_STAGE <- as.factor(TNM_2$PATH_M_STAGE)

frequency_data_M <- TNM_2 %>%
  count(PATH_M_STAGE, cluster) %>%
  group_by(cluster) %>%
  mutate(proportion = n / sum(n))

ggplot(frequency_data_M, aes(x = PATH_M_STAGE, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#1F77B4", "#2CA02C", "#FF7F0E")) +
  labs(x = "PATH_M_STAGE", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH_M_STAGE by Cluster") +
  theme_bw()

counts_M <- TNM_2 %>%
  group_by(cluster, PATH_M_STAGE) %>%
  count()
library(tidyverse)
counts_M <- as.data.frame(counts_M)

counts_M <- counts_M %>%
  pivot_wider(names_from = cluster, values_from = n)

counts_M<- na.omit(counts_M)

counts_M <- subset(counts_M,select = c(2:4))
table(counts_M$`1-gbm/lgg`)


counts_M$PATH_M_STAGE[counts_M$PATH_M_STAGE == 'M1b'] <- 2


colnames(counts_M)

result_M <- kruskal.test(counts_M)
############

table(unique(TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE))
TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE[TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage X'] <- NA
TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE[TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE == '[Not Available]'] <- NA
TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE[TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE == '[Discrepancy]'] <- NA
TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE[TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage I'] <- 'STAGE I'
TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE[TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage II'] <- 'STAGE II'
TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE[TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage III'] <- 'STAGE III'
TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE[TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IV'] <- 'STAGE IV'


TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE <- as.factor(TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE)

frequency_data_AJCC <- TNM_2 %>%
  count(AJCC_PATHOLOGIC_TUMOR_STAGE, cluster) %>%
  group_by(cluster) %>%
  mutate(proportion = n / sum(n))

ggplot(frequency_data_AJCC, aes(x = AJCC_PATHOLOGIC_TUMOR_STAGE, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#1F77B4", "#2CA02C", "#FF7F0E")) +
  labs(x = "AJCC_STAGE", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of AJCC_STAGE by Cluster") +
  theme_bw()


table(unique(TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE))


TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE[TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage IV'] <- 'STAGE IV'

?replace





table(unique(TNM_2$PATH_N_STAGE))

TNM_2$PATH_N_STAGE[TNM_2$PATH_N_STAGE == 'NX'] <- NA

N3_DF <- subset(TNM_2, PATH_N_STAGE == 'N3')%>%
  subset(cluster == '1')


T0_DF <- subset(TNM_2, PATH_T_STAGE == 'T0')%>%
  subset(cluster == '1')

Non_app_DF <- subset(TNM_2, PATH_M_STAGE == '[Not Applicable]')





table(TNM_2$cluster, TNM_2$PATH_T_STAGE)

kruskal.test(TNM_2$cluster, TNM_2$PATH_T_STAGE)
kruskal.test(TNM_2$cluster, TNM_2$PATH_N_STAGE)
kruskal.test(TNM_2$cluster, TNM_2$PATH_M_STAGE)
kruskal.test(TNM_2$cluster, TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE)
      













