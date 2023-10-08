Gene_sig_nh

#bind the diseases by sample ID. 
mes_sample <- subset(No_mes, select = c(2,3))

emt_sig <- merge(Gene_sig_nh, mes_sample, by = 'Sample')
dim(emt_sig)

emt_sig_1<- subset(emt_sig, select = -c(1134,1))
emt_sig_1 <- emt_sig_1[!(emt_sig_1$cluster == '1'),]

# load data
part_data_emt <- sample(nrow(emt_sig_1),nrow(emt_sig_1)*0.7)
df_train_emt <- emt_sig_1[part_data_emt,]
df_test_emt <- emt_sig_1[-part_data_emt,]

part_data_emt_val <- sample(nrow(df_test_emt), nrow(df_test_emt) * 0.5)
df_val_emt <- df_test_emt[part_data_emt_val, ]
df_test_emt_2 <- df_test_emt[-part_data_emt_val, ]



colnames(df_train)
emt_sig$di


# Loading labels of train data
labels_emt = df_train_emt['cluster']



df_train1_emt = subset(df_train_emt, select = -c(1132))#removes the cluster
# combine train and test data
#df_all = rbind(df_train,df_test)




param       = list("objective" = "binary:logistic",  		# Number of classes in the dependent variable.
                   "eval_metric" = "logloss",  	 # evaluation metric 
                   "nthread" = 8,   			 # number of threads to be used 
                   "max_depth" = 16,    		 # maximum depth of tree 
                   "eta" = 0.3,    			 # step size shrinkage 
                   "gamma" = 1,    			 # minimum loss reduction 
                   "subsample" = 0.7,    		 # part of data instances to grow tree 
                   "colsample_bytree" = 1, 		 # subsample ratio of columns when constructing each tree 
                   "min_child_weight" = 3  		 # minimum sum of instance weight needed in a child 
)


predictors_emt<-colnames(df_train1_emt)
predictors_emt
class(predictors_sig)
str(predictors_sig)

label_emt = as.numeric(df_train_emt$cluster)
print(table (label_emt))

label_emt =  as.numeric(df_train_emt$cluster)-2
print(table (label_emt))

set.seed(100)

cv.nround = 100;  # Number of rounds. This can be set to a lower or higher value, if you wish, example: 150 or 250 or 300  
bst.cv_emt = xgb.cv(
  param=param,
  data = as.matrix(df_train1_emt[,predictors_emt]),
  label = label_emt,
  nfold = 3,
  nrounds=cv.nround,
  prediction=TRUE)

print(bst.cv_emt, verbose = TRUE)
cb.evaluation.log()


min.loss.idx_emt = which.min(bst.cv_emt$evaluation_log[, test_logloss_mean]) 
#min.loss.idx = which.min(bst.cv$dt[, test.mlogloss.mean]) 
cat ("Minimum logloss occurred in round : ", min.loss.idx_emt, "\n")




?xgboost
print(bst.cv_emt$evaluation_log[min.loss.idx,])
#print(bst.cv$dt[min.loss.idx,])


bst_emt = xgboost(
  param=param,
  data = as.matrix(df_train1_emt), 
  label = label_emt,
  nrounds=min.loss.idx_emt)


# Make prediction on the testing data.
df_test_emt$prediction_prob = predict(bst_emt, as.matrix(df_test_emt[,predictors_emt]))

threshold <- 0.5

# Convert predicted probabilities to binary labels
df_test_emt$predicted_labels <- ifelse(df_test_emt$prediction_prob >= threshold, 1, 0)

#Translate the prediction to the original class or Species.
df_test_emt$predicted_labels = ifelse(df_test_emt$predicted_labels==0,"2",ifelse(df_test_emt$predicted_labels==1,"3",NA))

df_test_emt$cluster <- as.factor(df_test_emt$cluster)
df_test_emt$predicted_labels <- as.factor(df_test_emt$predicted_labels)

#Compute the accuracy of predictions.
confusionMatrix( df_test_emt$predicted_labels,df_test_emt$cluster)




importance_matrix_emt = xgb.importance(colnames(bst_emt), model = bst_emt)
importance_matrix_emt

xgb.plot.importance(importance_matrix_emt[1:20,], col = my_colors, main = 'EMT driven model')

library(pROC)

label_emt_roc =  as.numeric(df_test_emt$cluster)-1
print(table (label_emt_roc))


roc_obj <- roc(label_emt_roc, df_test_emt$prediction_prob)
auc <- auc(roc_obj)

# Print AUC
print(paste("AUC:", auc))

plot(roc_obj, main = "ROC Curve")
text(0.5, 0.3, paste("AUC =", round(auc, 3)), col = "blue", cex = 1.2)

# Add diagonal reference line
abline(a = 0, b = 1, col = "gray", lty = 2)


xgb_trcontrol_1 = trainControl(
  method = "repeatedcv",
  number = 10,
  verboseIter = TRUE,
  returnData = FALSE,
  returnResamp = "all",
  classProbs = TRUE,
  allowParallel = TRUE
)


xgb_grid_1 = expand.grid(
  nrounds = 10,
  # scale_pos_weight = 0.32, # uncommenting this line leads to the error
  eta = c( 0.5, 0.01, 0.001),
  max_depth = c(10),
  gamma = c(4), 
  subsample = c(0.5,0.75),
  min_child_weight = c(3,5,8,12), 
  colsample_bytree = c(0.3,0.6,1)
)
#run again with different colsample_bytree
#which factors impact the performance the model 
#understand the model biology 
#leave out the xgb grid 
#0.005, 0.01, 0.05,


# computing the accuracy of the model 

df_train_emt$cluster<- as.factor(df_train_emt$cluster)
df_val_emt$cluster<- as.factor(df_val_emt$cluster)

df_val_emt <- subset(df_val_emt, by = -c(1133:1135))

levels(df_val_emt$cluster) <- make.names(levels(df_val_emt$cluster))

RFE_<- rfe(x = as.data.frame(df_val_emt[, !(colnames(df_val_emt) %in% "cluster")]),
           
           y= as.factor(df_val_emt$cluster) ,
           
           #rerank =TRUE,
           
           #saveDetails=TRUE,
           
           sizes =  130,  
           
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

saveRDS(RFE, 'RFE_model.rds')

model_cv <- train(x = as.matrix(df_val_emt %>% select(-cluster)),              # Formula for the model, assuming all columns except the target are predictors
               y= as.factor(df_val_emt$cluster),             # Data frame containing the data
               method = "xgbTree",      # XGBoost algorithm
               trControl = xgb_trcontrol_1,        # Cross-validation control parameters
                tuneGrid = xgb_grid_1 )  




cv_predictions_emt <- predict(RFE_model, newdata = df_test_emt)

# Combine the predictions from all folds
combined_predictions_emt <- unlist(cv_predictions_emt)

new_model_df_emt<- cbind(df_test_emt,combined_predictions_emt)

new_model_df_emt$combined_predictions_emt = ifelse(new_model_df_emt$combined_predictions_emt=='1',"2",ifelse(new_model_df_emt$combined_predictions_emt=='2',"3",NA))


new_model_df_emt$combined_predictions_emt <- as.factor(new_model_df_emt$combined_predictions_emt)
new_model_df_emt$cluster <- as.factor(new_model_df_emt$cluster)

confusionMatrix( new_model_df_emt$cluster ,new_model_df_emt$combined_predictions_emt)


importance_matrix_emt = varImp(RFE_model)
importance_matrix_emt

imp_mat<- as.data.frame(importance_matrix_emt)

top_features_emt <- importance_matrix_emt$importance[1:130, ]
top_features_emt$Features <- rownames(importance_matrix_emt)

top_features_emt$Feature <- rownames(imp_mat)

ggplot(
  importance_matrix_emt,
  mapping = NULL,
  top = 20,
  ggtitle('RFE CV model')
)


importance_df <- as.data.frame(importance_matrix_emt)

# Add rownames as a column
importance_df$variable <- rownames(importance_df)

# Sort the data frame by importance values
importance_df <- importance_df[order(importance_df$Overall, decreasing = TRUE), ]

# Select the top 20 features
top_features <- head(importance_df, 130)

top_features <- top_features %>%
  mutate(variable = factor(variable, levels = variable[order(Overall, decreasing = TRUE)]))


# Create the bar plot
ggplot(top_features, aes(x = Overall, y = variable, fill = variable)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_viridis_d() +
  labs(x = "Importance", y = "Variable", title = "RFE CV Model") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")

sessionInfo()
###########################################################################
###########################################################################
#Final model after CV/RFE

cluster<- as.vector(emt_sig_1$cluster)

opt_variables <- scan("top_130.txt", what = character(), sep = "\n", quiet = TRUE)
#opt_variable <- RFE_model$optvariable

RFE_DF <- emt_sig_1[, opt_variables]
RFE_DF <- cbind(RFE_DF,cluster)
RFE_DF$cluster = ifelse(RFE_DF$cluster=='2-classical',"2",ifelse(RFE_DF$cluster=='3-no emt',"3",NA))




part_data_rfe <- sample(nrow(RFE_DF),nrow(RFE_DF)*0.7)
df_train_rfe <- RFE_DF[part_data_rfe,]
df_test_rfe <- RFE_DF[-part_data_rfe,]


# Loading labels of train data
labels_rfe = df_train_rfe['cluster']

which(colnames(df_train_rfe)== 'cluster')

df_train1_rfe = subset(df_train_rfe, select = -c(131))#removes the cluster
# combine train and test data
#df_all = rbind(df_train,df_test)




param_rfe       = list("objective" = "binary:logistic",  		# Number of classes in the dependent variable.
                   "eval_metric" = "logloss",  	 # evaluation metric 
                   "nthread" = 1,   			 # number of threads to be used 
                   "max_depth" = 10,    		 # maximum depth of tree 
                   "eta" = 0.5,    			 # step size shrinkage 
                   "gamma" = 4,    			 # minimum loss reduction 
                   "subsample" = 0.75,    		 # part of data instances to grow tree 
                   "colsample_bytree" = 0.6, 		 # subsample ratio of columns when constructing each tree 
                   "min_child_weight" = 3  		 # minimum sum of instance weight needed in a child 
)


predictors_rfe<-colnames(df_train1_rfe)
predictors_rfe
class(predictors_sig)
str(predictors_sig)

label_rfe = as.numeric(df_train_rfe$cluster)
print(table (label_rfe))

label_rfe =  as.numeric(df_train_rfe$cluster)-2
print(table (label_rfe))

set.seed(100)

cv.nround = 100;  # Number of rounds. This can be set to a lower or higher value, if you wish, example: 150 or 250 or 300  
bst.cv_rfe = xgb.cv(
  param=param_rfe,
  data = as.matrix(df_train1_rfe[,predictors_rfe]),
  label = label_rfe,
  nfold = 10,
  nrounds=cv.nround,
  prediction=TRUE)

print(bst.cv_rfe, verbose = TRUE)
cb.evaluation.log()


min.loss.idx_rfe = which.min(bst.cv_rfe$evaluation_log[, test_logloss_mean]) 
#min.loss.idx = which.min(bst.cv$dt[, test.mlogloss.mean]) 
cat ("Minimum logloss occurred in round : ", min.loss.idx_rfe, "\n")




?xgboost
print(bst.cv_emt$evaluation_log[min.loss.idx,])
#print(bst.cv$dt[min.loss.idx,])

set.seed(100)
bst_rfe = xgboost(
  param=param_rfe,
  data = as.matrix(df_train1_rfe), 
  label = label_rfe,
  nrounds=min.loss.idx_rfe)


?xgboost


#predicition on the training data 
df_train_rfe$prediction_prob = predict(bst_rfe, as.matrix(df_train_rfe[,predictors_rfe]))
threshold <- 0.5



# Make prediction on the testing data.
df_test_rfe$prediction_prob = predict(bst_rfe, as.matrix(df_test_rfe[,predictors_rfe]))

threshold <- 0.5

# Convert predicted probabilities to binary labels
df_test_rfe$predicted_labels <- ifelse(df_test_rfe$prediction_prob >= threshold, 1, 0)

#Translate the prediction to the original class or Species.
df_test_rfe$predicted_labels = ifelse(df_test_rfe$predicted_labels==0,"2",ifelse(df_test_rfe$predicted_labels==1,"3",NA))

df_test_rfe$cluster <- as.factor(df_test_rfe$cluster)
df_test_rfe$predicted_labels <- as.factor(df_test_rfe$predicted_labels)

#Compute the accuracy of predictions.
z<- confusionMatrix( df_test_rfe$predicted_labels,df_test_rfe$cluster)




importance_matrix_rfe = xgb.importance(colnames(bst_rfe), model = bst_rfe)
importance_matrix_rfe

xgb.plot.importance(importance_matrix_rfe[1:20,], col = my_colors, main = 'RFE driven model')

boosters_list <- bst.cv_rfe$all_booster

#creating the confusion matrix graph
conf_matrix_data <- as.data.frame(z$table)
conf_matrix_data$Prediction <- rownames(conf_matrix_data)

conf_matrix_data$Correct_Prediction <- ifelse(conf_matrix_data$Reference == conf_matrix_data$Prediction, "Correct", "Incorrect")

# Plot the confusion matrix
ggplot(data = conf_matrix_data, aes(x = Prediction, y = Reference)) +
  geom_tile(aes(fill = Correct_Prediction), color = "white") +
  scale_fill_manual(values = c("Correct" = "orange", "Incorrect" = "gray", "green")) +
  theme_minimal() +
  labs(title = "Confusion Matrix",
       x = "Predicted Label",
       y = "True Label")


conf_matrix_data <- data.frame(Reference = c("2", "2", "3", "3"),
                               Prediction = c("2", "3", "2", "3"),
                               Y = c(1302, 75, 58, 722),
                               Correct_Prediction = c("Correct", "Incorrect", "Incorrect", "Correct"))

# Plot the confusion matrix
ggplot(data = conf_matrix_data, aes(x = Prediction, y = Reference)) +
  geom_tile(aes(fill = Correct_Prediction))  +
  theme_minimal() +
  labs(title = "Confusion Matrix",
       x = "Predicted Label",
       y = "True Label")



ggplot(data = conf_matrix_data, aes(x = Prediction, y = Reference)) +
  geom_tile(aes(fill = Y), color = "white") +
  scale_fill_gradient(low = "gray", high = "blue", labels = scales::comma) +
  geom_text(aes(label = Y), vjust = 1) +
  theme_minimal() +
  labs(title = "Confusion Matrix",
       x = "Predicted Label",
       y = "True Label")

# Select the specific booster model you want to plot (e.g., the 5th booster)
selected_booster <- boosters_list[[5]]

# Load the necessary library
library(xgboost)

# Plot the specific tree
xgb.plot.tree(model = bst_rfe, trees = 1)  





library(pROC)

label_rfe_roc =  as.numeric(df_test_rfe$cluster)-1
print(table (label_rfe_roc))


roc_obj_rfe <- roc(label_rfe_roc, df_test_rfe$prediction_prob)
auc_rfe <- auc(roc_obj_rfe)

# Print AUC
print(paste("AUC:", auc_rfe))

plot(roc_obj_rfe, main = "ROC Curve RFE")
text(0.5, 0.3, paste("AUC =", round(auc_rfe, 3)), col = "blue", cex = 1.2)

# Add diagonal reference line
abline(a = 0, b = 1, col = "gray", lty = 2)

#prep for GO 
column_names_total <- colnames(df_train1_emt)
writeLines(column_names_total, "column_names.txt")

writeLines(opt_variables, 'top_130.txt')

#To do after I understand what is meant be reassign the clusters. Remember to do proportion for diseases. 
df_test_rfe$cluster_name <- NA


df_test_rfe$cluster_name[df_test_rfe$predicted_labels==2]<- '2- classical'
df_test_rfe$cluster_name[df_test_rfe$predicted_labels==3]<- '3- no EMT'

df_test_rfe$check <- df_test_rfe$cluster==df_test_rfe$predicted_labels
summary(df_test_rfe$check)


sample_id<- emt_sig[,c(1,504)]


wrong_class<- df_test_rfe[df_test_rfe$check == 'FALSE',]

wrong_class_ID_rfe <- merge(wrong_class, sample_id, by = 'LOXL2')


disease_df<- Clara_T3[1:2]

wrong_class_disease_rfe <- merge(wrong_class_ID_rfe, disease_df, by = 'Sample')
false_disease<- table(wrong_class_disease_rfe$Disease, wrong_class_disease_rfe$cluster_name, wrong_class_ID_rfe$cluster)

wrong_class_disease_rfe %>%
  group_by(Disease, cluster) %>%
  summarise(n = n())%>%
  arrange(desc(n))%>%
  print(n=50)

wrong_class_disease %>%
  select(Disease, cluster_name, cluster)

rfe_disease<- merge(sample_id,disease_df, by = 'Sample')
rfe_disease<- merge(df_test_rfe,rfe_disease, by = 'LOXL2')

wrong_prop<- rfe_disease%>%
  group_by(Disease, check)%>%
  #filter(check %in% 'FALSE')%>%
  summarise(n=n())%>%
  summarise(proportion = n / sum(n), check)%>%
  print(n=40)

write.csv(wrong_prop, file = 'wrong_prop.csv')

#testing on the data with mes 
#df_Test_sig

130/3
  

# Make prediction on the testing data.
df_test_sig$prediction_prob = predict(bst_rfe, as.matrix(df_test_sig[,predictors_rfe]))

threshold <- 0.5

# Convert predicted probabilities to binary labels
df_test_sig$predicted_labels <- ifelse(df_test_sig$prediction_prob >= threshold, 1, 0)

#Translate the prediction to the original class or Species.
df_test_sig$predicted_labels = ifelse(df_test_sig$predicted_labels==0,"2",ifelse(df_test_sig$predicted_labels==1,"3",NA))

df_test_sig$cluster <- as.factor(df_test_sig$cluster)
df_test_sig$predicted_labels <- as.factor(df_test_sig$predicted_labels)

#Compute the accuracy of predictions.
confusionMatrix( df_test_sig$predicted_labels,df_test_sig$cluster)


sig_check<- subset(df_test_sig, cluster == '3')
table(sig_check$predicted_labels)

(766/797)*100



# Survival curve to see if the OS/TNM is maintained
test_rfe<- df_test_rfe[,c(1,133,134)]

sam_ID<- rfe_disease[,c(136,1)]
test_rfe_2<- merge(test_rfe, sam_ID, by= 'LOXL2')

Survival__df<- combined_data[,c(2,25:30)]

surv_rfe<- merge(Survival__df,test_rfe_2, by = 'Sample')

surv_rfe$cluster_colour<- NA
surv_rfe$cluster_colour<-ifelse(surv_rfe$predicted_labels == 3, "#FF7F0E", ifelse(surv_rfe$predicted_labels == 2, '#2CA02C', NA ))

model_fit_surv_rfe<-survfit(Surv(OS_MONTHS,OS_status_binary)~predicted_labels , data = surv_rfe)
model_fit_surv_rfe


surv_rfe$DFS_STATUS[which(surv_rfe$DFS_STATUS =='1:Recurred/Progressed')] <- 'Recurred/Progressed'
surv_rfe$DFS_STATUS[which(surv_rfe$DFS_STATUS =='[Not Available]')] <- NA

table(surv_rfe$DFS_STATUS)


surv_rfe$DFS_binary <- ifelse(surv_rfe$DFS_STATUS == "Recurred/Progressed", 1,
                                   ifelse(surv_rfe$DFS_STATUS == "DiseaseFree", 0,
                                          NA))



surv_rfe_2 <- merge(surv_rfe, disease_df, by = 'Sample')

?survfit.formula
cox_surv_rfe <- coxph(Surv(OS_MONTHS,OS_status_binary)~predicted_labels+ Disease, data = surv_rfe_2)
summary(cox_surv_rfe)
cox_surv

# combined_data<-combined_data$cluster_color[which(cluster=='blue')]<-'#1F77B4'
# combined_data$cluster_color[which(cluster=='green')]<-'#2CA02C'
# combined_data$cluster_color[which(cluster=='orange')]<-'#FF7F0E'

survdiff(Surv(OS_MONTHS,OS_status_binary)~, data = )


# Plot the Kaplan-Meier survival curve with events marked, the main plot 

colours_rfe= unique(surv_rfe[order(surv_rfe$predicted_labels),]$cluster_colour)



surv_plot_rfe<- ggsurvplot(model_fit_surv_rfe, data = surv_rfe, pval = TRUE, pval.method= TRUE,conf.int = TRUE,
                        risk.table = TRUE, risk.table.title = "Risk table",
                        xlab = "Time (in months)", ylab = "Survival Probability",
                        title = "Kaplan-Meier Survival Curve with Events Marked",
                        tables.theme = theme_cleantable(),
                        surv.median.line = "hv",
                        break.time.by = 12,
                        ggtheme = theme_bw(),
                        palette = colours_rfe,
                        ncensor.plot = TRUE,      # is there data for this? 
                        ncensor.plot.height = 0.25)


select_bio<- prob_bio[,c(1,6:9)]
rfe_TNM <- merge(surv_rfe, tnm_final, by = 'Sample')

rfe_TNM$T0<- rfe_TNM$PATH_T_STAGE == "T0"
rfe_TNM$T1<- rfe_TNM$PATH_T_STAGE == "T1"
rfe_TNM$T2<- rfe_TNM$PATH_T_STAGE == "T2"
rfe_TNM$T3<- rfe_TNM$PATH_T_STAGE == "T3"
rfe_TNM$T4<- rfe_TNM$PATH_T_STAGE == "T4"
rfe_TNM$T5<- rfe_TNM$PATH_T_STAGE == "T5"

rfe_TNM$N0<- rfe_TNM$PATH_N_STAGE == "N-"
rfe_TNM$N1<- rfe_TNM$PATH_N_STAGE == "N+"


rfe_TNM$M0<- rfe_TNM$PATH_M_STAGE == "M-"
rfe_TNM$M1<- rfe_TNM$PATH_M_STAGE == "M+"

rfe_TNM$Stage1<- rfe_TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage I"
rfe_TNM$Stage2<- rfe_TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage II"
rfe_TNM$Stage3<- rfe_TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage III"
rfe_TNM$Stage4<- rfe_TNM$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV"


rfe_TNM_2 <- merge(rfe_TNM, disease_df, by = 'Sample')

cox_surv_rfe_tnm_2 <- coxph(Surv(OS_MONTHS, OS_status_binary) ~ predicted_labels+ as.factor(PATH_T_STAGE) + as.factor(PATH_N_STAGE) + as.factor(PATH_M_STAGE) + Disease, data = rfe_TNM_2)
cox_surv_rfe_tnm_3 <- coxph(Surv(OS_MONTHS, OS_status_binary) ~ predicted_labels+ as.factor(PATH_N_STAGE) + as.factor(PATH_M_STAGE), data = rfe_TNM_2)


model_fit_DFS_RFE <- survfit(Surv(DFS_MONTHS,DFS_binary)~predicted_labels, data = surv_rfe)
model_fit_DFS

DFS_plot<- ggsurvplot(model_fit_DFS_RFE, data = surv_rfe, pval = TRUE, conf.int = TRUE,
                      risk.table = TRUE, risk.table.title = "Risk table",
                      xlab = "Time (in months)", ylab = "DFS Probability",
                      title = "Kaplan-Meier DFS Curve with Events Marked",
                      tables.theme = theme_cleantable(),
                      surv.median.line = "hv",
                      break.time.by = 12,
                      ggtheme = theme_bw(),
                      palette = colours_rfe,
                      ncensor.plot = TRUE,      # is there data for this? 
                      ncensor.plot.height = 0.25)


tab_rfe <- rfe_TNM %>%
  filter(cluster_name%in% c("2- classical", "3- no emt")) %>%
  na.omit()%>%
  group_by(PATH_T_STAGE, cluster_name) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(cluster_name) %>%
  mutate(prop = n/sum(n))

rfe_TNM[!(is.na(rfe_TNM$PATH_M_STAGE)),]%>%
  group_by(cluster_name, PATH_M_STAGE) %>%
  summarise(n=n())%>%
  mutate(prop = n/sum(n))%>%
  ggplot(aes(x = PATH_M_STAGE, y = prop, fill = cluster_name)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_rfe) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH M Stage by Cluster") +
  theme_bw()


library(pROC)
library(ggplot2)

label_rfe_train =  as.numeric(df_train_rfe$cluster)-2



# Calculate ROC curve and AUC for training set
train_roc <- roc(label_rfe_train, df_train_rfe$prediction_prob, thresholds = thresholds_common, direction = "auto")
train_auc <- auc(train_roc)

# Calculate ROC curve and AUC for testing set
test_roc <- roc(label_rfe_roc, df_test_rfe$prediction_prob, thresholds = thresholds_common, direction = "auto")
test_auc <- auc(test_roc)

thresholds_common <- seq(0, 1, length.out = 6000)

# Interpolate the ROC curves to the common thresholds
train_roc_interp <- roc(label_rfe_train, df_train_rfe$prediction_prob, threshold= 6000, direction = "auto")
test_roc_interp <- roc(label_rfe_roc, df_test_rfe$prediction_prob, threshold = 6000, direction = "auto")

# Create the roc_data data frame with the same number of rows for both datasets
roc_data <- data.frame(
  FPR = c(train_roc_interp$specificities, test_roc_interp$specificities),
  TPR = c(train_roc_interp$sensitivities, test_roc_interp$sensitivities),
  Dataset = rep(c("Training", "Testing"), each = length(6000))
)

# Plot the ROC curves with different colors for training and testing
ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + # Add diagonal line for reference
  labs(title = "ROC Curve",
       x = "False Positive Rate",
       y = "True Positive Rate") 

ggplot() +
  geom_line(data = coords(train_roc, "all"), aes(x = 1 - specificity, y = sensitivity), color = "blue") +
  geom_line(data = coords(test_roc, "all"), aes(x = 1 - specificity, y = sensitivity), color = "red") +
  labs(title = "ROC Curve Comparison",
       x = "False Positive Rate",
       y = "True Positive Rate",
       color = "Dataset") +
  scale_color_manual(values = c("blue", "red"), labels = c("Training", "Testing")) +
  theme_minimal()



rfe_TNM_2$PATH_T_STAGE[which(rfe_TNM_2$PATH_T_STAGE == 'TX')]<-NA
rfe_TNM_2$PATH_T_STAGE[which(rfe_TNM_2$PATH_T_STAGE == '[Not Available]')]<-NA
rfe_TNM_2$PATH_T_STAGE[which(rfe_TNM_2$PATH_T_STAGE == '[Discrepancy]')]<-NA



rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T1')]<-'T1'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T1a')]<-'T1'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T1a1')]<-'T1'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T1b')]<-'T1'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T1b1')]<-'T1'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T1b2')]<-'T1'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T1c')]<-'T1'


rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE== 'T2')]<-'T2'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE== 'T2a1')]<-'T2'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE== 'T2a2')]<-'T2'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T2b')]<-'T2'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE== 'T2c')]<-'T2'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE== 'T2a')]<-'T2'



rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE== 'T3')]<-'T3'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE== 'T3a')]<-'T3'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T3b')]<-'T3'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T3c')]<-'T3'


rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T4')]<-'T4'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T4a')]<-'T4'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE == 'T4b')]<-'T4'
rfe_TNM$PATH_T_STAGE[which(rfe_TNM$PATH_T_STAGE== 'T4d')]<-'T4'


kruskal.test(rfe_TNM$PATH_N_STAGE, rfe_TNM$predicted_labels)
kruskal.test(rfe_TNM$PATH_M_STAGE, rfe_TNM$predicted_labels)
table(is.na(rfe_TNM$PATH_N_STAGE))
table(is.na(rfe_TNM$PATH_M_STAGE))




result_3<-rfe_TNM_2%>%
  analyse_multivariate(vars(OS_MONTHS, OS_status_binary),
                       vars(predicted_labels, PATH_T_STAGE, PATH_N_STAGE, PATH_M_STAGE, Disease))


mean_OS_MONTHS_paad <- mean(emt_no_emt$OS_MONTHS[emt_no_emt$Disease == "ov"], na.rm = TRUE)
mean_OS_MONTHS_prad <- mean(emt_no_emt$OS_MONTHS[emt_no_emt$Disease == "prad"], na.rm = TRUE)

forest_plot(result_3)
forest_plot(result_3,
            endpoint_labeller = c(time="OS"),
            labels_displayed = c("endpoint", "factor", "n"),
            ggtheme = ggplot2::theme_bw(base_size = 10),
            relative_widths = c(1, 1.5, 1),
            HR_x_breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2))


restable(unique(rfe_TNM_2$Disease))
names(result_3)
