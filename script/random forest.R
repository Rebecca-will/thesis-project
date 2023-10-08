set.seed(100)

library(AppliedPredictiveModeling)
library(DAAG)
library(rpart)
library(rpart.plot)
library(randomForest)
library(caret)
library(ggplot2)
library(lattice)


df_3 <- merge(combined_data, new_df, by = "PATIENT_ID")

colnames(combined_data)
colnames(new_df)

colnames(df_3)


df_3 <- subset(df_3, select = -c(32:37))
df_3 <- subset(df_3, select = -c(36))
df_3 <- subset(df_3, select = -(4:24))
df_3 <- subset(df_3, select = -(2))

TNM_DF <- subset(df_3, select = c(11:13,36))
detach(df_3)
attach(TNM_DF)

TNM_DF$cluster.y <- as.factor(TNM_DF$cluster.y)
class(TNM_DF$cluster.y)


train_2 <- createDataPartition(TNM_DF$cluster.y, p= 0.7, list = FALSE)
trainset_2 <- TNM_DF[train_2,]
validset_2 <- TNM_DF[-train_2,]

trainset_2 <- na.omit(trainset_2)

rf_1 <- randomForest(cluster.y ~ ., data = trainset_2, ntree= 500, importance = TRUE)
importance(rf_1)        
varImpPlot(rf_1, pch=16, col="black")

imp = as.data.frame(importance(rf_1))
imp
imp = cbind(vars=rownames(imp), imp)
imp
imp = imp[order(imp$MeanDecreaseGini),]  
imp$vars = factor(imp$vars, levels=unique(imp$vars))




barplot(imp$MeanDecreaseGini, names.arg=imp$vars,main = 'Importance of the variable using "gini"')

accuracy <- c()
for (i in 1:3) {
  sample_rf <- randomForest(cluster.y ~ ., data = trainset_2, ntree = 100, mtry = i)
  predValid <- predict(sample_rf,validset_2, type = 'class')
  accuracy[i] <- mean(predValid == validset_2$cluster.y)
}
accuracy



table(as.factor(test$diagnosis), predict(object = rf, newdata = test, type = 'class'))




Dec_tree <- rpart(cluster ~ ., data = trainset_2, method = 'class') #is this the right dependent variable to put through?
rpart.plot(Dec_tree)
summary(Dec_tree)
printcp(Dec_tree)
barplot(summary(Dec_tree)$variable.importance, main="Variable Importance Plot ")
plotcp(Dec_tree)
######################


str(TNM_2)

TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE <- as.factor(TNM_2$AJCC_PATHOLOGIC_TUMOR_STAGE)
TNM_2$PATH_T_STAGE <- as.factor(TNM_2$PATH_T_STAGE)
TNM_2$PATH_N_STAGE <- as.factor(TNM_2$PATH_N_STAGE)
TNM_2$PATH_M_STAGE <- as.factor(TNM_2$PATH_M_STAGE)
TNM_2$cluster <- as.factor(TNM_2$cluster)



train_TNM <- createDataPartition(TNM_2$cluster, p= 0.7, list = FALSE)
trainset_TNM <- TNM_2[train_TNM,]
validset_TNM <- TNM_2[-train_TNM,]

trainset_TNM <- na.omit(trainset_TNM)
colSums(is.na(trainset_TNM))

?na.omit

rf_TNM <- randomForest(cluster ~ ., data = trainset_TNM, ntree= 500, importance = TRUE)
importance(rf_TNM)        
varImpPlot(rf_TNM, pch=16, col="black")

imp_TNM = as.data.frame(importance(rf_TNM))
imp
imp_TNM = cbind(vars=rownames(imp_TNM), imp_TNM)
imp
imp_TNM = imp_TNM[order(imp_TNM$MeanDecreaseGini),]  
imp_TNM$vars = factor(imp_TNM$vars, levels=unique(imp_TNM$vars))




barplot(imp_TNM$MeanDecreaseGini, names.arg=imp_TNM$vars,main = 'Importance of the variable using "gini"')



accuracy <- c()
for (i in 1:3) {
  sample_rf <- randomForest(cluster ~ ., data = trainset_TNM, ntree = 100, mtry = i)
  predValid <- predict(sample_rf,validset_TNM, type = 'class')
  accuracy[i] <- mean(predValid == validset_TNM$cluster)
}
accuracy



table(as.factor(test$diagnosis), predict(object = rf, newdata = test, type = 'class'))




Dec_tree <- rpart(cluster ~ ., data = trainset_2, method = 'class') #is this the right dependent variable to put through?
rpart.plot(Dec_tree)
summary(Dec_tree)
printcp(Dec_tree)
barplot(summary(Dec_tree)$variable.importance, main="Variable Importance Plot ")
plotcp(Dec_tree)


