cluster_analysis_end <- cluster_2_analysis[!(cluster_2_analysis$cluster == '1'),]

model_fit_final<-survfit(Surv(OS_MONTHS,OS_status_binary)~cluster, data = cluster_analysis_end)

cluster_analysis_end$cluster <- ifelse(cluster_analysis_end$cluster == '3','3-no emt',
                                     ifelse(cluster_analysis_end$cluster == "2", '2-classical',
                                                   NA))

surv_plot_end<- ggsurvplot(model_fit_final, data = cluster_analysis_end, pval = TRUE, pval.method= TRUE,conf.int = TRUE,
                        risk.table = TRUE, risk.table.title = "Risk table",
                        xlab = "Time (in months)", ylab = "Survival Probability",
                        title = "Kaplan-Meier Survival Curve with Events Marked",
                        tables.theme = theme_cleantable(),
                        surv.median.line = "hv",
                        break.time.by = 12,
                        ggtheme = theme_bw(),
                        palette = c(  "#2CA02C", '#FF7F0E'),
                        ncensor.plot = TRUE,      # is there data for this? 
                        ncensor.plot.height = 0.25)


cox_final <-  coxph(Surv(OS_MONTHS,OS_status_binary)~cluster, data = cluster_analysis_end)
summary(cox_final)
cox.zph(cox_final) %>%
  plot()

subset_data_final <- subset(cluster_analysis_end, OS_MONTHS <= 60)

# Fit the Cox proportional hazards model using the subsetted data
model_final_2 <- survfit(Surv(OS_MONTHS, OS_status_binary) ~ cluster, data = subset_data_final)

surv_plot_end2<- ggsurvplot(model_final_2, data = subset_data_final, pval = TRUE, pval.method= TRUE,conf.int = TRUE,
                           risk.table = TRUE, risk.table.title = "Risk table",
                           xlab = "Time (in months)", ylab = "Survival Probability",
                           title = "Kaplan-Meier Survival Curve with Events Marked",
                           tables.theme = theme_cleantable(),
                           surv.median.line = "hv",
                           break.time.by = 12,
                           ggtheme = theme_bw(),
                           palette = c(  "#2CA02C", '#FF7F0E'),
                           ncensor.plot = TRUE,      # is there data for this? 
                           ncensor.plot.height = 0.25)

cox_final2 <-  coxph(Surv(OS_MONTHS,OS_status_binary)~cluster, data = subset_data_final)

cox.zph(cox_final2) %>%
  print()


surv_est <- survfit(cox_final)


model_fit_DFS <- survfit(Surv(DFS_MONTHS,DFS_binary) ~ cluster, data = subset_data)

# Create the survival plot for the current cluster
DFS_plot <- ggsurvplot(model_fit_DFS, data = subset_data,
                       pval = TRUE, conf.int = TRUE,risk.table = TRUE, 
                       xlab = "Time (in months)", ylab = "DFS Probability",
                       title = paste("DFS Curve - Disease", i),
                       ggtheme = theme_bw(),
                       palette = cluster_colors_1)

##################################################################
##################################################################
#RF to see what contributes to most important in determining the cluster biologically 

select_end <- cluster_analysis_end[, c(1,2, 3, 24, 26:35)]
tnm_final<- select_end[,c(2,11:14)]

table(select_end$PATH_T_STAGE)
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'TX')]<-NA
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == '[Not Available]')]<-NA
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == '[Discrepancy]')]<-NA

select_end$PATH_T_STAGE<- as.factor(select_end$PATH_T_STAGE)

select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T1')]<-'T1'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T1a')]<-'T1'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T1a1')]<-'T1'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T1b')]<-'T1'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T1b1')]<-'T1'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T1b2')]<-'T1'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T1c')]<-'T1'


select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE== 'T2')]<-'T2'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE== 'T2a1')]<-'T2'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE== 'T2a2')]<-'T2'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T2b')]<-'T2'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE== 'T2c')]<-'T2'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE== 'T2a')]<-'T2'



select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE== 'T3')]<-'T3'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE== 'T3a')]<-'T3'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T3b')]<-'T3'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T3c')]<-'T3'


select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T4')]<-'T4'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T4a')]<-'T4'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'T4b')]<-'T4'
select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE== 'T4d')]<-'T4'

select_end$PATH_T_STAGE[which(select_end$PATH_T_STAGE == 'Tis')]<-5


tnm_final$T0<- select_end$PATH_T_STAGE == "T0"
tnm_final$T1<- select_end$PATH_T_STAGE == "1"
tnm_final$T2<- select_end$PATH_T_STAGE == "2"
tnm_final$T3<- select_end$PATH_T_STAGE == "3"
tnm_final$T4<- select_end$PATH_T_STAGE == "4"
tnm_final$T5<- select_end$PATH_T_STAGE == "5"


table(select_end$PATH_N_STAGE)
select_end$PATH_N_STAGE[which(select_end$PATH_N_STAGE == '[Discrepancy]')]<-NA
select_end$PATH_N_STAGE[which(select_end$PATH_N_STAGE == '[Not Available]')]<-NA

select_end$PATH_N_STAGE<- as.factor(select_end$PATH_N_STAGE)

select_end$PATH_N_STAGE <- ifelse(select_end$PATH_N_STAGE == 'N0', 'N-',
                              ifelse(select_end$PATH_N_STAGE == 'N0 (i-)', 'N-',
                                     ifelse(select_end$PATH_N_STAGE == 'N0 (i+)','N-',
                                            ifelse(select_end$PATH_N_STAGE == 'N0 (mol+)','N-', 
                                                   ifelse(select_end$PATH_N_STAGE == 'N1','N+',
                                                          ifelse(select_end$PATH_N_STAGE == 'N1a','N+',
                                                                 ifelse(select_end$PATH_N_STAGE == 'N1b','N+',
                                                                        ifelse(select_end$PATH_N_STAGE == 'N1c','N+',
                                                                               ifelse(select_end$PATH_N_STAGE == 'N1mi','N+',
                                                                                      ifelse(select_end$PATH_N_STAGE == 'N2', 'N+',
                                                                                             ifelse(select_end$PATH_N_STAGE == 'N2a','N+',
                                                                                                    ifelse(select_end$PATH_N_STAGE == 'N2b','N+',
                                                                                                           ifelse(select_end$PATH_N_STAGE== 'N2c','N+',
                                                                                                                  ifelse(select_end$PATH_N_STAGE == 'N3', 'N+',
                                                                                                                         ifelse(select_end$PATH_N_STAGE == 'N3a','N+',
                                                                                                                                ifelse(select_end$PATH_N_STAGE == 'N3b','N+',
                                                                                                                                       ifelse(select_end$PATH_N_STAGE == 'N3c','N+',
                                                                                                                                              ifelse(select_end$PATH_N_STAGE == 'NX',NA,
                                                                                                                                                     NA))))))))))))))))))
tnm_final$N0<- select_end$PATH_N_STAGE == "1"
tnm_final$N1<- select_end$PATH_N_STAGE == "2"
tnm_final$N2<- select_end$PATH_N_STAGE == "3"
tnm_final$N3<- select_end$PATH_N_STAGE == "4"



table(unique(select_end$PATH_M_STAGE))
select_end$PATH_M_STAGE[which(select_end$PATH_M_STAGE == '[Discrepancy]')]<-NA
select_end$PATH_M_STAGE[which(select_end$PATH_M_STAGE == '[Not Available]')]<-NA
select_end$PATH_M_STAGE[which(select_end$PATH_M_STAGE == '[Not Applicable]')]<-NA

select_end$PATH_M_STAGE<- as.factor(select_end$PATH_M_STAGE)

select_end$PATH_M_STAGE <- ifelse(select_end$PATH_M_STAGE == 'cM0 (i+)','M-',
                              ifelse(select_end$PATH_M_STAGE == 'M0','M-',
                                     ifelse(select_end$PATH_M_STAGE == 'M1','M+',
                                            ifelse(select_end$PATH_M_STAGE == 'M1a','M+',
                                                   ifelse(select_end$PATH_M_STAGE == 'M1b','M+',
                                                          ifelse(select_end$PATH_M_STAGE == 'M1c','M+',
                                                                 ifelse(select_end$PATH_M_STAGE == 'MX',NA,
                                                                        NA)))))))


tnm_final$M0<- select_end$PATH_M_STAGE == "1"
tnm_final$M1<- select_end$PATH_M_STAGE == "2"



table(unique(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE))
select_end$AJCC_PATHOLOGIC_TUMOR_STAGE[which(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == '[Discrepancy]')]<-NA
select_end$AJCC_PATHOLOGIC_TUMOR_STAGE[which(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == '[Not Available]')]<-NA
select_end$AJCC_PATHOLOGIC_TUMOR_STAGE[which(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == '[Not Applicable]')]<-NA
select_end$AJCC_PATHOLOGIC_TUMOR_STAGE[which(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE I')]<- 'Stage I'
select_end$AJCC_PATHOLOGIC_TUMOR_STAGE[which(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE II')]<- 'Stage II'
select_end$AJCC_PATHOLOGIC_TUMOR_STAGE[which(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE III')]<- 'Stage III'
select_end$AJCC_PATHOLOGIC_TUMOR_STAGE[which(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE IV')]<- 'Stage IV'

select_end$AJCC_PATHOLOGIC_TUMOR_STAGE <- ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage I','Stage I',
                                  ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IA','Stage I',
                                         ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IA1','Stage I',
                                                ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IA2','Stage I',
                                                       ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IB','Stage I',
                                                              ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IB1','Stage I',
                                                                     ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IB2','Stage I',
                                                                            ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IC','Stage I',
                                                                                   ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage II','Stage II',
                                                                                          ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIA','Stage II',
                                                                                                 ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIA1','Stage II',
                                                                                                        ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIA2','Stage II',
                                                                                                               ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIB','Stage II',
                                                                                                                      ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIC','Stage II',
                                                                                                                             ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage III','Stage III',
                                                                                                                                    ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIIA','Stage III',
                                                                                                                                           ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIIB','Stage III',
                                                                                                                                                  ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIIC','Stage III',
                                                                                                                                                         ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIIC1','Stage III',
                                                                                                                                                                ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IIIC2','Stage III',
                                                                                                                                                                       ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IV','Stage IV',
                                                                                                                                                                              ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IVA','Stage IV',
                                                                                                                                                                                     ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IVB','Stage IV',
                                                                                                                                                                                            ifelse(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE  == 'Stage IVC','Stage IV',
                                                                                                                                                                                                   NA))))))))))))))))))))))))

tnm_final$Stage1<- select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == "1"
tnm_final$Stage2<- select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == "2"
tnm_final$Stage3<- select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == "3"
tnm_final$Stage4<- select_end$AJCC_PATHOLOGIC_TUMOR_STAGE == "4"


select_end[!(is.na(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE)),]%>%
  group_by(cluster, AJCC_PATHOLOGIC_TUMOR_STAGE) %>%
  summarise(n=n())%>%
  mutate(prop = n/sum(n))%>%
  ggplot(aes(x = AJCC_PATHOLOGIC_TUMOR_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_emt) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of AJCC Stage by Cluster") +
  theme_bw()

1276+1400+1178+619

dim(select_end[!(is.na(select_end$PATH_M_STAGE)),])

1276/4473






select_end$AJCC_PATHOLOGIC_TUMOR_STAGE<- as.factor(select_end$AJCC_PATHOLOGIC_TUMOR_STAGE)

select_end$cluster<- as.factor(select_end$cluster)

TNM_cols<- tnm_final[,c(1,6:21)]

TNM_cols$cluster <- as.factor(TNM_cols$cluster)
TNM_cols$cluster


train_final <- createDataPartition(TNM_cols$cluster, p= 0.7, list = FALSE)
trainset_final <- TNM_cols[train_final,]
validset_final <- TNM_cols[-train_final,]

trainset_final <- na.omit(trainset_final)
validset_final<- na.omit(validset_final)

rf_final <- randomForest(cluster ~ ., data = trainset_final, ntree= 500, importance = TRUE)
importance(rf_final)        
varImpPlot(rf_final, pch=16, col="black")

imp = as.data.frame(importance(rf_final))
imp
imp = cbind(vars=rownames(imp), imp)
imp
imp = imp[order(imp$MeanDecreaseGini),]  
imp$vars = factor(imp$vars, levels=unique(imp$vars))


barplot(imp$MeanDecreaseGini, names.arg=imp$vars,main = 'Importance of the variable using "gini"')



accuracy <- c()
for (i in 1:16) {
  sample_rf <- randomForest(cluster ~ ., data = trainset_final, ntree = 100, mtry = i)
  predValid <- predict(sample_rf,validset_final, type = 'class')
  accuracy[i] <- mean(predValid == validset_final$cluster)
}
accuracy

Dec_tree <- rpart(cluster ~ ., data = trainset_final, method = 'class') #is this the right dependent variable to put through?
rpart.plot(Dec_tree)
summary(Dec_tree)
printcp(Dec_tree)
barplot(summary(Dec_tree)$variable.importance, main="Variable Importance Plot ")
plotcp(Dec_tree)


###############################################################
###############################################################
# going to try the pearson correlation for OS months 



prob_df <- merge(emt_sig, df_test_emt , by = "CFH" )

#cluster = 131, prediciton 132
df_1 <- df_test_rfe[,c(4,131,132)]
df_2 <- emt_sig[,c(1,399)]

prob_df<- merge(df_1, df_2, by= 'PLAU')

df_3 <- cluster_analysis_end[,c(2,28,30:35)]

prob_bio <- merge(prob_df,df_3, by = 'Sample')


prob_bio$prediction_prob <- as.numeric(prob_bio$prediction_prob)
prob_bio$OS_MONTHS <- as.numeric(prob_bio$OS_MONTHS)
complete_cases <- complete.cases(prob_bio$prediction_prob, prob_bio$OS_MONTHS)
correlation <- cor(prob_bio$prediction_prob[complete_cases], prob_bio$OS_MONTHS[complete_cases])

plot(prob_bio$prediction_prob, prob_bio$OS_MONTHS,
     pch = 16,
     xlab = "Prediction Probability",
     ylab = "OS in Months")

# Add a text box with the correlation coefficient
text(x = max(prob_bio$prediction_prob, na.rm = TRUE), y = max(prob_bio$OS_bi, na.rm = TRUE),
     paste("Correlation:", round(correlation, 2)), adj = c(1, 0))

# Add a trend line
abline(lm(prob_bio$OS_status_binary ~ prob_bio$prediction_prob), col = "blue")

correlation_above <- cor(above_threshold$prediction_prob, above_threshold$OS_MONTHS, use = "complete.obs")
correlation_below <- cor(below_threshold$prediction_prob, below_threshold$OS_MONTHS, use = "complete.obs")

# Print the correlation coefficients
cat("Correlation coefficient for data above threshold: ", correlation_above, "\n")
cat("Correlation coefficient for data below threshold: ", correlation_below, "\n")


spearman_corr <- cor.test(prob_bio$OS_status_binary, prob_bio$prediction_prob, method = "spearman")

# Print the correlation coefficient and p-value
cat("Spearman's rank correlation coefficient: ", spearman_corr$estimate, "\n")
cat("p-value: ", spearman_corr$p.value, "\n")



###filter for possible biomarker?

biomarker<- emt_sig[,c(1134, 29)]

biomarker %>%
  filter( ACP3 >=12)%>%
  count(Disease)%>%
  print()

# Create a vector to store colors based on the prediction_prob values
point_colors <- ifelse(prob_bio$prediction_prob > 0.5, "green", "orange")

# Plot the data with colored points
plot(prob_bio$prediction_prob, prob_bio$OS_MONTHS,
     pch = 16,
     col = point_colors,  # Use the point_colors vector for coloring
     xlab = "Prediction Probability",
     ylab = "OS in Months")

# Calculate the correlation coefficient
correlation <- cor(prob_bio$prediction_prob, prob_bio$OS_MONTHS, use = "complete.obs")

# Add a text box with the correlation coefficient
text(x = max(prob_bio$prediction_prob, na.rm = TRUE), y = max(prob_bio$OS_MONTHS, na.rm = TRUE),
     paste("Correlation:", round(correlation, 2)), adj = c(1, 0))
