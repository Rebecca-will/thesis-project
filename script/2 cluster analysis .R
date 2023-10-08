attach(TNM_2)
library(dplyr)

#total path + 

dim(Clara_T3)
colnames(Clara_T3)
dim(stacked_data)
colnames(big_df)
dim(big_df)


big_df <- merge(Clara_T3, stacked_data, by ='PATIENT_ID')
dim(big_df)
colnames(big_df)

missing_1 <- setdiff(Total_path$PATIENT_ID, big_df$PATIENT_ID)
print(missing_1)


big_df <- merge(big_df, Total_path, by = 'PATIENT_ID')

colnames(stacked_data)

df <- merge(stacked_data,Total_path, by = 'PATIENT_ID')




# dim(df)
# length(cutclcolhc)
# duplicated_rows<- duplicated(df)


# print(df[duplicated_rows,])
# 
# colnames(select_comb)
# 
# select_comb <- subset(combined_data, select = c(1,24,25))
# 
# df_final <- merge(df,select_comb, by = "PATIENT_ID")

table(unique(df_final$disease))


table(big_df$OS_STATUS)


big_df$OS_STATUS[which(big_df$OS_STATUS =='1:DECEASED')] <- 'DECEASED'
big_df$OS_STATUS[which(big_df$OS_STATUS =='0:LIVING')] <- 'LIVING'
big_df$OS_STATUS[which(big_df$OS_STATUS =='[Not Available]')] <- NA

big_df$OS_status_binary <- ifelse(big_df$OS_STATUS == "DECEASED", 1,
                                         ifelse(big_df$OS_STATUS == "LIVING", 0,
                                                NA))

big_df$cluster <- ifelse(big_df$cluster == '#1F77B4',1,
                                ifelse(big_df$cluster == '#2CA02C',2,
                                       3))


big_df$OS_MONTHS<- as.numeric(big_df$OS_MONTHS)
colnames(big_df)
table(df$disease)


cluster_2_analysis <- big_df[!(big_df$disease == 'lgg'),]
cluster_2_analysis <- cluster_2_analysis[!(cluster_2_analysis$disease == 'gbm'),]
cluster_2_analysis <- cluster_2_analysis[!(cluster_2_analysis$disease == 'skcm'),]
cluster_2_analysis <- cluster_2_analysis[!(cluster_2_analysis$disease == 'tgct'),]



table(cluster_2_analysis$disease)




#cluster_duo2 <- TNM_DF_final[!(TNM_DF_final$cluster=='#1F77B4'),]#check colors as data is subsetted
#dim(cluster_duo2)

write.csv(cluster_2_analysis, file = 'No_mes.csv', row.names = FALSE)
write.csv(big_df, file = 'All_cols.csv', row.names = FALSE)

#########
#load from here the df

cluster_2_analysis<- No_mes

table(cluster_2_analysis$PATH_T_STAGE)

cluster_2_analysis$PATH_T_STAGE[which(cluster_2_analysis$PATH_T_STAGE == '[Not Available]')]<-NA
cluster_2_analysis$PATH_T_STAGE[which(cluster_2_analysis$PATH_T_STAGE == '[Discrepancy]')]<-NA
cluster_2_analysis$PATH_T_STAGE[which(cluster_2_analysis$PATH_T_STAGE == 'TX')]<-NA



table(cluster_2_analysis$PATH_M_STAGE)

cluster_2_analysis$PATH_M_STAGE[which(cluster_2_analysis$PATH_M_STAGE == '[Not Available]')]<-NA
cluster_2_analysis$PATH_M_STAGE[which(cluster_2_analysis$PATH_M_STAGE == '[Discrepancy]')]<-NA
cluster_2_analysis$PATH_M_STAGE[which(cluster_2_analysis$PATH_M_STAGE == 'MX')]<-NA

table(cluster_2_analysis$PATH_N_STAGE)

cluster_2_analysis$PATH_N_STAGE[which(cluster_2_analysis$PATH_N_STAGE == '[Discrepancy]')]<-NA
cluster_2_analysis$PATH_N_STAGE[which(cluster_2_analysis$PATH_N_STAGE == '[Not Available]')]<-NA
cluster_2_analysis$PATH_N_STAGE[which(cluster_2_analysis$PATH_N_STAGE == 'NX')]<-NA

table(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE)

cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE[which(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE == '[Discrepancy]')] <-NA
cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE[which(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE == '[Not Available]')] <-NA
cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE[which(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE == 'Stage X')] <-NA
cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE[which(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE I')] <- 'Stage I'
cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE[which(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE II')] <- 'Stage II'
cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE[which(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE III')] <- 'Stage III'
cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE[which(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE == 'STAGE IV')] <- 'Stage IV'

########################################################################################################################################################################
########################################################################################################################################################################

#survival analysis without mes.

unique(cluster_2_analysis$cluster)

TNM_DF_final<- na.omit(TNM_DF_final)



dim(cluster_2_analysis)#have removed the mes diseases so don't cry 


sum(table(cluster_2_analysis$cluster))
cluster_2_analysis$cluster <- ifelse(cluster_2_analysis$cluster == "#1F77B4", '1-gbm/lgg',
                        ifelse(cluster_2_analysis$cluster == "#2CA02C", '2-classical',
                               ifelse(cluster_2_analysis$cluster == '#FF7F0E','3-no emt',
                                      NA)))

mes_left<- cluster_2_analysis[(cluster_2_analysis$cluster == '1-gbm/lgg'),]
as.vector(unique(mes_left$disease))
colnames(mes_left)

disease_proprtions <- prop.table(table(mes_left$Disease))*100




cluster_2_analysis$cluster <- as.factor(cluster_2_analysis$cluster)


cluster_2_analysis$OS_MONTHS<- as.numeric(cluster_2_analysis$OS_MONTHS)



duo_colours <- c('#2CA02C','#FF7F0E')



model_surv_duo<-survfit(Surv(OS_MONTHS,OS_status_binary)~cluster, data = cluster_2_analysis)

colours_2 = unique(cluster_2_analysis[order(cluster_2_analysis$cluster),]$cluster_colour)

surv_plot_duo<- ggsurvplot(model_surv_duo, data = cluster_2_analysis, pval = TRUE, pval.method= TRUE,conf.int = TRUE,
                           risk.table = TRUE, risk.table.title = "Risk table",
                           xlab = "Time (in months)", ylab = "Survival Probability",
                           title = "Kaplan-Meier Survival Curve with Events Marked",
                           tables.theme = theme_cleantable(),
                           surv.median.line = "hv",
                           break.time.by = 12,
                           ggtheme = theme_bw(),
                           palette = colours_2,
                           ncensor.plot = TRUE,       
                           ncensor.plot.height = 0.25)

model_surv_duo_cox <- coxph(Surv(OS_MONTHS,OS_status_binary)~cluster, data = cluster_2_analysis)

summary(model_surv_duo_cox)

cox.zph(model_surv_duo_cox)%>%
  plot()

dim(cluster_2_analysis)
dim(emt_no_emt)





#is the survival difference between cluster 2 and 3 significant i.e. emt and no emt 

emt_no_emt <- cluster_2_analysis[!(cluster_2_analysis$cluster=='1-gbm/lgg'),]
unique(emt_no_emt$cluster)

emt_no_emt$cluster <- droplevels(emt_no_emt$cluster)
emt_no_emt$cluster<- as.factor(emt_no_emt$cluster)

model_surv_emt2<-survfit(Surv(OS_MONTHS,OS_status_binary)~cluster, data = emt_no_emt)

colours_emt = unique(emt_no_emt[order(emt_no_emt$cluster),]$cluster_colour)

surv_plot_emt<- ggsurvplot(model_surv_emt2, data = emt_no_emt, pval = TRUE, pval.method= TRUE,conf.int = TRUE,
                           risk.table = TRUE, risk.table.title = "Risk table",
                           xlab = "Time (in months)", ylab = "Survival Probability",
                           title = "Kaplan-Meier Survival Curve with Events Marked",
                           tables.theme = theme_cleantable(),
                           surv.median.line = "hv",
                           break.time.by = 12,
                           ggtheme = theme_bw(),
                           palette = colours_emt,
                           ncensor.plot = TRUE,       
                           ncensor.plot.height = 0.25)

model_surv_emt<- coxph(Surv(OS_MONTHS,OS_status_binary)~cluster + Disease, data = emt_no_emt)


cox_fit <- survfit(model_surv_emt)
plot(cox_fit)

survivalAnalysis::write_survival(model_surv_emt)

library(survminer)
forest_plot(model_surv_emt)

ggforest(model_surv_emt,
         main = "Forest Plot of Cox Model",
         xlab = "Hazard Ratio (HR)",
         xlim = c(0.1, 5),   # Optional: Set the range for the x-axis
         col = "royalblue",  # Optional: Set the color of the HR boxes
         cex = 0.8)  



coxphsummary<- summary(model_surv_emt)

a<- cox_as_data_frame(coxphsummary,
                  unmangle_dict = NULL,
                  factor_id_sep = ":",
                  sort_by = NULL)
forest_plot(a)

result_cox <- cox_multivariate(
  formula = Surv(OS_MONTHS, OS_status_binary) ~ cluster + Disease,
  data = emt_no_emt
)

result_cox <- SurvivalAnalysisResult(model_cox)

b<- analyze_multivariate(
  emt_no_emt,
  OS_status_binary,
  c(cluster + Disease),
  strata = NULL,
  covariate_name_dict = NULL,
  covariate_label_dict = NULL,
  reference_level_dict = NULL,
  sort_frame_by = vars(HR)
)

library(magrittr)
library(dplyr)

# Perform the multivariate survival analysis
result_2<-emt_no_emt%>%
  analyse_multivariate(vars(OS_MONTHS, OS_status_binary),
                       vars(cluster))


mean_OS_MONTHS_paad <- mean(emt_no_emt$OS_MONTHS[emt_no_emt$Disease == "ov"], na.rm = TRUE)
mean_OS_MONTHS_prad <- mean(emt_no_emt$OS_MONTHS[emt_no_emt$Disease == "prad"], na.rm = TRUE)

forest_plot(result_2)
forest_plot(result,
            endpoint_labeller = c(time="OS"),
            labels_displayed = c("endpoint", "factor", "n"),
            ggtheme = ggplot2::theme_bw(base_size = 10),
            relative_widths = c(1, 1.5, 1),
            HR_x_breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2))


summary(model_surv_emt)
cox.zph(model_surv_emt)%>%
  plot()


dev.off()
table(emt_no_emt$DFS_STATUS)
#reasign the different variables to one 
emt_no_emt$DFS_STATUS[which(emt_no_emt$DFS_STATUS =='1:Recurred/Progressed')] <- 'Recurred/Progressed'
emt_no_emt$DFS_STATUS[which(emt_no_emt$DFS_STATUS =='0:DiseaseFree')] <- 'DiseaseFree'
emt_no_emt$DFS_STATUS[which(emt_no_emt$DFS_STATUS =='[Not Available]')] <- NA

table(emt_no_emt$DFS_STATUS)


emt_no_emt$DFS_binary <- ifelse(emt_no_emt$DFS_STATUS == "Recurred/Progressed", 1,
                                   ifelse(emt_no_emt$DFS_STATUS == "DiseaseFree", 0,
                                          NA))
emt_no_emt$DFS_MONTHS<- as.numeric(emt_no_emt$DFS_MONTHS)

model_fit_DFS_noemt <- survfit(Surv(DFS_MONTHS,DFS_binary) ~ cluster, data = emt_no_emt)

# Create the survival plot for the current cluster
DFS_plot <- ggsurvplot(model_fit_DFS_noemt, data = emt_no_emt,
                       pval = TRUE, conf.int = TRUE,risk.table = TRUE, 
                       xlab = "Time (in months)", ylab = "DFS Probability",
                       title = paste("DFS Curve - Disease", i),
                       ggtheme = theme_bw(),
                       palette = colours_emt)


print(plot_list_emt)


plot_list_duo <- list()


# Iterate through each disease cluster
for (i in unique(cluster_2_analysis$Disease)) {
  # Subset the data for the current cluster
  subset_data <- subset(cluster_2_analysis, Disease == i)
  
  colours_duo = unique(cluster_2_analysis[order(cluster_2_analysis$cluster),]$cluster_colour)
  
  # Fit the survival model for the current cluster
  model_fit_surv <- survfit(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = subset_data)
  
  # Create the survival plot for the current cluster
  surv_plot <- ggsurvplot(model_fit_surv, data = subset_data,
                          pval = TRUE, conf.int = TRUE,risk.table = TRUE,
                          xlab = "Time (in months)", ylab = "Survival Probability",
                          title = paste("Survival Curve - Disease", i),
                          ggtheme = theme_bw(),
                          palette = colours_duo)
  
  #print which disease we are in 
  print(i)
  
  
  # Add the survival plot to the list
  plot_list_duo[[i]] <- surv_plot
  
  
}
########################################################################################################################
#come back to later
list_cox <- list()

plot(model_fit_cox)

# Iterate through each disease cluster
for (i in unique(cluster_2_analysis$Disease)) {
  # Subset the data for the current cluster
  subset_data <- subset(cluster_2_analysis, Disease == i)
  
  
  # Fit the survival model for the current cluster
  model_fit_cox <-coxph(Surv(OS_MONTHS,OS_status_binary)~cluster, data = emt_no_emt)
  
  # Create the survival plot for the current cluster
  surv_plot <- ggsurvplot(model_fit_cox, data = subset_data,
                          pval = TRUE, conf.int = TRUE, risk.table = TRUE,
                          hazard = TRUE, # Set hazard argument to TRUE
                          xlab = "Time (in months)", ylab = "Hazard Ratio",
                          title = paste("Hazard Plot - Disease", i),
                          ggtheme = theme_bw(),
                          palette = colours_duo)
  
  #print which disease we are in 
  print(i)
  
  
  # Add the survival plot to the list
  list_cox[[i]] <- surv_plot
  
  
}

plot_list_emt <- list()

# Iterate through each disease cluster
for (i in unique(emt_no_emt$Disease)) {
  # Subset the data for the current cluster
  subset_data <- subset(emt_no_emt, Disease == i)
  
  colours_emt = unique(emt_no_emt[order(emt_no_emt$cluster),]$cluster_colour)
  
  # Fit the survival model for the current cluster
  model_fit_surv <- survfit(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = subset_data)
  
  # Create the survival plot for the current cluster
  surv_plot <- ggsurvplot(model_fit_surv, data = subset_data,
                          pval = TRUE, conf.int = TRUE,risk.table = TRUE,
                          xlab = "Time (in months)", ylab = "Survival Probability",
                          title = paste("Survival Curve - Disease", i),
                          ggtheme = theme_bw(),
                          palette = colours_emt)
  
  #print which disease we are in 
  print(i)
  
  
  # Add the survival plot to the list
  plot_list_emt[[i]] <- surv_plot
  
  
}



print(plot_list_emt)
########################################################################################################################
########################################################################################################################
#individual diseases

sub_acc_duo <- subset(cluster_2_analysis, disease == 'acc')



acc_model_duo <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_acc_duo)
acc_model_duo 
summary(acc_model_duo)  

cox.zph(acc_model_duo)%>%
  plot()



sub_kich_duo <- subset(cluster_2_analysis, disease == 'kich')

kich_model_duo <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_kich_duo)
kich_model_duo 
summary(kich_model_duo)  

cox.zph(kich_model_duo)



sub_kirc_duo <- subset(cluster_2_analysis, disease == 'kirc')

kirc_model_duo <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_kirc_duo)
kirc_model_duo 
summary(kirc_model_duo)  

cox.zph(kirc_model_duo)


sub_kirp_duo <- subset(cluster_2_analysis, disease == 'kirp')

kirp_model_duo <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_kirp_duo)
kirp_model_duo 
summary(kirp_model_duo)  




sub_paad_duo <- subset(cluster_2_analysis, disease == 'paad')

paad_model_duo <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_paad_duo)
paad_model_duo 
summary(paad_model_duo)  



sub_esca_duo <- subset(cluster_2_analysis, disease == 'esca')

esca_model_duo <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_esca_duo)
paad_model_duo 
summary(esca_model_duo)  

cox.zph(esca_model_duo)




sub_ov_duo <- subset(cluster_2_analysis, disease == 'ov')

ov_model_duo <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_ov_duo)
summary(ov_model_duo)  

cox.zph(ov_model_duo)


########################################################################################################################
########################################################################################################################
#individual in EMT vs no EMT


sub_lihc_emt <- subset(emt_no_emt, disease == 'lihc')

lihc_model_emt <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_lihc_emt)
summary(lihc_model_emt)  


sub_paad_emt <- subset(emt_no_emt, disease == 'paad')

paad_model_emt <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_paad_emt)
summary(paad_model_emt)  


sub_kirp_emt <- subset(emt_no_emt, disease == 'kirp')

kirp_model_emt <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_kirp_emt)
summary(kirp_model_emt)  

cox.zph(kirp_model_emt)

########################################################################################################################
########################################################################################################################


list(unique(cluster_2_analysis$cluster))

cluster_2_analysis<- na.omit(cluster_2_analysis)

#write.csv(Cluster_duo, file = '2 Clusters', row.names = FALSE)


cluster_duo2 %>%
  count(AJCC_PATHOLOGIC_TUMOR_STAGE, cluster) %>%
  mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = AJCC_PATHOLOGIC_TUMOR_STAGE, y = n, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#2CA02C", '#FF7F0E')) +
  labs(x = "Stage", y = "Frequency", fill = "Cluster") +
  ggtitle("Frequency of Stage by Cluster") +
  theme_bw()

colnames(Cluster_duo)
cluster_2_analysis [!(is.na(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE)),]%>%
  filter(cluster %in% c("1-gbm/lgg","2-classical", "3-no emt")) %>%
  group_by(cluster, AJCC_PATHOLOGIC_TUMOR_STAGE) %>%
  summarise(n = n()) %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = AJCC_PATHOLOGIC_TUMOR_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_duo) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of AJCC Stage by Cluster") +
  theme_bw()


not_app <- cluster_2_analysis[(cluster_2_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE == '[Not Applicable]'),]
testing <- cluster_2_analysis[(cluster_2_analysis$disease == 'prad'),]


table(not_app$disease)


cluster_2_analysis[!(is.na(cluster_2_analysis$PATH_T_STAGE)),] %>%
  filter(cluster %in% c("1-gbm/lgg","2-classical", "3-no emt")) %>%
  group_by(cluster, PATH_T_STAGE) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = PATH_T_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_duo) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH T Stage by Cluster") +
  theme_bw()

cluster_2_analysis[!(is.na(cluster_2_analysis$PATH_N_STAGE)),] %>%
  filter(cluster %in% c("1-gbm/lgg","2-classical", "3-no emt")) %>%
  group_by(PATH_N_STAGE, cluster) %>%
  summarise(n = n(),) %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = PATH_N_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_duo) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH N Stage by Cluster") +
  theme_bw()

cluster_2_analysis[!(is.na(cluster_2_analysis$PATH_M_STAGE)),] %>%
  filter(cluster %in% c("1-gbm/lgg","2-classical", "3-no emt")) %>%
  group_by(PATH_M_STAGE, cluster) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = PATH_M_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_duo) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH M Stage by Cluster") +
  theme_bw()

not_app_M <- cluster_2_analysis[(cluster_2_analysis$PATH_M_STAGE == '[Not Applicable]'),]

table(not_app_M$disease)

########################################################################################################################
########################################################################################################################




emt_no_emt[!(is.na(emt_no_emt$AJCC_PATHOLOGIC_TUMOR_STAGE)),] %>%
  filter(cluster %in% c("2-classical", "3-no emt")) %>%
  group_by(cluster, AJCC_PATHOLOGIC_TUMOR_STAGE) %>%
  summarise(n = n()) %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = AJCC_PATHOLOGIC_TUMOR_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_emt) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of AJCC Stage by Cluster") +
  theme_bw()


emt_no_emt[!(is.na(emt_no_emt$PATH_N_STAGE)),] %>%
  filter(cluster %in% c("2-classical", "3-no emt")) %>%
  na.omit()%>%
  group_by(cluster, PATH_N_STAGE) %>%
  summarise(n = n())  %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x = PATH_N_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_emt) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH N Stage by Cluster") +
  theme_bw()




emt_no_emt[!(is.na(emt_no_emt$PATH_M_STAGE)),]%>%
  group_by(cluster, PATH_M_STAGE) %>%
  summarise(n=n())%>%
  mutate(prop = n/sum(n))%>%
  ggplot(aes(x = PATH_M_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_emt) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH M Stage by Cluster") +
  theme_bw()
  

emt_no_emt[!(is.na(emt_no_emt$PATH_N_STAGE)),]%>%
  group_by(cluster, PATH_N_STAGE) %>%
  summarise(n=n())%>%
  mutate(prop = n/sum(n))%>%
  ggplot(aes(x = PATH_N_STAGE, y = prop, fill = cluster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = colours_emt) +
  labs(x = "Stage", y = "Proportion", fill = "Cluster") +
  ggtitle("Proportion of PATH N Stage by Cluster") +
  theme_bw()




table(emt_no_emt[!(is.na(emt_no_emt$PATH_T_STAGE)),]$cluster)

275/4339

tab_emt <- emt_no_emt %>%
  filter(cluster %in% c("2-classical", "3-no emt")) %>%
  na.omit()%>%
  group_by(PATH_T_STAGE, cluster) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(cluster) %>%
  mutate(prop = n/sum(n))

sum(tab_emt$n)

emt_no_emt$PATH_N_STAGE <- ifelse(emt_no_emt$PATH_N_STAGE == 'N0', 'N-',
                                  ifelse(emt_no_emt$PATH_N_STAGE == 'N0 (i-)', 'N-',
                                         ifelse(emt_no_emt$PATH_N_STAGE == 'N0 (i+)','N-',
                                                ifelse(emt_no_emt$PATH_N_STAGE == 'N0 (mol+)','N-', 
                                                       ifelse(emt_no_emt$PATH_N_STAGE == 'N1','N+',
                                                              ifelse(emt_no_emt$PATH_N_STAGE == 'N1a','N+',
                                                                     ifelse(emt_no_emt$PATH_N_STAGE == 'N1b','N+',
                                                                            ifelse(emt_no_emt$PATH_N_STAGE == 'N1c','N+',
                                                                                   ifelse(emt_no_emt$PATH_N_STAGE == 'N1mi','N+',
                                                                                          ifelse(emt_no_emt$PATH_N_STAGE == 'N2', 'N+',
                                                                                                 ifelse(emt_no_emt$PATH_N_STAGE == 'N2a','N+',
                                                                                                        ifelse(emt_no_emt$PATH_N_STAGE == 'N2b','N+',
                                                                                                               ifelse(emt_no_emt$PATH_N_STAGE== 'N2c','N+',
                                                                                                                      ifelse(emt_no_emt$PATH_N_STAGE == 'N3', 'N+',
                                                                                                                             ifelse(emt_no_emt$PATH_N_STAGE == 'N3a','N+',
                                                                                                                                    ifelse(emt_no_emt$PATH_N_STAGE == 'N3b','N+',
                                                                                                                                           ifelse(emt_no_emt$PATH_N_STAGE == 'N3c','N+',
                                                                                                                                                  ifelse(emt_no_emt$PATH_N_STAGE == 'NX',NA,
                                                                                                                                                         NA))))))))))))))))))



emt_no_emt$PATH_M_STAGE <- ifelse(emt_no_emt$PATH_M_STAGE == 'cM0 (i+)','M-',
                                  ifelse(emt_no_emt$PATH_M_STAGE == 'M0','M-',
                                         ifelse(emt_no_emt$PATH_M_STAGE == 'M1','M+',
                                                ifelse(emt_no_emt$PATH_M_STAGE == 'M1a','M+',
                                                       ifelse(emt_no_emt$PATH_M_STAGE == 'M1b','M+',
                                                              ifelse(emt_no_emt$PATH_M_STAGE == 'M1c','M+',
                                                                     ifelse(emt_no_emt$PATH_M_STAGE == 'MX',NA,
                                                                            NA)))))))



275/4646
188/4646
4646+88

188/4734

159/2542

188/x = 0.064


188/0.064

