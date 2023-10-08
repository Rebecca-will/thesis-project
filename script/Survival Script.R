#survival analysis 
library(survival)
library(survminer)
library(ggplot2)


combined_data <- merge(Clara_T3, stacked_data[2:7], by = 'PATIENT_ID')
attach(combined_data)
#im going to have to add a column with the associated clusters. Recode as levels. need to wait until I get the data from emma 
#combined_data$cluster <- cutclcolhc


# Create a new column for cluster colors in combined_data based on the cluster values
#combined_data$cluster_color <- cluster_colors[combined_data$cluster]

#1= deceased, 0= living 

table(combined_data$OS_STATUS)

combined_data$OS_STATUS[which(combined_data$OS_STATUS =='1:DECEASED')] <- 'DECEASED'
combined_data$OS_STATUS[which(combined_data$OS_STATUS =='0:LIVING')] <- 'LIVING'
combined_data$OS_STATUS[which(combined_data$OS_STATUS =='[Not Available]')] <- NA

combined_data$OS_status_binary <- ifelse(combined_data$OS_STATUS == "DECEASED", 1,
                    ifelse(combined_data$OS_STATUS == "LIVING", 0,
                           NA))

combined_data$OS_STATUS[which(combined_data$OS_STATUS =='[Not Available]')] <- NA
combined_data$OS_MONTHS[which(combined_data$OS_MONTHS =='[Not Available]')] <- NA

table(is.na(combined_data$OS_STATUS))

table(combined_data$OS_status_binary)


combined_data$OS_MONTHS <- as.numeric(combined_data$OS_MONTHS)

combined_data$cluster <- ifelse(combined_data$cluster == '#1F77B4',1,
                                ifelse(combined_data$cluster == '#2CA02C',2,
                                       3))

table(combined_data$cluster)


class(combined_data$cluster)
class(combined_data$OS_MONTHS)
class(combined_data$OS_status_binary)

combined_data$cluster <- as.factor(combined_data$cluster)



combined_data$OS_MONTHS<- as.numeric(combined_data$OS_MONTHS)

combined_data$cluster<- as.factor(combined_data$cluster)

cluster_colors_1 <- c('#1F77B4', '#2CA02C', '#FF7F0E')

model_fit_surv<-survfit(Surv(OS_MONTHS,OS_status_binary)~cluster, data = cluster_2_analysis)
model_fit_surv


cluster_2_analysis$cluster<- cluster_2_analysis[!(cluster_2_analysis$cluster == '1'),]

cluster_2_analysis$cluster <- ifelse(cluster_2_analysis$cluster == '#1F77B4',1,
                                ifelse(cluster_2_analysis$cluster == '#2CA02C',2,
                                       3))

Surv_disease <- coxph(Surv(OS_MONTHS, OS_status_binary)~ cluster+disease, data = cluster_2_analysis)

ggsurvplot(Surv_disease, data = combined_data, risk.table = TRUE, pval = TRUE, conf.int = TRUE, ggtheme = theme_minimal())



?survfit.formula
cox_surv <- coxph(Surv(OS_MONTHS,OS_status_binary)~cluster, data = combined_data)
summary(cox_surv)
cox_surv

# combined_data<-combined_data$cluster_color[which(cluster=='blue')]<-'#1F77B4'
# combined_data$cluster_color[which(cluster=='green')]<-'#2CA02C'
# combined_data$cluster_color[which(cluster=='orange')]<-'#FF7F0E'

survdiff(Surv(OS_MONTHS,OS_status_binary)~cluster, data = combined_data)


# Plot the Kaplan-Meier survival curve with events marked, the main plot 

surv_plot2<- ggsurvplot(Surv_disease, data = combined_data, pval = TRUE, pval.method= TRUE,conf.int = TRUE,
                       risk.table = TRUE, risk.table.title = "Risk table",
                       xlab = "Time (in months)", ylab = "Survival Probability",
                       title = "Kaplan-Meier Survival Curve with Events Marked",
                       tables.theme = theme_cleantable(),
                       surv.median.line = "hv",
                       break.time.by = 12,
                       ggtheme = theme_bw(),
                       palette = rainbow(n=60),
                       ncensor.plot = TRUE,      # is there data for this? 
                       ncensor.plot.height = 0.25)

?ggsurvplot


surv_plot2


unique(combined_data$cluster)
class(combined_data$cluster)
levels(combined_data$cluster)



main_tph <- cox.zph(cox_surv) %>%
  print()

plot(cox.zph(cox_surv))

dev.off()

ggforest(cox_surv, data = combined_data)

?surv_fit()
#################################################################################################################################
#################################################################################################################################

#plot for prgoression free survival/relapse.
attach(combined_data)

table(combined_data$DFS_STATUS)
#reasign the different variables to one 
combined_data$DFS_STATUS[which(combined_data$DFS_STATUS =='1:Recurred/Progressed')] <- 'Recurred/Progressed'
combined_data$DFS_STATUS[which(combined_data$DFS_STATUS =='0:DiseaseFree')] <- 'DiseaseFree'

table(combined_data$DFS_STATUS)


combined_data$DFS_binary <- ifelse(combined_data$DFS_STATUS == "Recurred/Progressed", 1,
                                         ifelse(combined_data$DFS_STATUS == "DiseaseFree", 0,
                                                NA))


table(combined_data$DFS_binary)
class(combined_data$DFS_MONTHS)
class(combined_data$DFS_binary)

combined_data$DFS_MONTHS <- as.numeric(combined_data$DFS_MONTHS)


model_fit_DFS <- survfit(Surv(DFS_MONTHS,DFS_binary)~cluster, data = combined_data)
model_fit_DFS

DFS_plot<- ggsurvplot(model_fit_DFS, data = combined_data, pval = TRUE, conf.int = TRUE,
                        risk.table = TRUE, risk.table.title = "Risk table",
                        xlab = "Time (in months)", ylab = "DFS Probability",
                        title = "Kaplan-Meier DFS Curve with Events Marked",
                        tables.theme = theme_cleantable(),
                        surv.median.line = "hv",
                        break.time.by = 12,
                        ggtheme = theme_bw(),
                        palette = cluster_colors_1,
                        ncensor.plot = TRUE,      # is there data for this? 
                        ncensor.plot.height = 0.25)

cox_DFS <- coxph(Surv(DFS_MONTHS,DFS_binary)~cluster, data = combined_data)
summary(cox_DFS)


HRDFS_tph <- cox.zph(cox_DFS) %>%
  print()

ggforest(cox_DFS, data = combined_data)

survdiff(Surv(DFS_MONTHS,DFS_binary)~cluster, data = combined_data)
DFS_plot

plot(cox.zph(cox_DFS))






#################################################################################################################################
#################################################################################################################################

#this is OS for the entire data set. 
model_fit_surv_disease<-survfit(Surv(OS_MONTHS,OS_status_binary)~Disease, data = combined_data)
model_fit_surv_disease

surv_plot_disease<- ggsurvplot(model_fit_surv_disease, data = combined_data, pval = TRUE, conf.int = TRUE,
                       risk.table = TRUE, risk.table.title = "Risk table",
                       xlab = "Time (in months)", ylab = "Survival Probability",
                       title = "Kaplan-Meier Survival Curve with Events Marked",
                       tables.theme = theme_cleantable(),
                       surv.median.line = "hv",
                       break.time.by = 12,
                       ggtheme = theme_bw(),
                       palette = my_colors)

surv_plot_disease # will we just choose the ones of interest? 

length(unique(combined_data$Disease))
length(diseasecolors_set)


#################################################################################################################################
#################################################################################################################################
#cluster 1 breakdown of all diseases

combined_data_clus_1 <- subset(combined_data, combined_data$cluster == '1')

model_fit_surv_disease_clus1<-survfit(Surv(OS_MONTHS,OS_status_binary)~Disease, data = combined_data)
#this is all the diseases, same as above
surv_plot_disease_clus1<- ggsurvplot(model_fit_surv_disease_clus1, data = combined_data, pval = TRUE, conf.int = TRUE,
                               risk.table = TRUE, risk.table.title = "Risk table",
                               xlab = "Time (in months)", ylab = "Survival Probability",
                               title = "Kaplan-Meier Survival Curve with Events Marked",
                               tables.theme = theme_cleantable(),
                               surv.median.line = "hv",
                               break.time.by = 12,
                               ggtheme = theme_bw(),
                               palette = my_colors)
surv_plot_disease_clus1

#this was to compare GBM and LGG
brain_comb_1 <- subset(combined_data, combined_data$disease =='gbm')
brain_comb_2 <- subset(combined_data, combined_data$disease =='lgg')
brain_comb <- rbind.data.frame(brain_comb_1,brain_comb_2)

brain_suv <- survfit(Surv(OS_MONTHS,OS_status_binary) ~Disease, data = brain_comb)
brain_suv

surv_plot_disease_brain<- ggsurvplot(brain_suv, data = brain_comb, pval = TRUE, conf.int = TRUE,
                                     risk.table = TRUE, risk.table.title = "Risk table",
                                     xlab = "Time (in months)", ylab = "Survival Probability",
                                     title = "Kaplan-Meier Survival Curve with Events Marked",
                                     tables.theme = theme_cleantable(),
                                     surv.median.line = "hv",
                                     break.time.by = 12,
                                     ggtheme = theme_bw(),
                                     palette = my_colors)
surv_plot_disease_brain

204/12

#################################################################################################################################
#################################################################################################################################

combined_data<- na.omit(combined_data)

subset_data <- subset(combined_data, Disease == 'paad')

plot_list <- list()

# Iterate through each disease cluster
for (i in unique(combined_data$Disease)) {
  # Subset the data for the current cluster
  
  subset_data <- subset(combined_data, Disease == i)
  
  # Fit the survival model for the current cluster
  model_fit_surv <- survfit(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = subset_data)
  
  colours= unique(subset_data[order(subset_data$cluster),]$cluster_colour)#reorder the df according to the cluster number
  
  # Create the survival plot for the current cluster
  surv_plot <- ggsurvplot(model_fit_surv, data = subset_data, conf.int = TRUE, risk.table = TRUE,
                          xlab = "Time (in months)", ylab = "Survival Probability",
                          title = paste("Survival Curve - Disease", i),
                          ggtheme = theme_bw(),
                          palette = colours)
  #print which disease we are in 
  print(i)
  
  # Add the survival plot to the list
  plot_list[[i]] <- surv_plot
  
  
}

str(subset_data$cluster)


plot_list <- list()

print(plot_list)

# Iterate through each disease cluster
for (i in unique(combined_data$disease)) {
  # Subset the data for the current cluster
  i= 'acc'
  subset_data <- subset(combined_data, disease == i)
  
  # Fit the survival model for the current cluster
  model_fit_surv <- survfit(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = subset_data)
  
  # Create the survival plot for the current cluster
  surv_plot <- ggsurvplot(model_fit_surv, data = subset_data,
                          pval = TRUE, conf.int = TRUE,
                          xlab = "Time (in months)", ylab = "Survival Probability",
                          title = paste("Survival Curve - Cluster", i),
                          ggtheme = theme_bw(),
                          palette = cluster_colors_1,
                          risk.table = TRUE,
                          ncensor.plot = FALSE,
                          ncenosr.plot.height = NULL,
                          color = NULL)
  
  # Add the survival plot to the list
  plot_list[[i]] <- surv_plot
}



?ggsurvplot




par(mfrow=c(4,4))

?par
combined_plot <- do.call("par", c(plot_list, mfrow=c(7,4)))

combined_plot

print(plot_list)




DFS_list <- list()

# Iterate through each disease cluster
for (i in unique(combined_data$Disease)) {
  # Subset the data for the current cluster
  i='acc'
  subset_data <- subset(combined_data, Disease == i)
  
  
  # Fit the survival model for the current cluster
  model_fit_DFS <- survfit(Surv(DFS_MONTHS,DFS_binary) ~ cluster, data = subset_data)
  
  # Create the survival plot for the current cluster
  DFS_plot <- ggsurvplot(model_fit_DFS, data = subset_data,
                          pval = TRUE, conf.int = TRUE,risk.table = TRUE, 
                          xlab = "Time (in months)", ylab = "DFS Probability",
                          title = paste("DFS Curve - Disease", i),
                          ggtheme = theme_bw(),
                          palette = cluster_colors_1)
  
  # Add the survival plot to the list
DFS_list[[i]] <- DFS_plot
}


DFS_list






#################################################################################################################################
#################################################################################################################################
#the diseases that were significant to begin with skin cutaneous melanoma, Kidney renal clear cell carcinoma, Kidney renal papillary cell , pancreatic, ovarian, oesophageal

#skcm 

sub_skcm <- subset(combined_data, disease == 'skcm')

skcm_model <- survfit(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_skcm)
skcm_model  
summary(skcm_model)  


plot_list$skcm

skcm_model <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_skcm)
summary(skcm_model)  
survdiff(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_skcm)

tph_skcm <- cox.zph(skcm_model)%>%
  print()


#kirc
plot_list$kirc

survdiff(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_kirc)
sub_kirc <- subset(combined_data, disease == 'kirc')
kirc_model <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_kirc)
summary(kirc_model)

tph_kirc <- cox.zph(kirc_model)%>%
  print()


#kirp
plot_list$kirp

sub_kirp <- subset(combined_data, disease == 'kirp')
kirp_model<- coxph(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_kirp)
summary(kirp_model)

survdiff(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_kirp)

sub_kirp2 <- subset(new_df, disease == 'kirp')
kirp_model2<- coxph(Surv(OS_MONTHS.x, OS_status_binary.x)~cluster+PATH_M_STAGE+PATH_N_STAGE+PATH_T_STAGE, data = sub_kirp2 )
kirp_model2

cox.zph(kirp_model2)

tph_kirp <- cox.zph(kirp_model)%>%
  print()

#paad
plot_list$paad
sub_paad <- subset(combined_data, disease == 'paad')
paad_model<- coxph(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_paad)
summary(paad_model)

survdiff(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_paad)


new_df$OS_status_binary.x<- as.numeric(new_df$OS_status_binary.x)
new_df$OS_MONTHS.x<- as.numeric(new_df$OS_MONTHS.x)

sub_paad2 <- subset(new_df, disease == 'paad')
paad_model2<- coxph(Surv(OS_MONTHS.x, OS_status_binary.x)~cluster+PATH_M_STAGE+PATH_N_STAGE+PATH_T_STAGE, data = sub_paad2)
paad_model2

tph_paad <- cox.zph(paad_model) %>%
  print()

#ov
sub_ov <- subset(combined_data, disease == 'ov')
ov_model<- coxph(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_ov)
summary(ov_model)

tph_ov <- cox.zph(ov_model)%>%
  print()

survdiff(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_ov)

#esca

plot_list$esca
sub_esca <- subset(combined_data, disease == 'esca')
esca_model<- coxph(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_esca)
summary(esca_model)

survdiff(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_esca)


tph_esca <- cox.zph(esca_model)%>%
  print()


#kich

plot_list$kich
sub_kich <- subset(combined_data, disease == 'kich') 
kich_model<- coxph(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_kich)
summary(kich_model)
survdiff(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_kich)

sub_kich2 <- subset(new_df, disease == 'kich')
kich_model2<- coxph(Surv(OS_MONTHS.x, OS_status_binary.x)~cluster+PATH_M_STAGE+PATH_N_STAGE+PATH_T_STAGE, data = sub_kich2)
kich_model2


tph_kich <- cox.zph(kich_model)%>%
  print()


#acc
plot_list$acc
sub_acc <- subset(combined_data, disease == 'acc')
acc_model<- coxph(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_acc)
summary(acc_model)
survdiff(Surv(OS_MONTHS, OS_status_binary)~cluster, data = sub_acc)

plot(cox.zph(acc_model))

tph_acc <- cox.zph(acc_model)%>%
  print()


table(Clara_T3$Disease)
#################################################################################################################################
#################################################################################################################################
#create the forest plots for each of the diseases: 
#ov
ggforest(ov_model, data = sub_ov)


#esca
ggforest(esca_model, data = sub_esca)

#paad
ggforest(paad_model, data = sub_paad)



#################################################################################################################################
#################################################################################################################################
#coxph based on tnm stage 

colnames(combined_data)
colnames(new_df)

select_comb <- subset(combined_data, select = c(1,24,25,31))

new_df <- merge(stacked_data, Total_path, by = 'PATIENT_ID')
new_df <- merge(new_df,select_comb, by = 'PATIENT_ID')

acc_data_clinical_patient<-acc_data_clinical_patient[-c(1), ]
new_df$PATH_T_STAGE<- as.numeric(PATH_T_STAGE)
combined_data$cluster[which(cluster=='#1F77B4')]<-1

attach(new_df)
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T1')]<-1
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T1a')]<-1
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T1a1')]<-1
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T1b')]<-1
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T1b1')]<-1
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T1b2')]<-1
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T1c')]<-1


new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T2')]<-2
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T2a1')]<-2
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T2a2')]<-2
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T2b')]<-2
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T2c')]<-2
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T2a')]<-2



new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T3')]<-3
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T3a')]<-3
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T3b')]<-3
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T3c')]<-3


new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T4')]<-4
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T4a')]<-4
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T4b')]<-4
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'T4d')]<-4

new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'Tis')]<-5
new_df$PATH_T_STAGE[which(PATH_T_STAGE == 'TX')]<-6

table(new_df$PATH_T_STAGE)

new_df$OS_STATUS[which(OS_STATUS =='1:DECEASED')] <- 'DECEASED'
new_df$OS_STATUS[which(OS_STATUS =='0:LIVING')] <- 'LIVING'

new_df$OS_status_binary <- ifelse(new_df$OS_STATUS == "DECEASED", 1,
                                         ifelse(new_df$OS_STATUS == "LIVING", 0,
                                                NA))


new_df$OS_MONTHS <- as.numeric(new_df$OS_MONTHS)



TNM_model <- coxph(Surv(OS_MONTHS,OS_status_binary)~as.factor(PATH_T_STAGE), data = new_df)
summary(TNM_model)

ggsurvplot(survfit(TNM_model), data = new_df, risk.table = TRUE, xlab = "Time", ylab = "Survival Probability",
           title = "Survival Curves by PATH_T_STAGE", legend.title = "PATH_T_STAGE")


#################################################################################################################################
#################################################################################################################################
table(new_df$PATH_N_STAGE)



new_df$PATH_N_STAGE <- ifelse(new_df$PATH_N_STAGE == 'N0', 1,
                              ifelse(new_df$PATH_N_STAGE == 'N0 (i-)', 1,
                                     ifelse(new_df$PATH_N_STAGE == 'N0 (i+)',1,
                                            ifelse(new_df$PATH_N_STAGE == 'N0 (mol+)',1, 
                                                   ifelse(new_df$PATH_N_STAGE == 'N1',2,
                                                          ifelse(new_df$PATH_N_STAGE == 'N1a',2,
                                                                 ifelse(new_df$PATH_N_STAGE == 'N1b',2,
                                                                        ifelse(new_df$PATH_N_STAGE == 'N1c',2,
                                                                               ifelse(new_df$PATH_N_STAGE == 'N1mi',2,
                                                                                      ifelse(new_df$PATH_N_STAGE == 'N2', 3,
                                                                                             ifelse(new_df$PATH_N_STAGE == 'N2a',3,
                                                                                                    ifelse(new_df$PATH_N_STAGE == 'N2b',3,
                                                                                                           ifelse(new_df$PATH_N_STAGE== 'N2c',3,
                                                                                                                  ifelse(new_df$PATH_N_STAGE == 'N3', 4,
                                                                                                                         ifelse(new_df$PATH_N_STAGE == 'N3a',4,
                                                                                                                                ifelse(new_df$PATH_N_STAGE == 'N3b',4,
                                                                                                                                       ifelse(new_df$PATH_N_STAGE == 'N3c',4,
                                                                                                                                              ifelse(new_df$PATH_N_STAGE == 'NX',5,
                                                                                                                                                     NA))))))))))))))))))

table(new_df$PATH_N_STAGE)

3724+154+28+1
1116+297+168+5+37


table(new_df$PATH_M_STAGE)

new_df$PATH_M_STAGE <- ifelse(new_df$PATH_M_STAGE == 'cM0 (i+)',1,
                              ifelse(new_df$PATH_M_STAGE == 'M0',1,
                                     ifelse(new_df$PATH_M_STAGE == 'M1',2,
                                            ifelse(new_df$PATH_M_STAGE == 'M1a',2,
                                                   ifelse(new_df$PATH_M_STAGE == 'M1b',2,
                                                          ifelse(new_df$PATH_M_STAGE == 'M1c',2,
                                                                 ifelse(new_df$PATH_M_STAGE == 'MX',3,
                                                                        NA)))))))


new_df$PATH_M_STAGE <- as.factor(new_df$PATH_M_STAGE)
new_df$PATH_N_STAGE <- as.factor(new_df$PATH_N_STAGE)                                                                 
new_df$PATH_T_STAGE <- as.factor(new_df$PATH_T_STAGE)



TNM_model <- coxph(Surv(OS_MONTHS,OS_status_binary)~PATH_T_STAGE+PATH_M_STAGE+PATH_N_STAGE, data = new_df)
summary(TNM_model)
plot(TNM_model)

TNM_model2 <- coxph(Surv(OS_MONTHS,OS_status_binary)~PATH_M_STAGE + PATH_N_STAGE, data = new_df)
summary(TNM_model2)


ggsurvplot(survfit(TNM_model2), data = new_df, risk.table = TRUE, xlab = "Time", ylab = "Survival Probability",
           title = "Survival Curves by PATH_T_STAGE", legend.title = "N/M")


#################################################################################################################################
#SKCM censored 

sub_skcm <- subset(combined_data, disease == 'skcm')

subset_skcm <- subset(sub_skcm, OS_MONTHS <= 100)


skcm_model <- survfit(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = subset_skcm)
skcm_model  
summary(skcm_model)  

ggsurvplot(skcm_model, data = subset_skcm, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, risk.table.title = "Risk table",
           xlab = "Time (in months)", ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve with Events Marked",
           tables.theme = theme_cleantable(),
           surv.median.line = "hv",
           break.time.by = 12,
           ggtheme = theme_bw(),
           palette = my_colors)


plot_list$skcm

skcm_model <- coxph(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = subset_skcm)
summary(skcm_model)  
survdiff(Surv(OS_MONTHS,OS_status_binary) ~ cluster, data = sub_skcm)

tph_skcm <- cox.zph(skcm_model)%>%
  print()



########
#########################################################################################################################



