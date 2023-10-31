library(survival)
library(survminer)
library(ggplot2)

#Adding a column to the data frame that identifies the related cluster and a second to identify the colour later
cutclcolhc <- readLines('cutclcolhc.txt')

Clara_T3$cluster <- cutclcolhc
Clara_T3$cluster_colour <- cutclcolhc

#merging the survival data contained in the stacked_data df(obtained from TCGA) and the RNA expression data by the patient ID 
combined_data <- merge(Clara_T3, stacked_data[2:7], by = 'PATIENT_ID')

#Check to ensure all data was matched
dim(combined_data)
dim(Clara_T3)

#ClaraT3- comb
8610-8567

#identifying where the missing patients are 
missing_pt<- setdiff(Clara_T3$PATIENT_ID, stacked_data$PATIENT_ID) 
length(missing_pt)

setdiff(missing_pt, stacked_data$PATIENT_ID)%>%
  print()#43 patients in the stacked data that is not in ClaraT

setdiff(missing_pt, Clara_T3$PATIENT_ID)%>%
  print()#Clara_T3 can be used, we have the RNA but not the survival data. 


#cleaning the data 
table(combined_data$OS_STATUS) #investigating th contents of the column 

#Reassigning values to consistent names
combined_data$OS_STATUS[which(combined_data$OS_STATUS =='1:DECEASED')] <- 'DECEASED'
combined_data$OS_STATUS[which(combined_data$OS_STATUS =='0:LIVING')] <- 'LIVING'
combined_data$OS_STATUS[which(combined_data$OS_STATUS =='[Not Available]')] <- NA

#creating a new column with a binary representation of survival 
combined_data$OS_status_binary <- ifelse(combined_data$OS_STATUS == "DECEASED", 1,
                                         ifelse(combined_data$OS_STATUS == "LIVING", 0,
                                                NA))

#converting value to recognisable NA and ensuring the values are seen as numeric
combined_data$OS_MONTHS[which(combined_data$OS_MONTHS =='[Not Available]')] <- NA
combined_data$OS_MONTHS <- as.numeric(combined_data$OS_MONTHS)

#changing the column with the hex cluster code to a number for ease of interpretation 
combined_data$cluster <- ifelse(combined_data$cluster == '#1F77B4',1,
                                ifelse(combined_data$cluster == '#2CA02C',2,
                                       3))
class(combined_data$cluster)
combined_data$cluster <- as.factor(combined_data$cluster) #we want the clusters to be recognised as groups and not numbers

#survival analysis using cox regression to see the effect of cluster while controlling for disease on survival
cox_surv <- coxph(Surv(OS_MONTHS,OS_status_binary)~cluster + disease, data = combined_data)
summary(cox_surv)
cox_surv

#checking the proportional hazards of the cox model- hazards must be constant over time
main_tph <- cox.zph(cox_surv) %>%
  print()

plot(cox.zph(cox_surv))

#Basic forest plot of the hazards of the variables in the model
ggforest(cox_surv, data = combined_data)

#creating the Kaplan Meier survival curves 
Surv_model<- survfit(Surv(OS_MONTHS,OS_status_binary)~cluster, data = combined_data)

surv_plot2<- ggsurvplot(Surv_plot_1, data = combined_data, pval = TRUE, pval.method= TRUE,conf.int = TRUE,
                        risk.table = TRUE, risk.table.title = "Risk table",
                        xlab = "Time (in months)", ylab = "Survival Probability",
                        title = "Kaplan-Meier Survival Curve with Events Marked",
                        tables.theme = theme_cleantable(),
                        surv.median.line = "hv",
                        break.time.by = 12,
                        ggtheme = theme_bw(),
                        palette = cluster_colors_1,
                        ncensor.plot = TRUE,  
                        ncensor.plot.height = 0.25)
surv_plot2

#Creating a loop to assess the individual diseases survivals and plot results 

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

plot_list <- list()

print(plot_list)



table(combined_data$DFS_STATUS)
#reasign the different variables to one 
combined_data$DFS_STATUS[which(combined_data$DFS_STATUS =='1:Recurred/Progressed')] <- 'Recurred/Progressed'
combined_data$DFS_STATUS[which(combined_data$DFS_STATUS =='0:DiseaseFree')] <- 'DiseaseFree'
combined_data$DFS_STATUS[which(combined_data$DFS_STATUS =='[Not Available]')] <- NA

table(combined_data$DFS_STATUS)


combined_data$DFS_binary <- ifelse(combined_data$DFS_STATUS == "Recurred/Progressed", 1,
                                   ifelse(combined_data$DFS_STATUS == "DiseaseFree", 0,
                                          NA))


table(combined_data$DFS_binary)
class(combined_data$DFS_MONTHS)
class(combined_data$DFS_binary)

combined_data$DFS_binary<- as.factor(combined_data$DFS_binary)
combined_data$DFS_MONTHS<- as.numeric(combined_data$DFS_MONTHS)


#fitting the survfit model to carry out the Kaplan Meier plot 
model_fit_DFS <- survfit(Surv(DFS_MONTHS,DFS_binary)~cluster, data = combined_data)

#creating the Plot of DFS between the three clusters
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

#Creating a cox model to tests the differences in hazards between the clusters 
cox_DFS <- coxph(Surv(DFS_MONTHS,DFS_binary)~cluster, data = combined_data)
summary(cox_DFS)

#Testing that the hazards remain constant 
HRDFS_tph <- cox.zph(cox_DFS) %>%
  print()

plot(cox.zph(cox_DFS))


#Similar loop created but to assess disease free survival 
DFS_list <- list()

# Iterate through each disease cluster
for (i in unique(combined_data$Disease)) {
  # Subset the data for the current cluster
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


#Creating a new df with just the columns of interest in TNM analysis. 
select_columns <- subset(combined_data, select = c(1,24:32))

#combining the reduced df with the df containing the TNM data using patient ID column.
Surv_TNM_df <- merge(select_columns, Total_path, by = 'PATIENT_ID')


#Cleaning the survival data from TCGA for consistent values
vals_to_replace_1 <- c('T1', 'T1a', 'T1a1', 'T1b', 'T1b1', 'T1b2', 'T1c','N1','N1a','N1b','N1c','N1mi','M1','M1a','M1b','M1c')
vals_to_replace_2 <- c('T2', 'T2a1', 'T2a2', 'T2b', 'T2c', 'T2a','N2','N2a', 'N2b', 'N2c')
vals_to_replace_3 <- c('T3', 'T3a', 'T3b', 'T3c','N3', 'N3a', 'N3b','N3c')
vals_to_replace_4 <- c('T4', 'T4a', 'T4b', 'T4d')
vals_to_replace_0 <- c('T0','N0','N0 (i-)', 'N0 (i+)','N0 (mol+)', 'cM0 (i+)', 'M0','MX')
vals_to_replace_x <- c('TX', 'NX')
other_vals_to_replace <- c('[Discrepancy]','[Not Available]')


cols_to_replace <- c('PATH_T_STAGE', 'PATH_N_STAGE', 'PATH_M_STAGE')

#Using an anonymous function to replace the values in the column 
Surv_TNM_df[cols_to_replace] <- sapply(Surv_TNM_df[cols_to_replace],
                                      function(x) replace(x, x %in% vals_to_replace_1, 1)) 

Surv_TNM_df[cols_to_replace] <- sapply(Surv_TNM_df[cols_to_replace],
                                      function(x) replace(x, x %in% vals_to_replace_2, 2)) 

Surv_TNM_df[cols_to_replace] <- sapply(Surv_TNM_df[cols_to_replace],
                                      function(x) replace(x, x %in% vals_to_replace_3, 3)) 

Surv_TNM_df[cols_to_replace] <- sapply(Surv_TNM_df[cols_to_replace],
                                      function(x) replace(x, x %in% vals_to_replace_4, 4)) 

Surv_TNM_df[cols_to_replace] <- sapply(Surv_TNM_df[cols_to_replace],
                                       function(x) replace(x, x %in% vals_to_replace_0, 5)) 

Surv_TNM_df[cols_to_replace] <- sapply(Surv_TNM_df[cols_to_replace],
                                       function(x) replace(x, x %in% vals_to_replace_x, 6)) 

Surv_TNM_df[cols_to_replace] <- sapply(Surv_TNM_df[cols_to_replace],
                                      function(x) replace(x, x %in% other_vals_to_replace, NA)) 


Surv_TNM_df$PATH_T_STAGE[which(Surv_TNM_df$PATH_T_STAGE == 'Tis')]<-7


TStage_TNM_model <- coxph(Surv(OS_MONTHS,OS_status_binary)~as.factor(PATH_T_STAGE) + cluster, data = Surv_TNM_df)
summary(TStage_TNM_model)

NStage_TNM_model <- coxph(Surv(OS_MONTHS,OS_status_binary)~as.factor(PATH_N_STAGE) + cluster, data = Surv_TNM_df)
summary(NStage_TNM_model)

MStage_TNM_model <- coxph(Surv(OS_MONTHS,OS_status_binary)~as.factor(PATH_M_STAGE) + cluster, data = Surv_TNM_df)
summary(MStage_TNM_model)
