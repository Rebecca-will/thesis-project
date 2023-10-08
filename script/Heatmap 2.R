set.seed(100)

set_1 <- Clara_T[,1:3]
set_2 <- Clara_T[,7:8]
set_3 <- Clara_T[,14:30]

Clara_T2 <- cbind(set_1,set_2,set_3)

colnames(Clara_T2)


Clara_T3 <- Clara_T2[!duplicated(Clara_T2$PATIENT_ID),]
dim(Clara_T3)#use this as it doesn't have the duplicated pts 


#scale data 
scaled_data_3<- apply(Clara_T3[,3:22], 2, scale)

#distances

Row_dist_2 <- dist(t(scaled_data_2), method='euclidean')

Col_dist_2 <- dist(scaled_data_2, method='euclidean')

#Clustering 
Hr_row_2 <- hclust(Row_dist_2, method= 'ward.D')
dd_row_2 <- as.dendrogram(Hr_row_2)

HR_col_2 <- hclust(Col_dist_2, method='ward.D')
dd_col_2 <- as.dendrogram(HR_col_2)

hclust(dist(scaled_data_2, method='euclidean'),method = 'ward.D')


#calculating the gap statistic 

# myclusGapcols_2 <- clusGap(scaled_data_2,
#                          
#                          FUN = mycluster,
#                          
#                          K.max = 15,
#                          
#                          B = 500)
# 
# 
# myclusGaprows_2 <- clusGap(t(scaled_data_2),
#                          
#                          FUN = mycluster,
#                          
#                          K.max = 15,
#                          
#                          B = 500)
#number of clusters:1


krows_2 <- maxSE(myclusGaprows_2$Tab[, "gap"], myclusGaprows_2$Tab[, "SE.sim"], method="globalSEmax")
krows_2<- 14
#says optimum is 14, possibly due to different cluster method
#signatures

?maxSE
kcols_2 <- maxSE(myclusGapcols_2$Tab[, "gap"], myclusGapcols_2$Tab[, "SE.sim"], method="globalSEmax")
kcols_2<-15
#says 15 clusters


colors=rainbow(15)#keep this in mind to the k.max, colour to each cluster


cutclhr_2 <- cutree(Hr_row_2, k = krows_2)#cut the tree at many clusters as gap has decided 

cutclcolhr_2 <-my_colors[1:krows_2]#assigning colours to samples to the cluster

cutclcolhr_2 <- cutclcolhr_2[as.vector(cutclhr_2)]


#samples
cutclhc_2 <- cutree(HR_col_2, k = kcols_2)

cutclcolhc_2 <-my_colors[1:kcols_2]

cutclcolhc_2 <- cutclcolhc_2[as.vector(cutclhc_2)]

#removing the upper and lower 2%

#quantile(mgeneexpressionnorm$value, 0.98)
mgeneexpressionnorm_2 <- melt(scaled_data_2[,!(colnames(scaled_data_2) %in% c("Sample", "Disease"))])
#find the 2% value

TCGAClaraTsigsheatmap_2 <- scaled_data_2[, !(colnames(scaled_data_2) %in% c("Sample", "Disease"))]

TCGAClaraTsigsheatmap_2[TCGAClaraTsigsheatmap_2 > quantile(mgeneexpressionnorm_2$value, 0.98)] <- quantile(mgeneexpressionnorm_2$value, 0.98)

TCGAClaraTsigsheatmap_2[TCGAClaraTsigsheatmap_2 < quantile(mgeneexpressionnorm_2$value, 0.02)] <- quantile(mgeneexpressionnorm_2$value, 0.02)






unique_diseases <- levels(as.factor(Clara_T$Disease))
unique_diseases

# Print disease-color pairs
for (i in 1:length(diseasecolors_set)) {
  print(paste(diseasecolors_set[i], "is", my_colors[i]))
}



diseasecolors_set


#heatmap generation
plot_4 <- heatmap3(as.matrix(t(TCGAClaraTsigsheatmap_2)), scale = "none", labCol = "",
                   Rowv = dd_row_2,RowSideColors = cutclcolhr_2, Colv = rev(dd_col_2),ColSideColors = cbind(diseasecolors_set,cutclcolhc_2),RowSideLabs = 'Signature Clusters', ColSideLabs = cbind('Disease' , 'Sample Clusters'),
                   col = greenred(500), balanceColor = TRUE, cexCol = 1.5, cexRow = 1.0, margins = c(10, 20))


#Have adjusted and looks okay on the zoomed in display
legend('left', legend=unique(as.factor(Clara_T$Disease)),
       
       fill=unique(diseasecolors_set), border=FALSE, bty="n", cex=0.5)

dev.off()

?legend



#heatmap with three clusters for coloumns and 2 for rows 
krows_3 <- 4

cutclhr_3 <- cutree(Hr_row_2, k = krows_3)#cut the tree at many clusters as gap has decided 

cutclcolhr_3 <-my_colors[1:krows_3]#assigning colours to samples to the cluster

cutclcolhr_3 <- cutclcolhr_3[as.vector(cutclhr_3)]


kcols_3<- 3

cutclhc_3 <- cutree(HR_col_2, k = kcols_3)

cutclcolhc_3 <-my_colors[1:kcols_3]

cutclcolhc_3 <- cutclcolhc_3[as.vector(cutclhc_3)]


plot_5 <- heatmap3(as.matrix(t(TCGAClaraTsigsheatmap_2)), scale = "none", labCol = "",
                   Rowv = dd_row_2,RowSideColors = cutclcolhr_3, Colv = rev(dd_col_2),ColSideColors = cbind(diseasecolors_set,cutclcolhc_3),RowSideLabs = 'Signature Clusters', ColSideLabs = cbind('Disease' , 'Sample Clusters'),
                   col = greenred(500), balanceColor = TRUE, cexCol = 1.5, cexRow = 1.0, margins = c(10, 20))


#Have adjusted and looks okay on the zoomed in display
legend('left', legend=unique(as.factor(Clara_T$Disease)),
       
       fill=unique(diseasecolors_set), border=FALSE, bty="n", cex=0.5)

unique(as.factor(Clara_T$Disease))

#percentage in each cluster 
#blue is cluster 1 1F77B4
#orange is cluster 3 FF7F0E
#green is cluster 2 2CA02C

print(cutclcolhc_3)
library(dplyr)
?group_by

table(Clara_T$Disease) #total number of each disease 

table(cutclcolhc_3)# how many of each sample is separated into each cluster, is this the correct way to tell the order?

Disease_df <- cbind(Clara_T$Disease, cutclcolhc_3)

?subset
Cluster_1 <- subset(Disease_df, cutclcolhc_3 == '#1F77B4')
table(Cluster_1)
dim(Cluster_1)

Cluster_2 <- subset(Disease_df, cutclcolhc_3 == '#2CA02C')
table(Cluster_2)
dim(Cluster_2)

Cluster_3 <- subset(Disease_df, cutclcolhc_3 == '#FF7F0E')
table(Cluster_3)


# count_percentage<- function(df) {
#   input <-readline(prompt = 'Which disease: ')
#   count <- sum(df == input)
#   ans<-(count/nrow(df))*100
#   print(ans)
# }
# 
# count_percentage(Cluster_3)  
# 
#  sum(Cluster_1 == 'gbm')
# print(count_1)


#enrichment analysis for each cluster 
counts <- table(Cluster_1[,1])
percentages_culster_1 <- prop.table(counts) * 100

percentages_culster_1


for (val in names(percentages_culster_1)) {
  message(paste(val, "is", round(percentages_culster_1[val], 2), "% of cluster 1"))
}

#cluster 2
counts_2 <- table(Cluster_2[,1])
percentages_culster_2 <- prop.table(counts_2,) * 100


 for (val in names(percentages_culster_2)) {
  message(paste(val, "is", round(percentages_culster_2[val], 2), "% of cluster 2"))
}

#cluster 3
counts_3 <- table(Cluster_3[,1])
percentages_culster_3 <- prop.table(counts_3) * 100


for (val in names(percentages_culster_3)) {
  message(paste(val, "is", round(percentages_culster_3[val], 2), "% of cluster 3"))
}

unique(Clara_T$Disease)

#disease distribution across clusters

counts_4 <- table(Disease_df)
percentages_df <- prop.table(counts_4)*100 

for (val in names(percentages_df)) {
  message(paste(val, "is", round(percentages_df[val], 2), "% of df"))
}

 count_percentage<- function(x) {
   input <-as.character(x)
   cul_1 <- sum(Cluster_1 == input)
   cul_2 <- sum(Cluster_2 == input)
   cul_3 <- sum(Cluster_3 == input)
   ans_1<-(cul_1/sum(Clara_T== input)*100)
   ans_2<-(cul_2/sum(Clara_T== input)*100)
   ans_3<-(cul_3/sum(Clara_T== input)*100)
   print(paste(x,'Cluster 1: ',round(ans_1,2), 'Cluster_2: ', round(ans_2,2), 'Cluster_3: ', round(ans_3,2)))
 }

 
 count_percentage('brca')

 result <- lapply(unique(Clara_T$Disease), count_percentage)
 
 #summary of the enrichment data 
 
 for (val in names(percentages_culster_1)) {
   message(paste(val, "is", round(percentages_culster_1[val], 2), "% of cluster 1"))
 }
 
 for (val in names(percentages_culster_2)) {
   message(paste(val, "is", round(percentages_culster_1[val], 2), "% of cluster 2"))
 }
 
 for (val in names(percentages_culster_3)) {
   message(paste(val, "is", round(percentages_culster_1[val], 2), "% of cluster 3"))
 }
 
 result
 
 
 #represented as pie charts 
 library(ggplot2)
 
 data_frame_cluster_2 <- data.frame(Disease = names(counts_2), Count = counts_2, Percentage = percentages_culster_2)
 
 
 pie_chart_cluster_2 <- ggplot(data_frame_cluster_2, aes(x = "", y = Percentage.Freq, fill = Disease)) +
   geom_bar(stat = "identity", width = 1) +
   coord_polar("y", start = 0) +
   scale_fill_manual(values =  data_frame_cluster_2$color) + 
   labs(title = "Disease Distribution in Cluster 2")
 
 print(pie_chart_cluster_2)
 

 data_frame_cluster_2 <- merge(data_frame_cluster_2, color_df, by.x = 'Disease', by.y = 'disease')
dat




 #data frame, r relies on indexes being the same order. 
disease_list <- levels(as.factor(unique(Clara_T$Disease)))

hex_codes <- as.character(unique(diseasecolors_set))

color_df <- data.frame(disease = unique(as.factor(Clara_T$Disease)), color = unique(diseasecolors_set))

color_df <-color_df[order(color_df$disease),]
color_df


 #survival analysis
 
#Clara_T<- subset(Clara_T, select = -short_name)

Clara_T2$PATIENT_ID <- substr(Clara_T2$Sample, 1, 12)

combined_data <- merge(Clara_T3, stacked_data, by = 'PATIENT_ID')

dim(combined_data)
dim(Clara_T3)

complete.cases(combined_data$PATIENT_ID)%>%
  sum()

complete.cases(stacked_data$PATIENT_ID)%>%
  sum()
dim(stacked_data)

complete.cases(Clara_T3$Sample)%>%
  sum()

duplicated(Clara_T3$Sample)%>%
  sum()


dup_df<- subset(Clara_T3,  PATIENT_ID %in% missing_pt)

duplicated(stacked_data$PATIENT_ID)%>%
  sum()

duplicated(combined_data$PATIENT_ID)%>%
  sum()


duplicated_patient<- Clara_T$PATIENT_ID[duplicated(Clara_T$PATIENT_ID)]
  
duplicated_data <- Clara_T[Clara_T$PATIENT_ID %in% duplicated_patient,]



#claraT complete cases- rows of combined 
8567-8610

duplicated_patient
dup<- stacked_data[stacked_data$PATIENT_ID %in% Clara_T3$PATIENT_ID, ]
print(dup)
dim(dup)



missing_pt<- setdiff(Clara_T3$PATIENT_ID, stacked_data$PATIENT_ID) 


intersect(duplicated_patient, stacked_data$PATIENT_ID) %>%
  print()

setdiff(missing_pt, Clara_T3$PATIENT_ID)%>%
  print()

setdiff(missing_pt, stacked_data$PATIENT_ID)%>%
  print()#43 patients in the stacked data that is not in ClaraT

setdiff(duplicated_patient, Clara_T3$PATIENT_ID)%>%
  print()#Clara_T3 can be used 


#repeated with non duplicated pts 
write.csv(Clara_T3, file = "Clara_T3.csv", row.names = FALSE)
dim(Clara_T3)

dim(Clara_T2)
