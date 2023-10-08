dim(Clara_T3)
set.seed(100)


#scale the data 
scaled_data<- apply(Clara_T3[,3:22], 2, scale)
scaled_data


#distances

Row_dist <- dist(t(scaled_data), method='euclidean')

Col_dist <- dist(scaled_data, method='euclidean')

#Clustering - heirar
Hr_row <- hclust(Row_dist, method= 'ward.D')

dd_row <- as.dendrogram(Hr_row)

HR_col <- hclust(Col_dist, method='ward.D')

dd_col <- as.dendrogram(HR_col)

hclust(dist(scaled_data, method='euclidean'),method = 'ward.D')

#heatmap
library(heatmap3)
library(gplots)
library(reshape2)
library(cluster)

#removing the upper and lower 2%

#quantile(mgeneexpressionnorm$value, 0.98)
mgeneexpressionnorm <- melt(scaled_data[,!(colnames(scaled_data) %in% c("Sample", "Disease"))])
#find the 2% value

TCGAClaraTsigsheatmap <- scaled_data[, !(colnames(scaled_data) %in% c("Sample", "Disease"))]

TCGAClaraTsigsheatmap[TCGAClaraTsigsheatmap > quantile(mgeneexpressionnorm$value, 0.98)] <- quantile(mgeneexpressionnorm$value, 0.98)

TCGAClaraTsigsheatmap[TCGAClaraTsigsheatmap < quantile(mgeneexpressionnorm$value, 0.02)] <- quantile(mgeneexpressionnorm$value, 0.02)




mycluster <- function(x, k) list(cluster= cutree(hclust(dist(x, method='euclidean'),method = 'ward.D'), k = k))




# Calculate the gap statistic for up to 6 clusters
#change k.max to 5-> to make it more interpretable with the survival 
# myclusGapcols <- clusGap(scaled_data,
#                          
#                          FUN = mycluster,
#                          
#                          K.max = 15,
#                          
#                          B = 500)
# #output:15
#cluster in the same way as before and tell how similar the clusters are to eachother and what is the optimal number of clusters


# #signature clusters = 2
# 
 myclusGaprows <- clusGap(t(scaled_data),
                         
                          FUN = mycluster,
                          
                          K.max = 15,
                          
                          B = 500)
#output:13


#Get number of clusters from gap statistic, maxSE is objective, gives a solid reason

# krows <- maxSE(myclusGaprows$Tab[, "gap"], myclusGaprows$Tab[, "SE.sim"], method="globalSEmax")
# krows <- 13
# #signatures
# 
# ?maxSE
# kcols <- maxSE(myclusGapcols$Tab[, "gap"], myclusGapcols$Tab[, "SE.sim"], method="globalSEmax")
# kcols <- 15
#number of clusters for your samples 


my_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
               "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
               "#1A1A1A", "#FF4A46", "#008B8B", "#A0522D", "#800080",
               "#FFD700", "#D2691E", "#228B22", "#6B8E23", "#2F4F4F",
               "#F08080", "#00FA9A", "#FF1493")

krows <- 4

kcols <-3


cutclhr <- cutree(Hr_row, k = krows)#cut the tree at many clusters as gap has decided 

cutclcolhr <-my_colors[1:krows]#assigning colours to samples to the cluster

cutclcolhr <- cutclcolhr[as.vector(cutclhr)]



#samples
cutclhc <- cutree(HR_col, k = kcols)

cutclcolhc <-my_colors[1:kcols]

cutclcolhc <- cutclcolhc[as.vector(cutclhc)]


#disease 
diseasecolors = my_colors #change to different pallette later 

diseasecolors_set <- my_colors[as.factor(Clara_T3$Disease)] #changed to as.factor since its characters



plot <- heatmap3(as.matrix(t(TCGAClaraTsigsheatmap)), scale = "none", labCol = "",
                   Rowv = dd_row,RowSideColors = cutclcolhr, Colv = rev(dd_col),ColSideColors = cbind(diseasecolors_set,cutclcolhc),RowSideLabs = 'Signature Clusters', ColSideLabs = cbind('Disease' , 'Sample Clusters'),
                   col = greenred(500), balanceColor = TRUE, cexCol = 1.5, cexRow = 1.0, margins = c(7.5, 20))


#Have adjusted and looks okay on the zoomed in display
legend('right', legend=unique(as.factor(Clara_T3$Disease)),
       
       fill=unique(diseasecolors_set), border=FALSE, bty="n", cex=1)


?legend



#################################################################################################################################
#################################################################################################################################
# trying to jpeg it 
install.packages("Cairo")
library(Cairo)


CairoJPEG("heatmap.jpg", width = 10, height = 8, res = 600)
par(mar = c(2, 2, 2, 2))

plot <- heatmap3(as.matrix(t(TCGAClaraTsigsheatmap)), scale = "none", labCol = "",
                 Rowv = dd_row,RowSideColors = cutclcolhr, Colv = rev(dd_col),ColSideColors = cbind(diseasecolors_set,cutclcolhc),RowSideLabs = 'Signature Clusters', ColSideLabs = cbind('Disease' , 'Sample Clusters'),
                 col = greenred(500), balanceColor = TRUE, cexCol = 1.5, cexRow = 1.0, margins = c(2, 2))


#Have adjusted and looks okay on the zoomed in display
legend(x = 0.02, y = 1, legend=unique(as.factor(Clara_T3$Disease)),
       
       fill=unique(my_colors), border=FALSE, bty="n", cex=0.5)

dev.off()

#still cant get this to work :(

#################################################################################################################################
#################################################################################################################################
#enrichment of diseases in new clusters


table(cutclcolhc)
#1F77B4:blue is cluster 1
#2CA02C:green is cluster 2
#FF7F0E:orange is cluster 3

Disease_df <- cbind(Clara_T3$Disease, cutclcolhc)

Cluster_1 <- subset(Disease_df, cutclcolhc == '#1F77B4')
table(Cluster_1)
dim(Cluster_1)


Cluster_2 <- subset(Disease_df, cutclcolhc == '#2CA02C')
table(Cluster_2)
dim(Cluster_2)

Cluster_3 <- subset(Disease_df, cutclcolhc == '#FF7F0E')
table(Cluster_3)

#################################################################################################################################
#################################################################################################################################
#enrichment analysis for disease within each cluster 
counts <- table(Cluster_1[,1])
percentages_culster_1 <- prop.table(counts,) * 100

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
#################################################################################################################################
#################################################################################################################################
#% of diseases across clusters

count_percentage<- function(x) {
  input <-as.character(x)
  cul_1 <- sum(Cluster_1 == input)
  cul_2 <- sum(Cluster_2 == input)
  cul_3 <- sum(Cluster_3 == input)
  ans_1<-(cul_1/sum(Clara_T3== input)*100)
  ans_2<-(cul_2/sum(Clara_T3== input)*100)
  ans_3<-(cul_3/sum(Clara_T3== input)*100)
  print(paste(x,'Cluster 1: ',round(ans_1,2), 'Cluster_2: ', round(ans_2,2), 'Cluster_3: ', round(ans_3,2)))
}


count_percentage('brca')

result <- lapply(unique(Clara_T3$Disease), count_percentage)

#create df for clusters
df_cluster_1 <- data.frame(Disease = names(counts), Count = counts, Percentage = percentages_culster_1)
df_cluster_2 <- data.frame(Disease = names(counts_2), Count = counts_2, Percentage = percentages_culster_2)
df_cluster_3 <- data.frame(Disease = names(counts_3), Count = counts_3, Percentage = percentages_culster_3)


#creating the colour list to match 
disease_list <- levels(as.factor(unique(Clara_T3$Disease)))
hex_codes <- as.character(unique(diseasecolors_set))
color_df <- data.frame(disease = unique(as.factor(Clara_T3$Disease)), color = unique(diseasecolors_set))
color_df <-color_df[order(color_df$disease),]
color_df

df_cluster_1 <- merge(df_cluster_1, color_df, by.x = 'Disease', by.y = 'disease')
df_cluster_2 <- merge(df_cluster_2, color_df, by.x = 'Disease', by.y = 'disease')
df_cluster_3 <- merge(df_cluster_3, color_df, by.x = 'Disease', by.y = 'disease')


#################################################################################################################################
#################################################################################################################################
#pie charts

pie_chart_cluster_1 <- ggplot(df_cluster_1, aes(x = "", y = Percentage.Freq, fill = Disease)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values =  df_cluster_1$color) + 
  labs(title = "Disease Distribution in Cluster 1")

print(pie_chart_cluster_1)
ggsave(filename = "Disease Dis 1.jpg", plot = pie_chart_cluster_1, dpi = 600)

pie_chart_cluster_2 <- ggplot(df_cluster_2, aes(x = "", y = Percentage.Freq, fill = Disease)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values =  df_cluster_2$color) + 
  labs(title = "Disease Distribution in Cluster 2")

print(pie_chart_cluster_2)
ggsave(filename = "Disease Dis 2.jpg", plot = pie_chart_cluster_2, dpi = 600)

pie_chart_cluster_3 <- ggplot(df_cluster_3, aes(x = "", y = Percentage.Freq, fill = Disease)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values =  df_cluster_3$color) + 
  labs(title = "Disease Distribution in Cluster 3")

print(pie_chart_cluster_3)
ggsave(filename = "Disease Dis 3.jpg", plot = pie_chart_cluster_3, dpi = 600)

#################################################################################################################################
#################################################################################################################################
#combine data and check they are similar if not why 
# the duplicates have already been removed and the patient ID created in the Clara_T3 dataset from previous 
Clara_T3$cluster <- cutclcolhc
Clara_T3$cluster_colour <- cutclcolhc

combined_data <- merge(Clara_T3, stacked_data[2:7], by = 'PATIENT_ID')

dim(combined_data)
dim(Clara_T3)

#ClaraT3- comb
8610-8567


missing_pt<- setdiff(Clara_T3$PATIENT_ID, stacked_data$PATIENT_ID) 
length(missing_pt)

setdiff(missing_pt, stacked_data$PATIENT_ID)%>%
  print()#43 patients in the stacked data that is not in ClaraT

setdiff(missing_pt, Clara_T3$PATIENT_ID)%>%
  print()#Clara_T3 can be used, we have the RNA but not the survival data. 


table(cutclcolhc)
