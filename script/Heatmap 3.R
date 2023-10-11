#packages needed
library(heatmap3)
library(gplots)
library(reshape2)
library(cluster)
dim(Clara_T3)
set.seed(100)

#scale the data 
scaled_data<- apply(Clara_T3[,3:22], 2, scale)
scaled_data


#Calculating the distances for clustering (the dissimilarity between rows and columns respectively)
Row_dist <- dist(t(scaled_data), method='euclidean')
Col_dist <- dist(scaled_data, method='euclidean')

#Hierarchcal clustering using wards linkage to minimise variance
Hr_row <- hclust(Row_dist, method= 'ward.D')
dd_row <- as.dendrogram(Hr_row)

HR_col <- hclust(Col_dist, method='ward.D')
dd_col <- as.dendrogram(HR_col)

hclust(dist(scaled_data, method='euclidean'),method = 'ward.D')


#removing the upper and lower 2% for displaying the heat map
#quantile(mgeneexpressionnorm$value, 0.98)
mgeneexpressionnorm <- melt(scaled_data[,!(colnames(scaled_data) %in% c("Sample", "Disease"))])
TCGAClaraTsigsheatmap <- scaled_data[, !(colnames(scaled_data) %in% c("Sample", "Disease"))]

#upper 2%
TCGAClaraTsigsheatmap[TCGAClaraTsigsheatmap > quantile(mgeneexpressionnorm$value, 0.98)] <- quantile(mgeneexpressionnorm$value, 0.98)
#lower2%
TCGAClaraTsigsheatmap[TCGAClaraTsigsheatmap < quantile(mgeneexpressionnorm$value, 0.02)] <- quantile(mgeneexpressionnorm$value, 0.02)



#creating a function to use to identify the optimal number of clusters
mycluster <- function(x, k) list(cluster= cutree(hclust(dist(x, method='euclidean'),method = 'ward.D'), k = k))

#Calculate the gap statistic for up to 15 clusters as per kmax, b= bootstrap value.
myclusGapcols <- clusGap(scaled_data,

                         FUN = mycluster,

                         K.max = 15,

                         B = 5)
#output:15

# #signature clusters = 2
 myclusGaprows <- clusGap(t(scaled_data),
                         
                          FUN = mycluster,
                          
                          K.max = 15,
                          
                          B = 500)
#output:13
#Get number of clusters from gap statistic, maxSE is objective, gives a solid reason

#signatures
#Taking out the Gap statistic and the standard error and using the max gap value to identify the optimal number of clusters
krows <- maxSE(myclusGaprows$Tab[, "gap"], myclusGaprows$Tab[, "SE.sim"], method="globalSEmax")

#number of clusters for your samples 
kcols <- maxSE(myclusGapcols$Tab[, "gap"], myclusGapcols$Tab[, "SE.sim"], method="globalSEmax")


#creating a vector of colours to use to identify samples and signatures
my_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
               "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
               "#1A1A1A", "#FF4A46", "#008B8B", "#A0522D", "#800080",
               "#FFD700", "#D2691E", "#228B22", "#6B8E23", "#2F4F4F",
               "#F08080", "#00FA9A", "#FF1493")

#setting the krows and kcols values based on group discussion 
krows <- 3
kcols <-3

#cutting the tree at the desired number of clusters
cutclhr <- cutree(Hr_row, k = krows)

#assigning colours to samples to the cluster
cutclcolhr <-my_colors[1:krows]
cutclcolhr <- cutclcolhr[as.vector(cutclhr)]

#repeated for the samples
cutclhc <- cutree(HR_col, k = kcols)

cutclcolhc <-my_colors[1:kcols]
cutclcolhc <- cutclcolhc[as.vector(cutclhc)]


#disease 
diseasecolors = my_colors 

#creating a vector that indexes the diseases in the dataframe so each unique disease has a specific colour
diseasecolors_set <- my_colors[as.factor(Clara_T3$Disease)] 


#plots heatmap
plot <- heatmap3(as.matrix(t(TCGAClaraTsigsheatmap)), scale = "none", labCol = "",
                   Rowv = dd_row,RowSideColors = cutclcolhr, Colv = rev(dd_col),ColSideColors = cbind(diseasecolors_set,cutclcolhc),RowSideLabs = 'Signature Clusters', ColSideLabs = cbind('Disease' , 'Sample Clusters'),
                   col = greenred(500), balanceColor = TRUE, cexCol = 1.5, cexRow = 1.0, margins = c(7.5, 20))


#adds legend
legend('right', legend=unique(as.factor(Clara_T3$Disease)),
       
       fill=unique(diseasecolors_set), border=FALSE, bty="n", cex=1)

#################################################################################################################################
#################################################################################################################################
#enrichment of diseases in new clusters


table(cutclcolhc)
#1F77B4:blue is cluster 1
#2CA02C:green is cluster 2
#FF7F0E:orange is cluster 3

#Creating a vector that has the disease and corresponding cluster membership
Disease_df <- cbind(Clara_T3$Disease, cutclcolhc)

Cluster_1 <- subset(Disease_df, cutclcolhc == '#1F77B4')
table(Cluster_1)
dim(Cluster_1)


Cluster_2 <- subset(Disease_df, cutclcolhc == '#2CA02C')
table(Cluster_2)
dim(Cluster_2)

Cluster_3 <- subset(Disease_df, cutclcolhc == '#FF7F0E')
table(Cluster_3)
dim(Cluster_3)

#################################################################################################################################
#################################################################################################################################
#enrichment analysis for disease within each cluster 

#creating % of how the diseases make up the cluster
counts <- table(Cluster_1[,1])
percentages_culster_1 <- prop.table(counts,) * 100


#printing out the percentage value
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

#creating a function to calculate the spread of a disease across the clusters as %
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

#Sample how the function would be used for a specific disease
count_percentage('gbm')

#using the function across each unique disease in the DF
result <- lapply(unique(Clara_T3$Disease), count_percentage)

#create df for clusters
counts <- as.data.frame(counts)
percentages_culster_1<- as.data.frame(percentages_culster_1)
df_cluster_1 <- data.frame(Disease = counts$Var1,Count = counts$Freq, Percentage = percentages_culster_1$Freq)

counts_2 <- as.data.frame(counts_2)
percentages_culster_2<- as.data.frame(percentages_culster_2)
df_cluster_2 <- data.frame(Disease = counts_2$Var1, Count = counts_2$Freq, Percentage = percentages_culster_2$Freq)

counts_3 <- as.data.frame(counts_3)
percentages_culster_3<- as.data.frame(percentages_culster_3)
df_cluster_3 <- data.frame(Disease = counts_3$Var1, Count = counts_3$Freq, Percentage = percentages_culster_3$Freq)


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
library(ggplot2)

pie_chart_cluster_1 <- ggplot(df_cluster_1, aes(x = "", y = Percentage, fill = Disease)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values =  df_cluster_1$color) + 
  labs(title = "Disease Distribution in Cluster 1")

print(pie_chart_cluster_1)
ggsave(filename = "Disease Dis 1.jpg", plot = pie_chart_cluster_1, dpi = 600)

pie_chart_cluster_2 <- ggplot(df_cluster_2, aes(x = "", y = Percentage, fill = Disease)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values =  df_cluster_2$color) + 
  labs(title = "Disease Distribution in Cluster 2")

print(pie_chart_cluster_2)
ggsave(filename = "Disease Dis 2.jpg", plot = pie_chart_cluster_2, dpi = 600)

pie_chart_cluster_3 <- ggplot(df_cluster_3, aes(x = "", y = Percentage, fill = Disease)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values =  df_cluster_3$color) + 
  labs(title = "Disease Distribution in Cluster 3")

print(pie_chart_cluster_3)
ggsave(filename = "Disease Dis 3.jpg", plot = pie_chart_cluster_3, dpi = 600)

#################################################################################################################################
#################################################################################################################################

