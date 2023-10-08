# the et al are papers from where? 
#we have the columns but are they actual measures of what?


# ACC	Adrenocortical carcinoma
# BLCA	Bladder Urothelial Carcinoma
# LGG	Brain Lower Grade Glioma
# BRCA	Breast invasive carcinoma
# CESC	Cervical squamous cell carcinoma and endocervical adenocarcinoma
# COAD	Colon adenocarcinoma
# CNTL	Controls
# ESCA	Esophageal carcinoma
# GBM	Glioblastoma multiforme
# HNSC	Head and Neck squamous cell carcinoma
# KICH	Kidney Chromophobe
# KIRC	Kidney renal clear cell carcinoma
# KIRP	Kidney renal papillary cell carcinoma
# LIHC	Liver hepatocellular carcinoma
# LUAD	Lung adenocarcinoma
# LUSC	Lung squamous cell carcinoma
# OV	Ovarian serous cystadenocarcinoma
# PAAD	Pancreatic adenocarcinoma
# PRAD	Prostate adenocarcinoma
# SKCM	Skin Cutaneous Melanoma
# STAD	Stomach adenocarcinoma
# TGCT	Testicular Germ Cell Tumors
# THCA	Thyroid carcinoma
# UCEC	Uterine Corpus Endometrial Carcinoma

#what exaclty is the claraT tool- mRNA
#there are et al. so i assume they're papers 
#there are adenomcarcinoas not included like rectum, sarcoma etc
#so do we separate the data into their respective groupings or analyse as a whole? 
# does the claraT tool take out the qualtiy check
#theres no controls so comparing enrichment to what? how are we defining a baseline
table(Clara_T$Disease)
dim(Clara_T)
library(cluster)


#extract ID
ID<- Clara_T[, 1:2]


#scale the data 
scaled_data<- apply(Clara_T[,3:30], 2, scale)
scaled_data

#combine 
final_scaled <- cbind(ID,scaled_data)

#distances

Row_dist <- dist(t(scaled_data), method='euclidean')

Col_dist <- dist(scaled_data, method='euclidean')

#Clustering 
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
myclusGapcols <- clusGap(scaled_data,
                         
                         FUN = mycluster,
                         
                         K.max = 15,
                         
                         B = 500)
#output:15
#cluster in the same way as before and tell how similar the clusters are to eachother and what is the optimal number of clusters


#signature clusters = 2

myclusGaprows <- clusGap(t(scaled_data),
                         
                         FUN = mycluster,
                         
                         K.max = 15,
                         
                         B = 500)
#output:13
##put into console
myclusGapcols
plot(myclusGapcols)
help(maxSE) #into console, swap in with rows, rerun the krows/kcols
#change the options in line 126 to see what optimum clusters come out 



#Get number of clusters from gap statistic, maxSE is objective, gives a solid reason

krows <- maxSE(myclusGaprows$Tab[, "gap"], myclusGaprows$Tab[, "SE.sim"], method="globalSEmax")
krows <- 13
#signatures

?maxSE
kcols <- maxSE(myclusGapcols$Tab[, "gap"], myclusGapcols$Tab[, "SE.sim"], method="globalSEmax")
kcols <- 15
#number of clusters for your samples 

install.packages('RColorBrewer')
library(RColorBrewer)
display.brewer.all(n=23)

my_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
               "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
               "#1A1A1A", "#FF4A46", "#008B8B", "#A0522D", "#800080",
               "#FFD700", "#D2691E", "#228B22", "#6B8E23", "#2F4F4F",
               "#F08080", "#00FA9A", "#FF1493")

colors1=brewer.pal(23, 'Set3')

colors=rainbow(15)#keep this in mind to the k.max, colour to each cluster


cutclhr <- cutree(Hr_row, k = krows)#cut the tree at many clusters as gap has decided 

cutclcolhr <-colors[1:krows]#assigning colours to samples to the cluster

cutclcolhr <- cutclcolhr[as.vector(cutclhr)]


#samples
cutclhc <- cutree(HR_col, k = kcols)

cutclcolhc <-colors[1:kcols]

cutclcolhc <- cutclcolhc[as.vector(cutclhc)]


#disease 
diseasecolors = my_colours #change to different pallette later 

diseasecolors_set <- my_colors[as.factor(Clara_T$Disease)] #changed to as.factor since its characters

plot_2 <- heatmap3(as.matrix(t(TCGAClaraTsigsheatmap)), scale = "none", labCol = "",
                   Rowv = dd_row,RowSideColors = cutclcolhr, Colv = rev(dd_col),ColSideColors = cbind(diseasecolors_set,cutclcolhc),RowSideLabs = 'Signature Clusters', ColSideLabs = cbind('Disease' , 'Sample Clusters'),
                   col = greenred(500), balanceColor = TRUE, cexCol = 1.5, cexRow = 1.0, margins = c(7.5, 20))


#Have adjusted and looks okay on the zoomed in display
legend(x = 0.02, y = 1, legend=unique(as.factor(Clara_T$Disease)),
       
       fill=unique(my_colors), border=FALSE, bty="n", cex=0.5)

dev.off()



