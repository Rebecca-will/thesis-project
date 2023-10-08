acc <- data_clinical_sample_acc[['#Patient Identifier', 'Tissue Sort Site']]

dim(data_clinical_sample_acc)
colnames(data_clinical_sample_acc)


location <- rbind(data_clinical_sample_acc, data_clinical_sample_blca, data_clinical_sample_brca, data_clinical_sample_cesc, data_clinical_sample_coadread, data_clinical_sample_gbm, data_clinical_sample_hnsc, data_clinical_sample_kich, data_clinical_sample_kirc, data_clinical_sample_kirp, data_clinical_sample_lgg, data_clinical_sample_lihc, data_clinical_sample_luad, data_clinical_sample_lusc, data_clinical_sample_ov, data_clinical_sample_paad, data_clinical_sample_prad, data_clinical_sample_skcm, data_clinical_sample_stad, data_clinical_sample_tgct, data_clinical_sample_thca, data_clinical_sample_ucec)

library(dplyr)

location<-bind_rows(location, data_clinical_sample_esca)

location <- slice(location, -(6978))

dim(location)


filtered_df <- location[ c(1,3, 18)]
#esca 

colnames(data_clinical_sample_esca)
colnames(data_clinical_sample_acc)


filtered_df <- rename(filtered_df, 'PATIENT_ID' = "#Patient Identifier" )


loc_df <- merge(Clara_T3, filtered_df, by = 'PATIENT_ID')


#################################################################################################################################
#################################################################################################################################


n <- 155  
palette <- colorRampPalette(c("#FF0000", "#0000FF"))  
loc_colors <- palette(n)

unique(loc_colors)


library(colorspace)

# Generate 155 unique hexadecimal colors
n <- 155  # Number of colors in the vector
hex_colors <- rainbow_hcl(n, start = 0, end = 360, c = 100, l = 50, alpha = 1)



loc_colours = hex_colors
loc_colours_set <- loc_colours[as.factor(loc_df$`Tissue Source Site`)]

unique(loc_colours_set)


plot_2 <- heatmap3(as.matrix(t(TCGAClaraTsigsheatmap)), scale = "none", labCol = "",
                 Rowv = dd_row,RowSideColors = cutclcolhr, Colv = rev(dd_col),ColSideColors = cbind(loc_colours_set,cutclcolhc),RowSideLabs = 'Signature Clusters', ColSideLabs = cbind('Tissue Source Site' , 'Sample Clusters'),
                 col = greenred(500), balanceColor = TRUE, cexCol = 1.5, cexRow = 1.0, margins = c(7.5, 20))

legend(x=0.02, y= 0.85, legend=unique(as.factor(loc_df$`Tissue Source Site`)),
       
       fill=unique(loc_colours_set), border=FALSE, bty="n", cex=0.5)


print(cutclcolhc)

dev.off()



#use signatures and see if any in common 


