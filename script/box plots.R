library(ggplot2)
library(dplyr)
mes_sample <- subset(No_mes, select = c(2,3))

emt_sig <- merge(Gene_sig_nh, mes_sample, by = 'Sample')
dim(emt_sig)

emt_sig_1<- subset(emt_sig, select = -c(1134,1))
emt_sig_1 <- emt_sig_1[!(emt_sig_1$cluster == '1'),]

emt_sig_1$cluster = ifelse(emt_sig_1$cluster==2,"2-classical",ifelse(emt_sig_1$cluster==3,"3-no emt",NA))

plot_emt <- list()
# Iterate through each disease cluster
for (i in top_5) {

  
 violin_plot<- emt_sig_1 %>%
    filter(cluster %in% c( '2-classical', '3-no emt')) %>%
    group_by(cluster, !!sym(as.character(i))) %>%
    ggplot( aes(x = cluster, y = !!sym(as.character(i)), fill = cluster)) +
    geom_violin()+
    scale_fill_manual(values = c("#2CA02C", '#FF7F0E')) +
    labs(x = "cluster", y = "Expression", fill = "Cluster") +
    ggtitle("Violon plot of ", i) +
    theme_bw()

  
  # Add the survival plot to the list
  plot_emt[[i]] <- violin_plot
  
}

print(plot_emt)

emt_sig_1 %>%
  filter(cluster %in% c( '2-classical', '3-no emt')) %>%
  group_by(cluster) %>%
  ggplot( aes(x = cluster, y = LOXL2, fill = cluster)) +
  geom_violin()+
  scale_fill_manual(values = c("#2CA02C", '#FF7F0E')) +
  labs(x = "cluster", y = "Expression", fill = "Cluster") +
  ggtitle("Violon plot of LOXL2") +
  theme_bw()


emt_sig_1$cluster

RFE_DF$cluster_name = ifelse(RFE_DF$cluster=="2",'2-classical',ifelse(RFE_DF$cluster=="3",'3-no emt',NA))


RFE_DF %>%
  filter(cluster_name %in% c( '2-classical', '3-no emt')) %>%
  group_by(cluster_name, SULF2) %>%
  ggplot( aes(x = cluster, y = SULF2, fill = cluster_name)) +
  geom_violin()+
  scale_fill_manual(values = c("#2CA02C", '#FF7F0E')) +
  labs(x = "cluster", y = "Expression", fill = "Cluster") +
  ggtitle("Violon plot of SULF2") +
  theme_bw()

emt_sig$cluster<- as.factor(emt_sig$cluster)

emt_sig %>%
  filter(cluster %in% c( '2', '3')) %>%
  group_by(cluster) %>%
  ggplot( aes(x = cluster, y = ACTA2, fill = cluster)) +
  geom_violin()+
  scale_fill_manual(values = c("#2CA02C", '#FF7F0E')) +
  labs(x = "cluster", y = "Expression", fill = "Cluster") +
  ggtitle("Violon plot of ACTA2") +
  theme_bw()
