Almac <- c('ANGPTL2', 'CDH11', "COL10A1", "COL5A1", "FAP", "FN1", "GFPT2", "GJB2", "INHBA", "MMP14", "PLAU", "RAB31", "THBS1", "THBS2", "VCAN")

Mak <- c('ADAM12',
         'ADAMTS12',
         'ADAMTS2',
         'AEBP1',
         'ANGPTL2',
         'ANTXR1',
         'AXL',
         'BNC2',
         'CALD1',
         'CDH2',
         'CMTM3',
         'CNRIP1',
         'COL10A1',
         'COL1A1',
         'COL1A2',
         'COL3A1',
         'COL5A1',
         'COL5A2',
         'COL6A1',
         'COL6A2',
         'COL6A3',
         'COL8A1',
         'DACT1',
         'EMP3',
         'FAP',
         'FBN1',
         'FN1',
         'FSTL1',
         'GPC6',
         'GYPC',
         'HTRA1',
         'INHBA',
         'ITGA11',
         'LOXL2',
         'LRRC15',
         'MMP2',
         'MSRB3',
         'NAP1L3',
         'NID2',
         'OLFML2B',
         'PCOLCE',
         'PDGFRB',
         'PMP22',
         'POSTN',
         'SPARC',
         'SPOCK1',
         'SULF1',
         'SYT11',
         'THBS2',
         'VCAN',
         'VIM',
         'ZEB2',
         'AP1G1',
         'ATP8B1',
         'CDH1',
         'CDS1',
         'CGN',
         'CLDN4',
         'CNOT1',
         'CTNND1',
         'DYNC1LI2',
         'ERBB3',
         'ESRP1',
         'ESRP2',
         'F11R',
         'GALNT3',
         'GPR56',
         'GRHL2',
         'HOOK1',
         'IRF6',
         'MAP7',
         'MARVELD2',
         'MARVELD3',
         'MYO5B',
         'OCLN',
         'PRSS8',
         'SPINT1')


Rokavec <- c( 'AFTPH',
              'CDH1',
              'CDS1',
              'EPCAM',
              'ESRP1',
              'ESRP2',
              'FA2H',
              'HOOK1',
              'LGALS1',
              'LNX1',
              'MAP7',
              'MAPK13',
              'MARVELD2',
              'MARVELD3',
              'MYO5B',
              'MYO5C',
              'P4HA3',
              'RBM47',
              'UROD',
              'VIM')

Dry_functional <-scan("gene_list.txt", what = "character", sep = "\n")



cluster_1 <- c(schell, Rokavec, apical_surface, gentles, apical_junction, dry_resistance, Lee, EMT_transition, becht, Almac)

print(cluster_1)

cluster_2 <- c(beta_catenin, notch, hedgehog, Tan, thompson, Mak, byers, wagle, dry_functional)

total <- c(cluster_1, cluster_2)



paste(total)

dput(total)


install.packages("clusterProfiler")
library(clusterProfiler)

?enrichGO

aa <- enrichGO(gene= beta_catenin, universe = names(total), 'org.Hs.eg.db')



