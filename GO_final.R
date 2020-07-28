BiocManager::install(c("AnnotationHub", "clusterProfiler", "org.Pf.plasmo.db", "AnnotationDbi"))
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, "Malaria")
library(clusterProfiler)
Malaria <- hub[["AH77330"]]
Malaria
library(org.Pf.plasmo.db)
library(AnnotationDbi)

#Data Input
setwd("C:\\Users\\35826\\Desktop\\Dr. Mindy Shi\\Malaria dsRNA Transcriptome Project")
res <- read.csv("ALL_DE.csv", head = T) 
rnk_EL <- res$EarlyvsLate_LogFC
names(rnk_EL) <- res$Gene
rnk_EM <- res$EarlyvsMid_LogFC
names(rnk_EM) <- res$Gene
rnk_ML <- res$MidvsLate_LogFC
names(rnk_ML) <- res$Gene

### GO analysis
#Condition comparison Early vs Late
rnk_EL <- rnk_EL[!is.na(names(rnk_EL))]
head(rnk_EL)
View(rnk_EL)
#Rank in descending order by fold change
rnk_EL <- rnk_EL[order(rnk_EL, decreasing = TRUE)]
rnk_EL

gse <- gseGO(geneList=rnk_EL, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             pvalueCutoff = 1,
             OrgDb = org.Pf.plasmo.db)
View(gse@result)
write.csv (x=gse@result, file="GO_EL_R.csv")

#Condition comparison Early vs Mid
rnk_EM <- rnk_EM[!is.na(names(rnk_EM))]
head(rnk_EM)
View(rnk_EM)
rnk_EM <- rnk_EM[order(rnk_EM, decreasing = TRUE)]
rnk_EM

gse1 <- gseGO(geneList=rnk_EM, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             pvalueCutoff = 1,
             OrgDb = org.Pf.plasmo.db)
View(gse1@result)
write.csv (x=gse1@result, file="GO_EM_R.csv")

#Condition comparison Mid vs Late
rnk_ML <- rnk_ML[!is.na(names(rnk_ML))]
head(rnk_ML)
View(rnk_ML)
rnk_ML <- rnk_ML[order(rnk_ML, decreasing = TRUE)]
View(rnk_ML)
gse2 <- gseGO(geneList=rnk_ML, 
              ont ="BP", 
              keyType = "SYMBOL", 
              nPerm = 10000, 
              pvalueCutoff = 1,
              OrgDb = org.Pf.plasmo.db)
View(gse2@result)
write.csv (x=gse2@result, file="GO_ML_R.csv")

### KEGG
search_kegg_organism('Plasmodium falciparu')

#Condition comparison Early vs Late
kegg <- gseKEGG(geneList = rnk_EL, organism = "pfa", nPerm = 10000, pvalueCutoff = 1)
View(kegg@result)
write.csv (x=kegg@result, file="KEGG_EL.csv")

#Condition comparison Early vs Mild
kegg1 <- gseKEGG(geneList = rnk_EM, organism = "pfa", nPerm = 10000, pvalueCutoff = 1)
View(kegg1@result)
write.csv (x=kegg1@result, file="KEGG_EM.csv")

#Condition comparison Mild vs Late
kegg2 <- gseKEGG(geneList = rnk_ML, organism = "pfa", nPerm = 10000, pvalueCutoff = 1)
View(kegg2@result)
write.csv (x=kegg2@result, file="KEGG_ML.csv")
