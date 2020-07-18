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
res <- read.csv("DE.csv", head = T) 
rnk_RS <- res$Ring.vsSchizont_LogFC
names(rnk_RS) <- res$Gene
rnk_RT <- res$Ring.VsTroph_LogFC
names(rnk_RT) <- res$Gene
rnk_TS <- res$Troph.vsSchizont_LogFC
names(rnk_TS) <- res$Gene

### GO analysis
#Condition comparison Ring vs Schizont
rnk_RS <- rnk_RS[!is.na(names(rnk_RS))]
head(rnk_RS)
View(rnk_RS)
#Rank in descending order by fold change
rnk_RS <- rnk_RS[order(rnk_RS, decreasing = TRUE)]
rnk_RS

gse <- gseGO(geneList=rnk_RS, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             pvalueCutoff = 1,
             OrgDb = org.Pf.plasmo.db)
View(gse@result)
write.csv (x=gse@result, file="GO_RS_R.csv")

#Condition comparison Ring vs Troph
rnk_RT <- rnk_RT[!is.na(names(rnk_RT))]
head(rnk_RT)
View(rnk_RT)
rnk_RT <- rnk_RT[order(rnk_RT, decreasing = TRUE)]
rnk_RT

gse1 <- gseGO(geneList=rnk_RT, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             pvalueCutoff = 1,
             OrgDb = org.Pf.plasmo.db)
View(gse1@result)
write.csv (x=gse1@result, file="GO_RT_R.csv")

#Condition comparison Troph vs Schizont
rnk_TS <- rnk_TS[!is.na(names(rnk_TS))]
head(rnk_TS)
View(rnk_TS)
rnk_TS <- rnk_TS[order(rnk_TS, decreasing = TRUE)]
View(rnk_TS)
gse2 <- gseGO(geneList=rnk_TS, 
              ont ="BP", 
              keyType = "SYMBOL", 
              nPerm = 10000, 
              pvalueCutoff = 1,
              OrgDb = org.Pf.plasmo.db)
View(gse2@result)
write.csv (x=gse2@result, file="GO_TS_R.csv")

### KEGG
search_kegg_organism('Plasmodium falciparu')

#Condition comparison Ring vs Schizont
kegg <- gseKEGG(geneList = rnk_RS, organism = "pfa", nPerm = 10000, pvalueCutoff = 1)
View(kegg@result)
write.csv (x=kegg@result, file="KEGG_RS.csv")

#Condition comparison Ring vs Troph
kegg1 <- gseKEGG(geneList = rnk_RT, organism = "pfa", nPerm = 10000, pvalueCutoff = 1)
View(kegg1@result)
write.csv (x=kegg1@result, file="KEGG_RT.csv")

#Condition comparison Troph vs Schizont
kegg2 <- gseKEGG(geneList = rnk_TS, organism = "pfa", nPerm = 10000, pvalueCutoff = 1)
View(kegg2@result)
write.csv (x=kegg2@result, file="KEGG_TS.csv")
