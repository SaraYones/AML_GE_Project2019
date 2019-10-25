if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("topGO")


BiocManager::install("clusterProfiler", version = "3.8")
de <- classifierMCFSgenes
background<-colnames(TARGET_GE_Classifier)[-dim(TARGET_GE_Classifier)[2]]

gene.df <- bitr(de, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

x<-org.Hs.egENSEMBL
yy <- enrichGO(de,'org.Hs.eg.db',keyType= 'ENSEMBL', ont="BP",universe=background, pvalueCutoff=0.1)
filteredBoruta=getAnnotatedGenes(genes,BorutaFeatures)



library(org.Hs.eg.db)
data(geneList)
gene <- names(geneList)[abs(geneList) >= 1]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)



