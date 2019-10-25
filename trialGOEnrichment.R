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
gene.df <- bitr(gene, fromType ="ENTREZID",
                toType = c("ENSEMBL"),
                OrgDb = org.Hs.eg.db)


GOenrichment=function(features,background,flag)
{
  
  if (flag=="ENSEMBL")
  {
    
  ego <- enrichGO(gene          = features,
                  universe      = background,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.01)
  }
  return(ego)
  
  
  
}



de=Linda_MCFSFeatures[1:700]
background=colnames(Linda_GE_Classifier)[1:length(Linda_GE_Classifier)-1]
ego <- enrichGO(gene          = de,
                universe      = background,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)


de=Integrated_Matched_MCFSFeatures[1:20]
background=colnames(Integrated_Matched_Classifier)[1:length(Integrated_Matched_Classifier)-1]
ego <- enrichGO(gene          = de,
                universe      = background,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

de=MCFSFeatures[1:500]
background=colnames(TARGET_GE_Classifier)[1:length(TARGET_GE_Classifier)-1]
ego <- enrichGO(gene          = de,
                universe      = background,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)



data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

geneList.df <- bitr(names(geneList), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)


head(gene.df)
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 universe      = geneList.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

yy <- enrichGO(gene         = de,
              universe      = background,
               OrgDb         = org.Hs.eg.db,
               keyType       = 'ENSEMBL',
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)

ggo <- groupGO(gene     = de,
              # universe      = background,
               OrgDb    = org.Hs.eg.db,
               keyType       = 'ENSEMBL',
               ont      = "BP",
               level    = 3)

gene.df <- bitr(de, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

#Try with EntrezID (was done by bioMart the conversion)

genes.with.id=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=de, mart= ensembl) 
trial=genes.with.id[which(!is.na(genes.with.id[which(genes.with.id$ensembl_gene_id %in% de),])),]
trial=trial[which(!is.na(trial$entrezgene)),]

background=colnames(Linda_GE_Classifier)[-dim(Linda_GE_Classifier)[2]]

genes.with.id.bg=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=background, mart= ensembl) 
trial.bg=genes.with.id.bg[which(!is.na(genes.with.id.bg[which(genes.with.id.bg$ensembl_gene_id %in% background),])),]
trial.bg=trial.bg[which(!is.na(trial.bg$entrezgene)),]


genes.with.id=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=de, mart= ensembl) 
trial=genes.with.id[which(!is.na(genes.with.id[which(genes.with.id$ensembl_gene_id %in% de),])),]
trial=trial[which(!is.na(trial$entrezgene)),]

genes.with.id.bg=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=background, mart= ensembl) 
trial.bg=genes.with.id.bg[which(!is.na(genes.with.id.bg[which(genes.with.id.bg$ensembl_gene_id %in% background),])),]
trial.bg=trial.bg[which(!is.na(trial.bg$entrezgene)),]


yy <- enrichGO(gene         = as.character(trial$entrezgene),
               universe      =as.character(trial.bg$entrezgene),
               OrgDb         = org.Hs.eg.db,
             #  keyType       = 'ENSEMBL',
               ont           = "CC",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)

dotplot(yy, showCategory=30)
#Without using BG
x <- enrichDO(gene          = as.character(trial$entrezgene),
              ont           = "DO",
            #  keyType       = 'ENSEMBL',
              pvalueCutoff  = 0.01,
              pAdjustMethod = "BH",
              universe      = as.character(trial.bg$entrezgene),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.01,
              readable      = FALSE)

#Try with gene symbols alone without any transformation


#Protein Coding Genes only

geneProtein=append(filteredMCFS[[1]]$geneID,filteredMCFS[[2]])
geneProteinbg=append(filter[[1]]$geneID,filter[[2]])

ego <- enrichGO(gene          = geneProtein,
             #   universe      = geneProteinbg,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)


#I didnt get any significant BP when using protein coding genes


#--------------------------Trying GOSim and similarties between Go enrichment

result=rosetta(Linda_GE_Classifier[,append(Linda_MCFSFeatures[1:50],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",reducer="Genetic",ruleFiltration=TRUE, discreteParam=3)
result=recalculateRules(Linda_GE_Classifier[,append(Linda_MCFSFeatures[1:50],"decision")],result$main)
features=getFeatures(Linda_GE_Classifier[,append(Linda_MCFSFeatures[1:50],"decision")],result,1.0)
features=unique(append(levels(features[[1]]),levels(features[[2]])))



de=Linda_classifierMCFSgenes
genes.with.id.bg=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=de, mart= ensembl) 
trial.bg=genes.with.id.bg[which(!is.na(genes.with.id.bg[which(genes.with.id.bg$ensembl_gene_id %in% background),])),]
trial.bg=trial.bg[which(!is.na(trial.bg$entrezgene)),]


#--------Checking the similarity of genes in the rules using GOSemSim package (According to their GO terms)
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="BP")
hsGO2 <- godata('org.Hs.eg.db', keytype = "ENSEMBL", ont="MF", computeIC=FALSE)
hsGO3 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)



features=getFeatures(Linda_GE_Classifier_MCFS,Linda_resultRosettaMCFS$main,1.0)

features=unique(append(levels(features[[1]]),levels(features[[2]])))

trialAnnotated=mgeneSim(features, semData=hsGO3, measure="Wang", combine="BMA", verbose=FALSE)

trialAnnotated[which(trialAnnotated<0.5,arr.ind=TRUE)]=0

unique(sort(trial))


result=rosetta(Linda_GE_Classifier[,append(Linda_MCFSFeatures[1:10],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3)

featuresENSMBL=getFeatures(Linda_GE_Classifier[,append(Linda_MCFSFeatures[1:10],"decision")],result$main,1.0)

featuresENSMBL=unique(append(levels(features[[1]]),levels(features[[2]])))

trialUnAnnotated=mgeneSim(featuresENSMBL  , semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)



my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green


heatmap.2(as.matrix(trialAnnotated),
       #   cellnote = as.matrix(trialAnnotated),  # same data set for cell labels
          main = "similarity", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks)   # enable color transition at specified limits
         # dendrogram="row",     # only draw a row dendrogram
         # Colv="NA")            # turn off column clustering

#--------------Feature selection based on Genes and priortizing Go terms

BiocManager::install("GOexpress", version = "3.8")


set.seed(4543)

listMarts(host='feb2014.archive.ensembl.org')
ensembl75 = useMart( host='feb2014.archive.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl')
allgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'external_gene_id', 'description'),mart=ensembl75)
colnames(allgenes.Ensembl)[1] = 'gene_id'
allGO.Ensembl = getBM(attributes=c('go_id', 'name_1006', 'namespace_1003'),mart=ensembl75)
GOgenes.Ensembl = getBM(attributes=c('ensembl_gene_id', 'go_id'),mart=ensembl75)
colnames(GOgenes.Ensembl)[1] = 'gene_id'
GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$go_id != '',]
GOgenes.Ensembl = GOgenes.Ensembl[GOgenes.Ensembl$gene_id != '',]

#-------------The real analysis-------------------------------

pd <- new("AnnotatedDataFrame", data =as.data.frame(Linda_GE_Classifier$decision))
rownames(pd)<-rownames(Linda_GE_Classifier)
exprs=t(as.matrix(Linda_GE_Classifier[,1:dim(Linda_GE_Classifier)[2]-1]))
all(rownames(pData)==colnames(exprs))
exampleSet <- ExpressionSet(assayData=exprs,phenoData=pd)


set.seed(4543) # set random seed for reproducibility
results <- GO_analyse(eSet = exampleSet, f = "Linda_GE_Classifier$decision", GO_genes=GOgenes.Ensembl, all_GO=allGO.Ensembl, all_genes=allgenes.Ensembl)
names(results)

head(results$GO[, c(1:5, 7)], n=5) #
head(results$genes[, c(1:3)], n=10)


#Binding Genes to GOterms and their annotation
head(results$mapping)
annotationGenes=lapply(results$mapping$go_id, function(x) (which(allGO.Ensembl$go_id %in% x)))
annotationGenes=unlist(annotationGenes)
resultsGenesDF=cbind(results$mapping,allGO.Ensembl[annotationGenes,"name_1006"])

Linda_MCFSFeatures_GOannotations1=lapply(Linda_MCFSFeatures[1:10], function(x) (results$mapping$go_id[which(results$mapping$gene_id %in% x)]))
Linda_MCFSFeatures_GOannotations2=lapply(Linda_MCFSFeatures_GOannotations1, function(x) (allGO.Ensembl$name_1006[which(allGO.Ensembl$go_id %in% x)]))
Linda_MCFSFeatures_GOannotations=cbind(Linda_MCFSFeatures[1:10],Linda_annotation,Linda_MCFSFeatures_GOannotations1,Linda_MCFSFeatures_GOannotations2) 
Linda_MCFSFeatures_GOannotations[,3]=lapply(Linda_MCFSFeatures_GOannotations[,3],function(x) x=paste(unlist(x), collapse=';'))
Linda_MCFSFeatures_GOannotations[,4]=lapply(Linda_MCFSFeatures_GOannotations[,4],function(x) x=paste(unlist(x), collapse=';'))



results.pVal = pValue_GO(results, N=100)

BP.5 <- subset_scores(result = results.pVal,namespace = "biological_process",total = 5, p.val=0.05)

MF.10 <- subset_scores(result = results.pVal, namespace = "molecular_function",total = 10,p.val=0.05)


CC.15 <- subset_scores(result =results.pVal,namespace = "cellular_component",total = 15,p.val=0.05)

heatmap_GO(go_id = "GO:0002328", result = BP.5, eSet=exampleSet, cexRow=0.4, cexCol=1, cex.main=1, main.Lsplit=30)
