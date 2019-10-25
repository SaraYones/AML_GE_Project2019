library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)


#Protein coding genes and those that are unannotated (To run feature selection only on protein coding genes)
filter=getAnnotatedGenes(genes,colnames(TARGET_GE_Classifier))
protein_coding=filter[[1]][which(filter[[1]]$gene_type =="protein_coding"),geneID]
non_pseudo_genes=filter[[1]][which(!(grepl("*\\_pseudogene",filter[[1]]$gene_type))),]$geneID
#Try Feature selection with all RNAs or genes except psudo genes and IGV


#Have to run TARGETAMLPreprocess and TARGETAMLfunctions first 
TARGET_GE_Classifier_ProteinCoding=cbind.data.frame(TARGET_GE_Classifier[,append(protein_coding,filter[[2]])],decision)
TARGET_GE_Classifier_nonpseudo=cbind.data.frame(TARGET_GE_Classifier[,append(non_pseudo_genes,filter[[2]])],decision)
#All Genes
TARGET_GE_Classifier=cbind.data.frame(TARGET_GE_Classifier,decision)
#MCFS and Boruta are run on ULAM
write.table(TARGET_GE_Classifier,"MCFS_AML_GE_Classifier",row.names=TRUE)
write.table(TARGET_GE_Classifier_nonpseudo,"MCFS_AML_GE_Classifier",row.names=TRUE)
#After filtering out according to CPM normalization
#matched 3  folder is output of MCFS for all genes after filtering according edgeR
#matched 4 folder is output of MCFS for protein coding and unannotated genes after filtering according edgeR
#matched 5 is output of MCFS for non pseudo genes after filtering according edgeR
MCFSFeatures=FilterFeatures("MCFSoutput/matched 4/matched__RI.csv",200)
MCFSFeatures=FilterFeatures("MCFSoutput/matched 3/matched__RI.csv",1000)
MCFSFeatures=FilterFeatures("MCFSoutput/matched5/matched5__RI.csv",200)

MCFSFeaturesnonpseudo=FilterFeatures("MCFSoutput/matched5/matched5__RI.csv",200)

#---Before applying the CPM filter using edge R method
boruta.train=readRDS("TARGET_AML_GE_BORUTA.rds")
#After filtering out according to CPM using edge R method in normalization
boruta.train=readRDS("TARGET_AML_GE_BORUTA1.rds")
BorutaFeatures=extractFeaturesBoruta(boruta.train)
#Only confirmed features
  
geneLS=list(MCFS=MCFSFeatures[1:30],Boruta=MCFSFeaturesnonpseudo[1:40])
intersection=plotVenn(geneLS,list(c("darkmagenta", "darkblue"),c(0.5,0.5),2,2, c("ProteinCoding","NonPseudo")))

intersection=plotVenn(geneLS,list(c("darkmagenta", "darkblue"),c(0.5,0.5),2,2, c("MCFS","Boruta")))
# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)



ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensmblid <- getBM(attributes=c('hgnc_symbol','chromosome_name','start_position','end_position'),
  filters=('ensembl_gene_id'), values =MCFSFeatures, mart = ensembl)

ensmblid <- getBM(attributes=c('hgnc_symbol','chromosome_name','start_position','end_position'),
                  filters=('ensembl_gene_id'), values = MCFSFeatures[which(!(MCFSFeatures %in% genes$geneID  ))], mart = ensembl)
filteredAll=getAnnotatedGenes(genes,MCFSFeatures)

genes <- readGeneDB("gencode.v20.annotation.gtf")
genes$geneID=gsub("(*)\\.([0-9]+)","\\1",genes$gene_id)
#20 Features when using all genes 
#30 Features when running protein coding genes on MCFS
#40 Features when running MCFS on non psuedo genes
filteredMCFS=getAnnotatedGenes(genes,MCFSFeatures[1:20])
filteredintersection=getAnnotatedGenes(genes,intersection$`MCFS:Boruta`)
filteredMCFS=getAnnotatedGenes(genes,MCFSFeatures[1:30])
filteredMCFS=getAnnotatedGenes(genes,MCFSFeatures[1:40])
#MCFS genes after removing all pseudo genes and comparing accuracies than getting the optimal number of features
filteredMCFSNoPseudo=getAnnotatedGenes(genes,classifierMCFSgenesNoPseudo[1:17])
write.table(filteredMCFS[[1]][,-c("attributes")],"filteredGenesMCFS",row.names=TRUE)
filteredBoruta=getAnnotatedGenes(genes,BorutaFeatures)
write.table(filteredBoruta[[1]][,-c("attributes")],"filteredGenesBoruta",row.names=TRUE)

plotNames=gsub("(TARGET-20-((.)+-(R|T)))","\\2",rownames(TARGET_GE_Classifier))

plotDistributionCategory(list(filteredMCFS[[1]]$gene_type,filteredBoruta[[1]]$gene_type),c("MCFS","Boruta"),"Distribution_Classes_FS",getwd(),"")

#Only MCFS after removing the pseudo genes from the first 20 genes
plotDistributionCategory(list(filteredBoruta[[1]]$gene_type),c("Boruta"),"Distribution_Classes_FS_MCFS_WithoutPseudoGenes","ResultsRules/AllGenes","")

#Plot distribution of expirement 2 with removing psuedo genes from 200 MCFS list then getting 17 as the optimal number of features
plotDistributionCategory(list(filteredMCFSNoPseudo[[1]]$gene_type),c("MCFS"),"Distribution_Classes_FS","ResultsRules/AllGenes","")


#Without PseudoGenes MCFS
classifierMCFSgenes=append(filteredMCFS[[1]][which(!(grepl("*\\_pseudogene",filteredMCFS[[1]]$gene_type))),]$geneID,filteredMCFS[[2]])
#Without PseudoGenes MCFS
classifierMCFSgenesNoPseudo=append(filteredAll[[1]][which(!(grepl("*\\_pseudogene",filteredAll[[1]]$gene_type))),]$geneID,filteredAll[[2]])

#All Genes MCFS
classifierMCFSgenes=append(filteredMCFS[[1]]$geneID,filteredMCFS[[2]])
#Ordering according to MCFS Features
classifierMCFSgenes=classifierMCFSgenes[!is.na(match(MCFSFeatures, classifierMCFSgenes))]

#All Genes Boruta
classifierBorutagenes=append(filteredBoruta[[1]]$geneID,filteredBoruta[[2]])
#Ordering according to MCFS Features
classifierBorutagenes=classifierBorutagenes[!is.na(match(BorutaFeatures, classifierBorutagenes))]




annotation=annotateInOrder(filteredMCFS,classifierMCFSgenes)

#classifierMCFSgenes=MCFSFeatures[1:30]
#Without PseudoGenes
classifierBorutagenes=append(filteredBoruta[[1]][which(!(grepl("*\\_pseudogene",filteredBoruta[[1]]$gene_type))),]$geneID,filteredBoruta[[2]])
#All Genes
classifierBorutagenes=append(filteredBoruta[[1]]$geneID,filteredBoruta[[2]])

my.plots <- vector(2, mode='list')

plotPCA(getwd(),TARGET_GE_Classifier[,classifierMCFSgenes],decision,"PCA_MCFS")
my.plots[[1]]=recordPlot()
plotPCA(getwd(),TARGET_GE_Classifier[,classifierBorutagenes],decision,"PCA_Boruta")
my.plots[[2]]=recordPlot()
#savePDF(my.plots,"PCA_AllGenes",getwd())
savePDF(my.plots,"PCA_WithoutPseudogenes",getwd())

TARGET_GE_Classifier_MCFS=TARGET_GE_Classifier[,append(classifierMCFSgenes,"decision")]
#colnames(TARGET_GE_Classifier_MCFS)=append(filteredMCFS[[1]]$gene_name[which(filteredMCFS[[1]]$geneID %in% classifierMCFSgenes)],append(filteredMCFS[[2]],"decision"))
colnames(TARGET_GE_Classifier_MCFS)=append(annotation,"decision") 

TARGET_GE_Classifier_Boruta=TARGET_GE_Classifier[,append(classifierBorutagenes,"decision")]
colnames(TARGET_GE_Classifier_Boruta)=append(filteredBoruta[[1]]$gene_name[which(filteredBoruta[[1]]$geneID %in% classifierBorutagenes)],append(filteredBoruta[[2]],"decision"))
#Compare accuracies on all MCFS features even with the pseudoGenes included
Accuracies=compareAccuracies(TARGET_GE_Classifier,200,MCFSFeatures)

Accuracies=compareAccuracies(TARGET_GE_Classifier,50,BorutaFeatures)
#Compare accuracies on all MCFS features but without pseudoGenes included
Accuracies=compareAccuracies(TARGET_GE_Classifier,length(classifierMCFSgenesNoPseudo),classifierMCFSgenesNoPseudo)
#Compare accuracies on all MCFS features when running MCFS on only protein coding genes 
Accuracies=compareAccuracies(TARGET_GE_Classifier,length(classifierMCFSgenesNoPseudo),classifierMCFSgenesNoPseudo)

Accuracies=compareAccuracies(TARGET_GE_Classifier,200,MCFSFeatures)

resultRosettaMCFS=rosetta(TARGET_GE_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3)

recalculatedResultRosettaMCFS=recalculateRules(TARGET_GE_Classifier_MCFS,resultRosettaMCFS$main)

#write.csv(recalculatedResultRosettaMCFS,paste("ResultsRules/AllGenes/RulesAllGenes-",Sys.Date(),".csv",sep=""))
#saveLineByLine(recalculatedResultRosettaMCFS,  paste("ResultsRules/AllGenes/NetworksAllGenes-",Sys.Date(),".txt",sep=""))

write.csv(recalculatedResultRosettaMCFS,paste("ResultsRules/ProteinCoding/RulesAllGenes-",Sys.Date(),".csv",sep=""))
saveLineByLine(recalculatedResultRosettaMCFSGenetic, paste("ResultsRules/AllGenes/NetworksAllGenes-",Sys.Date(),".txt",sep=""))

resultRosettaMCFSGenetic=rosetta(TARGET_GE_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3,reducer="Genetic",ruleFiltration=TRUE, ruleFiltrAccuracy=c(0,0.8),ruleFiltrSupport=c(1,3))
#resultRosettaMCFSGenetic=rosetta(TARGET_GE_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3,reducer="Genetic")
recalculatedResultRosettaMCFSGenetic=recalculateRules(TARGET_GE_Classifier_MCFS,resultRosettaMCFSGenetic$main)
write.csv(recalculatedResultRosettaMCFS,"ResultsRules/AllGenes/RulesGenetic.csv")
saveLineByLine(recalculatedResultRosettaMCFSGenetic, "ResultsRules/AllGenes/NetworksAllGenesGenetic.txt")

resultRosettaBoruta=rosetta(TARGET_GE_Classifier_Boruta,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3)
clusteredRulesMCFS=clusterRules(recalculatedResultRosettaMCFS,rownames(TARGET_GE_Classifier_MCFS))
colnames(clusteredRulesMCFS)<-plotNames


#my.plots=vector(1, mode='list');
#svg(paste("ResultsRules/AllGenes/HeatMapAllGenes-",Sys.Date(),".svg",sep=""))
#clusters=heatmap.F(t(clusteredRulesMCFS), colors=c('white','white','white','white','blue','blue'),distmethod='pearson')
#dev.off()
#write.csv(clusters,paste("ResultsRules/AllGenes/Clusters-",Sys.Date(),".csv",sep = ""))
#my.plots[[1]]=recordPlot()
#savePDF(my.plots,paste("HeatMapAllGenes-",Sys.Date(),"ResultsRules/AllGenes/"))

ordered_metadata_for_exploratory=writeOutput("ResultsRules/",Sys.Date(),"AllGenes",clusteredRulesMCFS)

#Merge output clusters with metadata exploratory
#paste "TARGET-20-  to the names of the clusters then match it with the name of the clusters, get the values and add it to a new coloumn in metadata_exploratory
metadata_for_exploratory$clusters=clusters[match(rownames(metadata_for_exploratory),paste(rep("TARGET-20-",length(clusters)),names(clusters),sep=""))]
#Order based on the order of the clusters 
i1=match(paste(rep("TARGET-20-",length(clusters)),names(clusters),sep=""),rownames(metadata_for_exploratory))
ordered_metadata_for_exploratory=metadata_for_exploratory[i1,]

plotClustersMetadataPlot(ordered_metadata_for_exploratory,unique(clusters),c("Protocol","FAB.Category","t.8.21.","inv.16.","MLL","FLT3.ITD.positive.","NPM.mutation","CEBPA.mutation","WT1.mutation","CR.status.at.end.of.course.1"),"ResultsRules/","2019-01-07","AllGenes")


#Gene Enrichment Analysis

enrichmentTARGET=GOenrichment(MCFSFeatures,colnames(TARGET_GE_Classifier)[1:length(TARGET_GE_Classifier)-1],"ENSEMBL")
TARGET_ReleventBP=findReleventBP(enrichmentTARGET,append(MCFSFeatures[1:20],"decision"),TARGET_GE_Classifier,1.0)

BP=unique(sapply(as.data.frame(enrichmentTARGET)$Description, function(x) x))

plotEnrichment(enrichmentTARGET,"Enrichment TARGET cohort","ResultsRules/AllGenes/.csv")
plotGeneRulesEnrichment(TARGET_ReleventBP,"ResultsRules/","2019-01-07","AllGenes")




plotEnrichment(enrichmentTARGET,"Enrichment TARGET All Genes","recalculatedResultRosettaMCFS","ResultsRules/AllGenes/.csv")
#--------------------------------------------------------------------------------------------------
#-------------------------Classifier On Metadata-------------------------------------------------
i2=match(rownames(TARGET_GE_Classifier_MCFS),rownames(metadata_for_exploratory))
TARGET_GE_Classifier_Meta=metadata_for_exploratory[i2,]
TARGET_GE_Classifier_Meta=TARGET_GE_Classifier_Meta[,c("Protocol","FAB.Category","t.8.21.","inv.16.","MLL","FLT3.ITD.positive.","NPM.mutation","CEBPA.mutation","WT1.mutation","CR.status.at.end.of.course.1")]

TARGET_GE_Classifier_Meta=cbind.data.frame(TARGET_GE_Classifier_Meta,decision)

revalue(TARGET_GE_Classifier_Meta$Protocol, c("AAML03P1" = 1,"AAML0531"=2,"CCG-2961"=3)) -> TARGET_GE_Classifier_Meta$Protocol

#,ruleFiltration=TRUE, ruleFiltrAccuracy=c(0,0.8),ruleFiltrSupport=c(1,3))

revalue(TARGET_GE_Classifier_Meta$Protocol, c("AAML03P1" = 1,"AAML0531"=2,"CCG-2961"=3)) -> TARGET_GE_Classifier_Meta$Protocol
revalue(TARGET_GE_Classifier_Meta$FAB.Category, c("." = 0,"M0"=1,"M1"=2,"M2"=3,"M4"=4,"M5"=5,"M6"=6,"M7"=7,"NOS"=8,"Unknown"=9)) -> TARGET_GE_Classifier_Meta$FAB.Category
revalue(TARGET_GE_Classifier_Meta$t.8.21., c("No" = 1,"Yes"=2,"Unknown"=3)) -> TARGET_GE_Classifier_Meta$t.8.21.
revalue(TARGET_GE_Classifier_Meta$inv.16., c("No" = 1,"Yes"=2,"Unknown"=3)) -> TARGET_GE_Classifier_Meta$inv.16.
revalue(TARGET_GE_Classifier_Meta$MLL, c("No" = 1,"Yes"=2,"Unknown"=3)) -> TARGET_GE_Classifier_Meta$MLL
revalue(TARGET_GE_Classifier_Meta$FLT3.ITD.positive., c("No" = 1,"NO"= 1,"Yes"=2,"YES"=2,"Unknown"=3)) -> TARGET_GE_Classifier_Meta$FLT3.ITD.positive.
revalue(TARGET_GE_Classifier_Meta$NPM.mutation, c("No" = 1,"NO"= 1,"Yes"=2,"YES"=2,"Unknown"=3)) -> TARGET_GE_Classifier_Meta$NPM.mutation
revalue(TARGET_GE_Classifier_Meta$CEBPA.mutation, c("No" = 1,"NO"= 1,"Yes"=2,"YES"=2,"Unknown"=3)) -> TARGET_GE_Classifier_Meta$CEBPA.mutation
revalue(TARGET_GE_Classifier_Meta$WT1.mutation, c("No" = 1,"NO"= 1,"Yes"=2,"YES"=2,"Unknown"=3)) -> TARGET_GE_Classifier_Meta$WT1.mutation
revalue(TARGET_GE_Classifier_Meta$CR.status.at.end.of.course.1, c("CR" = 1,"Death" = 2,"Not in CR" =3,"Unevaluable"=4)) -> TARGET_GE_Classifier_Meta$CR.status.at.end.of.course.1

TARGET_GE_Classifier_Meta=as.data.frame(apply(TARGET_GE_Classifier_Meta[,-dim(TARGET_GE_Classifier_Meta)[2]-1],2,as.numeric))
TARGET_GE_Classifier_Meta$decision=decision

resultRosettaMeta=rosetta(TARGET_GE_Classifier_Meta,discrete=TRUE)