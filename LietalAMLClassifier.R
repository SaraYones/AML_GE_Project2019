#Lietal Classifier

#Adult Cohort processing
#You have to run Normalization adults first

library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)

#remove the . in the ENSEMBL ID
colnames(lietalcountData)<-gsub("(.*)\\.(.*)","\\1",colnames(lietalcountData))

write.table(lietalcountData,"MCFS_AML_Lietal_GE_Classifier",row.names=TRUE)


Lietal_MCFSFeatures=FilterFeatures("Linda_Adults_Cohort_Results/MCFSoutput/LietalmatchedAdults1/LietalmatchedAdults1__RI.csv",1000)
Lietal_Accuracies=compareAccuracies(lietalcountData,200,as.character(Lietal_MCFSFeatures))
Lietal_Accuracies_Genetic=compareAccuracies(lietalcountData,200,as.character(Lietal_MCFSFeatures))


Lietal_filtered_MCFS=getAnnotatedGenes(genes,Lietal_MCFSFeatures[1:30])
#All Genes MCFS
Lietal_classifierMCFSgenes=append(Lietal_filtered_MCFS[[1]]$geneID,Lietal_filtered_MCFS[[2]])


#Ordering according to MCFS Features
Lietal_classifierMCFSgenes=Lietal_classifierMCFSgenes[!is.na(match(Lietal_MCFSFeatures, Lietal_classifierMCFSgenes))]

Lietal_annotation=annotateInOrder(Lietal_filtered_MCFS,Lietal_classifierMCFSgenes)

Lietal_annotation=make.unique(Lietal_annotation)


Lietal_Classifier_MCFS=lietalcountData[,append(Lietal_classifierMCFSgenes,"decision")]

colnames(Lietal_Classifier_MCFS)<-append(Lietal_annotation,"decision")
rownames_Lietal=paste("Lietal_",gsub("((Patient_AML)_([0-9][0-9][0-9]-(DX|RX)))","\\3",rownames(Lietal_Classifier_MCFS)),sep="")

rownames(Lietal_Classifier_MCFS)<-NULL
rownames(Lietal_Classifier_MCFS)<-make.names(rownames_Lietal, unique=TRUE)

Lietal_resultRosettaMCFSGenetic=rosetta(Lietal_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=10,reducer="Genetic",ruleFiltration=TRUE, discreteParam=3,ruleFiltrSupport=c(0,3))

Lietal_recalculatedResultRosettaMCFSGenetic=recalculateRules(Lietal_Classifier_MCFS, Lietal_resultRosettaMCFSGenetic$main)

#Lietal_recalculatedResultRosettaMCFSGenetic=Lietal_recalculatedResultRosettaMCFSGenetic[Lietal_recalculatedResultRosettaMCFSGenetic$supportLHS!=0,]

Lietal_clusteredRulesMCFS_Genetic=clusterRules(Lietal_recalculatedResultRosettaMCFSGenetic,rownames(Lietal_Classifier_MCFS))

Lietal_resultRosettaMCFS=rosetta(Lietal_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=10, discreteParam=3)

Lietal_recalculatedResultRosettaMCFS=recalculateRules(Lietal_Classifier_MCFS,Lietal_resultRosettaMCFS$main)

Lietal_clusteredRulesMCFS=clusterRules(Lietal_recalculatedResultRosettaMCFS,rownames(Lietal_Classifier_MCFS))

writeOutput("Adults_Results/Lietal_Cohort_Results",Sys.Date(),"AllGenes","Genetic",append(Lietal_annotation,Lietal_filtered_MCFS[[2]]),Lietal_clusteredRulesMCFS_Genetic,Lietal_recalculatedResultRosettaMCFSGenetic,Lietal_resultRosettaMCFSGenetic$main,enrichmentLietal,"enrichment",FALSE)

writeOutput("Adults_Results/Lietal_Cohort_Results",Sys.Date(),"AllGenes","Johnson",append(Lietal_annotation,Lietal_filtered_MCFS[[2]]),Lietal_clusteredRulesMCFS,Lietal_recalculatedResultRosettaMCFS,Lietal_resultRosettaMCFS$main,enrichmentLietal,"enrichment",FALSE)

#------------------------------#What are the common genes between TARGET and Linda_Classifier---------------------------------------------------------------------------

intersection=reportCommonGenes("Adults_Results","500",Sys.Date(),list(Lietal=unique(append(Lietal_annotation,Lietal_filtered_MCFS[[2]])),Linda=unique(append(Linda_Adults_annotation,Linda_filtered_Adults_MCFS[[2]]))),list(c("darkmagenta", "darkblue"),c(0.5,0.5),2,2, c("Lietal","Linda_Adults"),alpha = c(0.5, 0.5), cat.pos = 3))


plotDistributionCategory(list(Lietal_filtered_MCFS[[1]]$gene_type),list(c("MCFS")),"Distribution_Classes_FS",paste(getwd(),"/Adults_Results/Lietal_Cohort_Results",sep=""),"")


#Enrichment Analysis and finding which features from Rosetta is enriched from the significant biological processes

enrichmentLietal=GOenrichment(Lietal_MCFSFeatures,colnames(lietalcountData)[1:length(lietalcountData)-1],"ENSEMBL")


Linda_Adults_RosettaEnrichment=findReleventBP(enrichmentLinda_Adults,append(Linda_Adult_MCFSFeatures[1:50],"decision"),Linda_GE_Classifier2,1.0)

BP=unique(sapply(as.data.frame(enrichmentLinda)$Description, function(x) x))

plotEnrichment(enrichmentLietal,"Enrichment Lietal cohort All Genes","/Linda_Cohort_Results/2019-04-01/AllGenes/.csv")
plotGeneRulesEnrichment(Linda_RosettaEnrichment,"Linda_Cohort_Results/","2019-04-02","AllGenes")

#---------------------------Trying to predict the 36 samples with the TARGET_Classifier--------------------------------------

Linda_GE_Classifier_temp=Linda_GE_Classifier[,c(MCFSFeatures[1:20],"decision_Linda")]
filteredMCFS_Linda=getAnnotatedGenes(genes,colnames(Linda_GE_Classifier_temp)[1:20])
annotation_Linda=annotateInOrder(filteredMCFS_Linda,classifierMCFSgenes)
colnames(Linda_GE_Classifier_temp)=append(annotation_Linda,"decision") 
predictClass(Linda_GE_Classifier_temp, resultRosettaMCFS$main)

#-------------------------------Predict Normals------------------------------------------
#You have to run normalizationPediatrics first
Normals=Linda_GE_Classifier[which(grepl("BM.*",rownames (Linda_GE_Classifier))),]
decision_Normals=rep("Diagnosis",dim(Normals)[1])
Normals=cbind(Normals,decision_Normals)
Linda_GE_Classifier_temp=Normals[,c(Linda_MCFSFeatures[1:10],"decision_Normals")]
filteredMCFS_Linda=getAnnotatedGenes(genes,colnames(Linda_GE_Classifier_temp)[1:10])
annotation_Linda=annotateInOrder(filteredMCFS_Linda,Linda_classifierMCFSgenes)
colnames(Linda_GE_Classifier_temp)=append(annotation_Linda,"decision") 
predictClass(Linda_GE_Classifier_temp, Linda_resultRosettaMCFSGenetic$main)

#Checking runs of Lindaunmatched6 and Lindamatched3
x=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindamatched3/Lindamatched3__RI.csv",50)
y=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindaunmatched6/Lindaunmatched6__RI.csv",50)
