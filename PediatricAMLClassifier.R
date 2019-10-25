library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)


#All Genes
Integrated_Matched_filter=getAnnotatedGenes(genes,colnames(Integrated_Matched))
#All Genes
Linda_GE_Classifier=cbind.data.frame(Linda_GE_Classifier,decision_Linda)
#Trying MCFS after normalizing for the whole cohort
Linda_GE_Classifier=cbind.data.frame(Linda_GE_Classifier[index,],decision_LindaAll[index])
names(Linda_GE_Classifier)[names(Linda_GE_Classifier)=="decision_LindaAll[index]"]="decision"
write.table(Linda_GE_Classifier_trial,"MCFS_AML_GE_Classifier",row.names=TRUE)
#MCFS and Boruta are run on ULAM
write.table(Linda_GE_Classifier,"MCFS_AML_GE_Classifier",row.names=TRUE)
Target_Genes=colnames(TARGET_GE_Classifier[,1:length(TARGET_GE_Classifier)-1])
Linda_Genes=colnames(Linda_GE_Classifier[,1: length(Linda_GE_Classifier)-1])
Target_Genes=Target_Genes[which(!(Target_Genes %in% Linda_Genes))]
Linda_Genes=Linda_Genes[which(!(Linda_Genes %in% Target_Genes ))]
Target_Genes=getAnnotatedGenes(genes,Target_Genes)
#Lindamatched3 -->Unmatched cases 46
#Lindamatched2 -->Matched cases 36
#Lindamatched4 --> Unmatched cases 48
#Lindamatched6-> Unmatched cases 36 (after Normalizing the whole pediatric cohort)

Linda_MCFSFeatures=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindamatched3/Lindamatched3__RI.csv",1000)

Linda_MCFSFeatures=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindaunmatched6/Lindaunmatched6__RI.csv",1000)
names(Linda_GE_Classifier)[names(Linda_GE_Classifier)=="decision_Linda"]="decision"
Linda_Accuracies=compareAccuracies(Linda_GE_Classifier,200,as.character(Linda_MCFSFeatures))

#Linda_MCFSFeatures=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindamatched2/Lindamatched1__RI.csv",10)

#Linda_MCFSFeatures=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindamatched2/Lindamatched1__RI.csv",10)

Linda_filteredMCFS=getAnnotatedGenes(genes,Linda_MCFSFeatures[1:10])
#Build a classifier with the intersection genes between TARGET and Linda Cohort
Linda_filteredMCFS=getAnnotatedGenes(genes,intersection$`MCFS:Boruta`)
Linda_MCFSFeatures=intersection$`MCFS:Boruta`
#All Genes MCFS
Linda_classifierMCFSgenes=append(Linda_filteredMCFS[[1]]$geneID,Linda_filteredMCFS[[2]])

Linda_classifierMCFSgenes=intersection_annotated

#Ordering according to MCFS Features
Linda_classifierMCFSgenes=Linda_classifierMCFSgenes[!is.na(match(Linda_MCFSFeatures, Linda_classifierMCFSgenes))]


Linda_annotation=annotateInOrder(Linda_filteredMCFS,Linda_classifierMCFSgenes)

Linda_GE_Classifier_MCFS=Linda_GE_Classifier[,append(Linda_classifierMCFSgenes,"decision")]

colnames(Linda_GE_Classifier_MCFS)<-append(Linda_annotation,"decision")

Linda_resultRosettaMCFSGenetic=rosetta(Linda_GE_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5,reducer="Genetic",ruleFiltration=TRUE, discreteParam=3)

Linda_recalculatedResultRosettaMCFSGenetic=recalculateRules(Linda_GE_Classifier_MCFS,Linda_resultRosettaMCFSGenetic$main)


Linda_resultRosettaMCFS=rosetta(Linda_GE_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5, discreteParam=3)

Linda_recalculatedResultRosettaMCFS=recalculateRules(Linda_GE_Classifier_MCFS,Linda_resultRosettaMCFS$main)


Linda_clusteredRulesMCFS=clusterRules(Linda_recalculatedResultRosettaMCFS,rownames(Linda_GE_Classifier_MCFS))

Linda_clusteredRulesMCFSGenetic=clusterRules(Linda_recalculatedResultRosettaMCFSGenetic,rownames(Linda_GE_Classifier_MCFS))

writeOutput("Linda_Cohort_Results/",Sys.Date(),"AllGenes","Genetic",Linda_clusteredRulesMCFSGenetic,Linda_recalculatedResultRosettaMCFSGenetic,Linda_resultRosettaMCFSGenetic$main,enrichmentLinda,"enrichment",FALSE)

writeOutput("Linda_Cohort_Results/",Sys.Date(),"AllGenes","Johnson",Linda_clusteredRulesMCFS,Linda_recalculatedResultRosettaMCFS,Linda_resultRosettaMCFS$main,enrichmentLinda,"enrichment",FALSE)



writeOutput("Integrated_Cohort_Results/",Sys.Date(),"AllGenes","Genetic",Integrated_Matched_clusteredRulesMCFSGenetic,Integrated_Matched_recalculatedResultRosettaMCFSGenetic,Integrated_Matched_resultRosettaMCFSGenetic$main,
            enrichmentLinda,"enrichment",FALSE)

writeOutput("Integrated_Cohort_Results/",Sys.Date(),"AllGenes","Johnson",Integrated_Matched_clusteredRulesMCFS,Integrated_Matched_recalculatedResultRosettaMCFS,Integrated_Matched_resultRosettaMCFS$main,
            enrichmentLinda,"enrichment",FALSE)


plotDistributionCategory(list(Linda_filteredMCFS[[1]]$gene_type),list(c("MCFS")),"Distribution_Classes_FS",paste(getwd(),"/Linda_Cohort_Results/2019-04-01/AllGenes/",sep=""),"")


#Enrichment Analysis and finding which features from Rosetta is enriched from the significant biological processes

enrichmentLinda=GOenrichment(Linda_MCFSFeatures,colnames(Linda_GE_Classifier)[1:length(Linda_GE_Classifier)-1],"ENSEMBL")

Linda_RosettaEnrichment=findReleventBP(enrichmentLinda,append(Linda_MCFSFeatures[1:50],"decision"),Linda_GE_Classifier,1.0)

BP=unique(sapply(as.data.frame(enrichmentLinda)$Description, function(x) x))

plotEnrichment(enrichmentLinda,"Enrichment Pediatric cohort All Genes","/Linda_Cohort_Results/2019-04-01/AllGenes/.csv")
plotGeneRulesEnrichment(Linda_RosettaEnrichment,"Linda_Cohort_Results/","2019-04-02","AllGenes")

#---------------------------Trying to predict the 36 samples with the TARGET_Classifier--------------------------------------

Linda_GE_Classifier_temp=Linda_GE_Classifier[,c(MCFSFeatures[1:20],"decision_Linda")]
filteredMCFS_Linda=getAnnotatedGenes(genes,colnames(Linda_GE_Classifier_temp)[1:20])
annotation_Linda=annotateInOrder(filteredMCFS_Linda,classifierMCFSgenes)
colnames(Linda_GE_Classifier_temp)=append(annotation_Linda,"decision") 
predictClass(Linda_GE_Classifier_temp, resultRosettaMCFS$main)
#------------------------------#What are the common genes between TARGET and Linda_Classifier---------------------------------------------------------------------------
geneLS=list(MCFS=MCFSFeatures,Boruta=Linda_MCFSFeatures)
intersection=plotVenn(geneLS,list(c("darkmagenta", "darkblue"),c(0.5,0.5),2,2, c("TARGET-MCFS","LC-MCFS")))
intersection_annotated=getAnnotatedGenes(genes,intersection$`MCFS:Boruta`)

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
