#Integrated Adult Cohort processing
#You have to run Normalization adults first

library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)
#Protein coding genes and those that are unannotated (To run feature selection only on protein coding genes)
filter=getAnnotatedGenes(genes,colnames(IntegratedAdultsClassifiertemp))
protein_coding=filter[[1]][which(filter[[1]]$gene_type =="protein_coding"),geneID]
non_pseudo_genes=filter[[1]][which(!(grepl("*\\_pseudogene",filter[[1]]$gene_type))),]$geneID
#-----------------------------------------------------------------------------------
IntegratedAdultsClassifiertemp=as.data.frame(cbind(IntegratedAdultsClassifier,decisionIntegratedAdultsDF))
#IntegratedAdultsClassifiertemp_PC=cbind.data.frame(IntegratedAdultsClassifiertemp[,append(protein_coding,filter[[2]])],decisionIntegratedAdultsDF)

IntegratedAdultsClassifiertemp=apply(IntegratedAdultsClassifiertemp[,1:dim(IntegratedAdultsClassifiertemp)[2]-1],2,as.numeric)
IntegratedAdultsClassifiertemp=as.data.frame(IntegratedAdultsClassifiertemp)
IntegratedAdultsClassifiertemp$decision=as.character(decisionIntegratedAdultsDF)
rownames(IntegratedAdultsClassifiertemp)=rownames(IntegratedAdultsClassifier)
#Splitting both data after normalization and batch correction then running the models again
index=which(rownames(IntegratedAdultsClassifiertemp) %in% sampleUsageAdults)
Linda_GE_classifer2AN= IntegratedAdultsClassifiertemp[which(rownames(IntegratedAdultsClassifiertemp) %in% sampleUsageAdults),]
Lietal_countDataAN=IntegratedAdultsClassifiertemp[which(rownames(IntegratedAdultsClassifiertemp) %in% rownames(lietalcountDataBN)),]
#two pass filter After Integration then splitting
Linda_GE_classifer2AN=Linda_GE_classifer2AN[,append(protein_coding,"decision")]
temp=Linda_GE_Classifer2BN[str_replace_all(sampleUsageAdults, "-", "\\."),which(colnames(Linda_GE_Classifer2BN) %in% colnames(Linda_GE_classifer2AN))]
result_genes=normalizeGE(as.matrix(temp[,protein_coding]),as.factor(Linda_GE_classifer2AN$decision),FALSE,FALSE,FALSE,TRUE)

#temp=Linda_GE_classifer2AN[,-result_genes$remove_genes]
#temp=temp[,append(result_genes$keep_genes,"decision")]

if(length(result_genes$remove_genes)!=0)
{temp=Linda_GE_classifer2AN[,-result_genes$remove_genes]
}else{
        temp=Linda_GE_classifer2AN
}
#temp=Linda_GE_classifer2AN[,-(append(normalizeGE(temp,as.factor(Linda_GE_classifer2AN$decision),FALSE,FALSE,FALSE,TRUE),dim(Linda_GE_classifer2AN)[2]))]
temp=temp[,append(result_genes$keep_genes,"decision")]
Linda_GE_classifer2AN=temp
Linda_GE_classifer2AN$decision=as.character(decisionIntegratedAdultsDF[index])
#----------------------For Lietal-------------------------------------------
temp=colnames(Lietal_countDataAN)[1:length(colnames(Lietal_countDataAN))-1]
#proteincoding
result_genes=normalizeGE(as.matrix(lietalcountDataBN[,protein_coding]),as.factor(Lietal_countDataAN$decision),FALSE,FALSE,FALSE,TRUE)
if(length(result_genes$remove_genes)!=0)
{temp=Lietal_countDataAN[,-result_genes$remove_genes]
}else{
        temp=Lietal_countDataAN
}
temp=temp[,append(result_genes$keep_genes,"decision")]
Lietal_countDataAN=temp
Lietal_countDataAN$decision=as.character(decisionLietal)

#--------------------------------------------------------------------

write.table(Linda_GE_classifer2AN,"MCFS_AML_Adults_GE_Classifier",row.names=TRUE)
write.table(Lietal_countDataAN,"MCFS_AML_Adults_GE_Classifier",row.names=TRUE)
write.table(IntegratedAdultsClassifiertemp,"MCFS_AML_Adults_GE_Classifier",row.names=TRUE)


# IntegratedAdultsClassifierJohnson <- Classifier(classifier = IntegratedAdultsClassifiertemp,flagAccuracy="Johnson",
#                 MCFSFeatures=FilterFeatures("Adults_Results/Integrated_Adults_Results/IntegratedAdultsMCFS/IntegratedAdults__RI.csv",1000)
# 
#                 )
IntegratedAdultsClassifier<- Classifier(classifier = IntegratedAdultsClassifiertemp,flagAccuracy="Johnson",path="/Adults_Results/Integrated_Adults_Results",
                                                MCFSFeatures=FilterFeatures("Adults_Results/Integrated_Adults_Results/IntegratedAdultsMCFS/IntegratedAdults1/IntegratedAdults1__RI.csv",1000)
                                                
)
IntegratedAdultsClassifierGenetic<- Classifier(classifier = IntegratedAdultsClassifiertemp,flagAccuracy="Genetic",path="/Adults_Results/Integrated_Adults_Results/Genetic",
                                        MCFSFeatures=FilterFeatures("Adults_Results/Integrated_Adults_Results/IntegratedAdultsMCFS/IntegratedAdults1/IntegratedAdults1__RI.csv",1000)
                                        
)

Linda_Adults_Johnson<- Classifier(classifier = Linda_GE_classifer2AN,flagAccuracy="Johnson",path="/Adults_Results/Linda_Adults_Cohort_Results/AfterIntegration",
                                               MCFSFeatures=FilterFeatures("Adults_Results/Linda_Adults_Cohort_Results/AfterIntegration/LocalAN1/LocalAN1__RI.csv",1000)[[1]]
                                               
)
Linda_Adults_Johnson_PC<- Classifier(classifier = Linda_GE_classifer2AN,flagAccuracy="Johnson",path="/Adults_Results/Linda_Adults_Cohort_Results/AfterIntegration",
                                  MCFSFeatures=FilterFeatures("Adults_Results/Linda_Adults_Cohort_Results/AfterIntegration/Linda_Adults_AN_PC-5All/Linda_Adults_AN_PC-5All__RI.csv",1000)[[1]]
                                  
,keyType="ENSEMBL",underSample=TRUE)


Lietal_Adults_Johnson<- Classifier(classifier = Lietal_countDataAN,flagAccuracy="Johnson",path="/Adults_Results/Lietal_Cohort_Results/AfterIntegration",
                                  MCFSFeatures=FilterFeatures("Adults_Results/Lietal_Cohort_Results/AfterIntegration/LietalAN1/LietalAN1__RI.csv",1000)[[1]]
                                  
)
Lietal_Adults_Johnson_PC<- Classifier(classifier = Lietal_countDataAN,flagAccuracy="Johnson",path="/Adults_Results/Lietal_Cohort_Results/AfterIntegration",
                                   MCFSFeatures=FilterFeatures("Adults_Results/Lietal_Cohort_Results/AfterIntegration/LietalANPC-5/LietalANPC-5__RI.csv",1000)[[1]]
                                   
)


#If i update in the class
IntegratedAdultsClassifierGenetic<- IntegratedAdultsClassifierGenetic$copy()
IntegratedAdultsClassifier<- IntegratedAdultsClassifier$copy()
Linda_Adults_Johnson<-Linda_Adults_Johnson$copy()
Linda_Adults_Johnson<-Linda_Adults_Johnson$copy()
Linda_Adults_Johnson_PC<-Linda_Adults_Johnson_PC$copy()
Lietal_Adults_Johnson<-Lietal_Adults_Johnson$copy()
Lietal_Adults_Johnson_PC<-Lietal_Adults_Johnson_PC$copy()
# IntegratedAdultsClassifier$flagAccuracy="Johnson"
 # IntegratedAdultsClassifier$AccuraciesGenetic=tempAccuracyGenetic
 # IntegratedAdultsClassifier$Accuracies=tempAccuracyJohnson
 # 
 # IntegratedAdultsClassifierGenetic$AccuraciesGenetic=tempAccuracyGenetic
 # IntegratedAdultsClassifierGenetic$Accuracies=tempAccuracyJohnson
 #
IntegratedAdultsClassifier$findAccuracies()
IntegratedAdultsClassifierGenetic$findAccuracies()
Linda_Adults_Johnson$findAccuracies()
Linda_Adults_Johnson_PC$findAccuracies()
Lietal_Adults_Johnson$findAccuracies()
Lietal_Adults_Johnson_PC$findAccuracies()
Linda_Adults_Johnson$createAnnotatedClassifier()
Linda_Adults_Johnson_PC$createAnnotatedClassifier()
Lietal_Adults_Johnson$createAnnotatedClassifier()
Lietal_Adults_Johnson_PC$createAnnotatedClassifier()
IntegratedAdultsClassifier$createAnnotatedClassifier()
IntegratedAdultsClassifierGenetic$createAnnotatedClassifier()
# IntegratedAdultsClassifierGenetic$createAnnotatedClassifier()
IntegratedAdultsClassifier$clusterRulesandWriteoutput()
Linda_Adults_Johnson$clusterRulesandWriteoutput()
Linda_Adults_Johnson_PC$clusterRulesandWriteoutput()
Linda_Adults_Johnson_PC$clusterRulesandWriteoutput()
Linda_Adults_Johnson_PC$computeEnrichment()
Lietal_Adults_Johnson$clusterRulesandWriteoutput()
Lietal_Adults_Johnson_PC$clusterRulesandWriteoutput()
IntegratedAdultsClassifierGenetic$clusterRulesandWriteoutput()
IntegratedAdultsClassifierGenetic$computeEnrichment()

# tempAccuracyGenetic=IntegratedAdultsClassifier$AccuraciesGenetic
# tempAccuracyJohnson=IntegratedAdultsClassifier$Accuracies

#Run VisuNet
#Remember to click DONE once you finish working on VisuNet
visunet(Linda_Adults_Johnson_PC$resultRosettaMCFSJohnson$main, type = "RDF")
visunet(Linda_Adults_Johnson_PC$resultRosettaMCFSGenetic$main, type = "RDF")

 

plot(Linda_Adult_Accuracies,type='l',xaxt="n",main="Accuracies")
axis(1,at=seq(1,as.numeric(200/10), by = 1),labels=as.character(seq(10, 200, by = 10)))

