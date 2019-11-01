#Integrated Adult Cohort processing
#You have to run Normalization adults first

library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)

IntegratedAdultsClassifiertemp=as.data.frame(cbind(IntegratedAdultsClassifier,decisionIntegratedAdults))
IntegratedAdultsClassifiertemp=as.data.frame(apply(IntegratedAdultsClassifiertemp[,1:dim(IntegratedAdultsClassifier)[2]-1],2,as.numeric))
IntegratedAdultsClassifiertemp$decision=as.character(decisionIntegratedAdults)
rownames(IntegratedAdultsClassifiertemp)=rownames(IntegratedAdultsClassifier)
write.table(IntegratedAdultsClassifier,"MCFS_AML_Adults_GE_Classifier",row.names=TRUE)


# IntegratedAdultsClassifierJohnson <- Classifier(classifier = IntegratedAdultsClassifiertemp,flagAccuracy="Johnson",
#                 MCFSFeatures=FilterFeatures("Adults_Results/Integrated_Adults_Results/IntegratedAdultsMCFS/IntegratedAdults__RI.csv",1000)
# 
#                 )
IntegratedAdultsClassifier<- Classifier(classifier = IntegratedAdultsClassifiertemp,flagAccuracy="Johnson",path="/Adults_Results/Integrated_Adults_Results",
                                                MCFSFeatures=FilterFeatures("Adults_Results/Integrated_Adults_Results/IntegratedAdultsMCFS/IntegratedAdults__RI.csv",1000)
                                                
)
IntegratedAdultsClassifierGenetic<- Classifier(classifier = IntegratedAdultsClassifiertemp,flagAccuracy="Genetic",path="/Adults_Results/Integrated_Adults_Results/Genetic",
                                        MCFSFeatures=FilterFeatures("Adults_Results/Integrated_Adults_Results/IntegratedAdultsMCFS/IntegratedAdults__RI.csv",1000)
                                        
)

# IntegratedAdultsClassifier$flagAccuracy="Johnson"
 IntegratedAdultsClassifier$AccuraciesGenetic=tempAccuracyGenetic
 IntegratedAdultsClassifier$Accuracies=tempAccuracyJohnson
 
 IntegratedAdultsClassifierGenetic$AccuraciesGenetic=tempAccuracyGenetic
 IntegratedAdultsClassifierGenetic$Accuracies=tempAccuracyJohnson
 #
##IntegratedAdultsClassifier$findAccuracies()
IntegratedAdultsClassifier$createAnnotatedClassifier()
IntegratedAdultsClassifierGenetic$createAnnotatedClassifier()
IntegratedAdultsClassifier$clusterRulesandWriteoutput()
IntegratedAdultsClassifierGenetic$clusterRulesandWriteoutput()
# tempAccuracyGenetic=IntegratedAdultsClassifier$AccuraciesGenetic
# tempAccuracyJohnson=IntegratedAdultsClassifier$Accuracies


Classifier <- setRefClass("Classifier",
              #          fields = list(classifier = "data.frame",MCFSFeatures = "list", Accuracies = "list",classifierMCFSgenes="data.frame",
              #             path="string",Accuracies_Genetic="list",classifierMCFSgenes="list",Classifier_MCFS="as.data.frame",resultRosettaMCFSGenetic="as.data.frame"
              # ),
              fields = list(classifier = "data.frame",flagAccuracy="character",MCFSFeatures = "character", Accuracies = "numeric",AccuraciesGenetic="numeric", path="character",classifierMCFSgenes="character",Classifier_MCFS="data.frame",annotation="character",resultRosettaMCFSGenetic="list",resultRosettaMCFSJohnson="list",resultRosettaMCFSMerged="list",recalculatedResultRosettaMCFSJohnson="data.frame",recalculatedResultRosettaMCFSGenetic="data.frame"),
                       methods =  list(findAccuracies = function() {
                         classifier <<- classifier
                        MCFSFeatures<<-MCFSFeatures
                        flagAccuracy<<-flagAccuracy
                        
                       if(flagAccuracy=="Genetic")
                        {
  
                     AccuraciesGenetic<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy)
                       }else if(flagAccuracy=="Johnson")
                        {
                     Sys.sleep(4)
                      Accuracies<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy)
                          
                          
                          
                     }
         
                      
                         
                       },
                    createAnnotatedClassifier = function() {
                      
                      Accuracies<<-Accuracies
                      AccuraciesGenetic<<-AccuraciesGenetic
                      classifierMCFSgenes<<-classifierMCFSgenes
                      Classifier_MCFS<<-Classifier_MCFS
                      resultRosettaMCFSJohnson<<- resultRosettaMCFSJohnson
                      resultRosettaMCFSGenetic<<- resultRosettaMCFSGenetic
                      recalculatedResultRosettaMCFSJohnson<<-recalculatedResultRosettaMCFSJohnson
                      recalculatedResultRosettaMCFSGenetic<<-recalculatedResultRosettaMCFSGenetic
                   #   recalculatedResultRosettaMCFS<<-recalculatedResultRosettaMCFS
                      classifier<<-classifier
                      flagAccuracy<<-flagAccuracy
                      resultRosettaMCFSMerged<<-resultRosettaMCFSMerged
                      genes<<-genes
                      if(flagAccuracy=="Johnson")
                      {
                      maxAccuracy<-which(Accuracies==max(Accuracies))*10
                      }else if(flagAccuracy=="Genetic") {
                        
                        maxAccuracy<-which(AccuraciesGenetic==max(AccuraciesGenetic))*10
                        
                      }
                        
                      filtered_MCFS<<-getAnnotatedGenes(genes,MCFSFeatures[1:maxAccuracy])
                      #All Genes MCFS
                      
                      classifierMCFSgenes<<-append(filtered_MCFS[[1]]$geneID,filtered_MCFS[[2]])
                      
                      
                      #Ordering according to MCFS Features
                     classifierMCFSgenes<<-classifierMCFSgenes[!is.na(match(MCFSFeatures, classifierMCFSgenes))]
                      
                      
                     annotation<<-annotateInOrder(filtered_MCFS,classifierMCFSgenes)
                      
                    annotation<<-make.unique(annotation)
                      
                      
                      Classifier_MCFS<<-classifier[,append(classifierMCFSgenes,"decision")]
                      
                      
                      colnames(Classifier_MCFS)<<-append(annotation,"decision")
                              
                      
                       resultRosettaMCFSGenetic<<-rosetta(Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=10,reducer="Genetic",ruleFiltration=TRUE,ruleFiltrSupport=c(1,3) , discreteParam=3)
                       Sys.sleep(10)
                      # resultRosettaMCFSJohnson<<-NULL
                      resultRosettaMCFSJohnson<<-rosetta(Classifier_MCFS)
                      resultRosettaMCFSMerged<<-list("Johnson"=resultRosettaMCFSJohnson,"Genetic"=resultRosettaMCFSGenetic)

                      recalculatedResultRosettaMCFSJohnson<<-recalculateRules(Classifier_MCFS,resultRosettaMCFSJohnson$main)
                      Sys.sleep(10)
                     # recalculatedResultRosettaMCFSGenetic<<-NULL
                     recalculatedResultRosettaMCFSGenetic<<-recalculateRules(Classifier_MCFS,resultRosettaMCFSGenetic$main)
                     resultRosettaMCFSMerged<<-list("JohnsonResult"=resultRosettaMCFSJohnson,"GeneticResult"=resultRosettaMCFSGenetic,"recalculatedResultJohnson"=recalculatedResultRosettaMCFSJohnson,"recalculatedResultGenetic"=recalculatedResultRosettaMCFSGenetic)
                     plotDistributionCategory(list(filtered_MCFS[[1]]$gene_type),list(c("MCFS")),"Distribution_Classes_FS",paste(getwd(),"/Adults_Results/Linda_Adults_Cohort_Results",sep=""),"")
                     
                    },
                   
                   
                   clusterRulesandWriteoutput = function() {
                     #recalculatedResultRosettaMCFS<<-recalculatedResultRosettaMCFS
                     resultRosettaMCFSGenetic<<-resultRosettaMCFSGenetic
                     recalculatedResultRosettaMCFSGenetic<<-recalculatedResultRosettaMCFSGenetic
                     resultRosettaMCFSJohnson<<-resultRosettaMCFSJohnson
                     recalculatedResultRosettaMCFSJohnson<<-recalculatedResultRosettaMCFSJohnson
                     path<<-path
                     
                    clusteredRulesMCFS<<-clusterRules(recalculatedResultRosettaMCFSJohnson,rownames(Classifier_MCFS))
                     
                    clusteredRulesMCFSGenetic<<-clusterRules(recalculatedResultRosettaMCFSGenetic,rownames(Classifier_MCFS))
                    writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Genetic",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFSGenetic,recalculatedResultRosettaMCFSGenetic,resultRosettaMCFSGenetic$main,enrichment = NULL,"enrichment",FALSE)
                    
                    writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Johnson",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFS, recalculatedResultRosettaMCFSJohnson,resultRosettaMCFSJohnson$main,enrichment = NULL,"enrichment",FALSE)
                    
                     
                   }
                     
                  
                     )
                         
                       )




plot(Linda_Adult_Accuracies,type='l',xaxt="n",main="Accuracies")
axis(1,at=seq(1,as.numeric(200/10), by = 1),labels=as.character(seq(10, 200, by = 10)))


Linda_Adults_clusteredRulesMCFS=clusterRules(Linda_Adults_recalculatedResultRosettaMCFS,rownames(Linda_GE_Adults_Classifier_MCFS))

Linda_Adults_clusteredRulesMCFSGenetic=clusterRules(Linda_Adults_recalculatedResultRosettaMCFSGenetic,rownames(Linda_GE_Adults_Classifier_MCFS))

writeOutput("Adults_Results/Linda_Adults_Cohort_Results/",Sys.Date(),"AllGenes","Genetic",append(Linda_Adults_annotation,Linda_filtered_Adults_MCFS[[2]]),Linda_Adults_clusteredRulesMCFSGenetic,Linda_Adults_recalculatedResultRosettaMCFSGenetic,Linda_Adults_resultRosettaMCFSGenetic$main,enrichmentLinda_Adults,"enrichment",FALSE)

writeOutput("Adults_Results/Linda_Adults_Cohort_Results/",Sys.Date(),"AllGenes","Johnson",append(Linda_Adults_annotation,Linda_filtered_Adults_MCFS[[2]]),Linda_Adults_clusteredRulesMCFS,Linda_Adults_recalculatedResultRosettaMCFS,Linda_Adults_resultRosettaMCFS$main,enrichmentLinda_Adults,"enrichment",FALSE)


plotDistributionCategory(list(Linda_filtered_Adults_MCFS[[1]]$gene_type),list(c("MCFS")),"Distribution_Classes_FS",paste(getwd(),"/Adults_Results/Linda_Adults_Cohort_Results",sep=""),"")


#Enrichment Analysis and finding which features from Rosetta is enriched from the significant biological processes

enrichmentLinda_Adults=GOenrichment(Linda_Adult_MCFSFeatures,colnames(Linda_GE_Classifier2)[1:length(Linda_GE_Classifier2)-1],"ENSEMBL")

Linda_Adults_RosettaEnrichment=findReleventBP(enrichmentLinda_Adults,append(Linda_Adult_MCFSFeatures[1:50],"decision"),Linda_GE_Classifier2,1.0)

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
