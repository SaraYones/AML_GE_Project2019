library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)


#All Genes
Integrated_Matched_filter=getAnnotatedGenes(genes,colnames(Integrated_Matched_Corrected))
#All Genes
Integrated_Matched_Classifier=cbind.data.frame(Integrated_Matched_Corrected,decision_Integrated_Matched)
#MCFS and Boruta are run on ULAM
write.table(Integrated_Matched_Classifier,"MCFS_AML_GE_Classifier",row.names=TRUE)
Integrated_Matched_MCFSFeatures=FilterFeatures("Integrated_Cohort_Results/MCFSoutput/Integratedmatched1/Integratedmatched1__RI.csv",1000)
names(Integrated_Matched_Classifier)[names(Integrated_Matched_Classifier)=="decision_Integrated_Matched"]="decision"
Integrated_Matched_Accuracies=compareAccuracies(Integrated_Matched_Classifier,200,as.character(Integrated_Matched_MCFSFeatures))

IntegratedgeneLS=list(TARGET=colnames(TARGET_GE_Classifier_MCFS)[1:dim(TARGET_GE_Classifier_MCFS)[2]-1],Linda=colnames(Linda_GE_Classifier_MCFS)[1:dim(Linda_GE_Classifier_MCFS)[2]-1],Integrated=colnames(Integrated_Matched_Classifier_MCFS)[1:dim(Integrated_Matched_Classifier_MCFS)[2]-1])
Integratedintersection=plotVenn(IntegratedgeneLS,list(c("darkmagenta", "darkblue","red"),c(0.5,0.5,0.5),2,2, c("TARGET","Linda","Integrated")))



#Linda_MCFSFeatures=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindamatched2/Lindamatched1__RI.csv",10)
Integrated_Matched_filteredMCFS=getAnnotatedGenes(genes,Integrated_Matched_MCFSFeatures[1:20])
#Build a classifier with the intersection genes between TARGET and Linda Cohort
Linda_filteredMCFS=getAnnotatedGenes(genes,intersection$`MCFS:Boruta`)
Linda_MCFSFeatures=intersection$`MCFS:Boruta`
#All Genes MCFS
Integrated_Matched_classifierMCFSgenes=append(Integrated_Matched_filteredMCFS[[1]]$geneID,Integrated_Matched_filteredMCFS[[2]])

Linda_classifierMCFSgenes=intersection_annotated

#Ordering according to MCFS Features
Integrated_Matched_classifierMCFSgenes=Integrated_Matched_classifierMCFSgenes[!is.na(match(Integrated_Matched_MCFSFeatures, Integrated_Matched_classifierMCFSgenes))]

plotDistributionCategory(list(Integrated_Matched_filteredMCFS[[1]]$gene_type),list(c("MCFS")),"Distribution_Classes_FS",paste(getwd(),"/Integrated_Cohort_Results/2019-04-02/AllGenes/",sep=""),"")



Integrated_Matched_annotation=annotateInOrder(Integrated_Matched_filteredMCFS,Integrated_Matched_classifierMCFSgenes)

Integrated_Matched_Classifier_MCFS=Integrated_Matched_Classifier[,append(Integrated_Matched_classifierMCFSgenes,"decision")]

colnames(Integrated_Matched_Classifier_MCFS)<-append(Integrated_Matched_annotation,"decision")


plotNames_Integrated=gsub("(TARGET-20-((.)+-(R|T)))","\\2",rownames(Integrated_Matched_Classifier_MCFS))

rownames(Integrated_Matched_Classifier_MCFS)=plotNames_Integrated
Integrated_Matched_resultRosettaMCFS=rosetta(Integrated_Matched_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3)

Integrated_Matched_resultRosettaMCFSGenetic=rosetta(Integrated_Matched_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",reducer="Genetic",ruleFiltration=TRUE, discreteParam=3)

Integrated_Matched_recalculatedResultRosettaMCFS=recalculateRules(Integrated_Matched_Classifier_MCFS,Integrated_Matched_resultRosettaMCFS$main)

Integrated_Matched_recalculatedResultRosettaMCFSGenetic=recalculateRules(Integrated_Matched_Classifier_MCFS,Integrated_Matched_resultRosettaMCFSGenetic$main)

Integrated_Matched_clusteredRulesMCFS=clusterRules(Integrated_Matched_recalculatedResultRosettaMCFS,rownames(Integrated_Matched_Classifier_MCFS))

Integrated_Matched_clusteredRulesMCFSGenetic=clusterRules(Integrated_Matched_recalculatedResultRosettaMCFSGenetic,rownames(Integrated_Matched_Classifier_MCFS))


writeOutput("Integrated_Cohort_Results/",Sys.Date(),"AllGenes","Genetic",Integrated_Matched_clusteredRulesMCFSGenetic,Integrated_Matched_recalculatedResultRosettaMCFSGenetic,Integrated_Matched_resultRosettaMCFSGenetic$main,
            enrichmentIntegrated_Matched,"enrichment",FALSE)

writeOutput("Integrated_Cohort_Results/",Sys.Date(),"AllGenes","Johnson",Integrated_Matched_clusteredRulesMCFS,Integrated_Matched_recalculatedResultRosettaMCFS,Integrated_Matched_resultRosettaMCFS$main,
            enrichmentIntegrated_Matched,"enrichment",FALSE)

my.plots=vector(1, mode='list');
#svg(paste(temp,"/HeatMapAllGenes-",Sys.Date(),".svg",sep=""))
#pdf(paste("Integrated_Cohort_Results/AllGenes","/HeatMapAllGenes-",Sys.Date(),".pdf",sep=""), onefile=TRUE)

Integrated_Matched_clusters=heatmap.F(t(Integrated_Matched_clusteredRulesMCFS), colors=c('white','blue'),distmethod='pearson')
#write.csv(clusters,paste(temp,"/Clusters-",Sys.Date(),".csv",sep = ""))


plotNames_Integrated=gsub("(TARGET-20-((.)+-(R|T)))","\\2",rownames(Integrated_Matched_Classifier))


#Enrichment Analysis and finding which features from Rosetta is enriched from the significant biological processes

enrichmentIntegrated_Matched=GOenrichment(Integrated_Matched_MCFSFeatures,colnames(Integrated_Matched_Classifier)[1:length(Integrated_Matched_Classifier)-1],"ENSEMBL")

Integrated_Matched_RosettaEnrichment=findReleventBP(enrichmentIntegrated_Matched,append(Integrated_Matched_MCFSFeatures[1:20],"decision"),Integrated_Matched_Classifier,1.0)

plotGeneRulesEnrichment(Integrated_Matched_RosettaEnrichment,"Integrated_Cohort_Results/","2019-02-04","AllGenes")


plotEnrichment(enrichmentIntegrated_Matched,"Enrichment Integrated cohort All Genes","ResultsRules/AllGenes/.csv")

