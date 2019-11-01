#Read Lietal data without normalization
#Read Linda cohort2 data before normalization
#Integrate whole cohort 
#remove contorls and then run PCA and remove batches

Linda_GE_Classifer2BN=Linda_GE_Classifier2
rownames(lietalcountData)<-gsub("(.*)\\.(.*)","\\1",rownames(lietalcountData))
lietalcountDataBN=t(lietalcountData)
Common_GenesAdults=colnames(Linda_GE_Classifer2BN)[which(colnames(Linda_GE_Classifer2BN) %in% colnames(lietalcountDataBN))]
IntegratedAdultsClassifier=rbind(Linda_GE_Classifer2BN[,Common_GenesAdults],lietalcountDataBN[,Common_GenesAdults])
decisionIntegratedAdults=append(decision_Linda2,decisionLietal)

#Normalize the whole cohort first then 

IntegratedAdultsClassifier=normalizeGE(IntegratedAdultsClassifier,as.factor(decisionIntegratedAdults),FALSE,TRUE,TRUE)
#Plot PCA before removal of batch 
plotPCA(getwd(),IntegratedAdultsClassifier,decisionIntegratedAdults,"PCA_Adults_Integrated")
#Create grouping for batches 
batchesAdultsIntegrated=rep("batch1",dim(Linda_GE_Classifer2BN[,Common_GenesAdults])[1])
batchesAdultsIntegrated=append(batchesAdultsIntegrated,rep("batch2",dim(lietalcountDataBN[,Common_GenesAdults])[1]))
#Correcting batches
IntegratedAdultsCorrected=removeBatchEffect(IntegratedAdultsClassifier,as.factor(batchesAdultsIntegrated),batchesAdultsIntegrated,1)
IntegratedAdultsCorrected=t(IntegratedAdultsCorrected)
plotPCA(getwd(),IntegratedAdultsCorrected,decisionIntegratedAdults,"PCA_Adults_Integrated")
IntegratedAdultsClassifier=IntegratedAdultsCorrected
#Split the cohort to only include R1 and D 

cases_Adults_Relapse1_unmatchedtrial<-read.xlsx("Linda_Cohort_Metadata/AdultsData/AdultsMetadata.xlsx", sheetName = "Relapse1")
#The coloumn names are marked with "." so we need to change to -
new_col_names=str_replace_all(rownames(IntegratedAdultsClassifier), "\\.","-")
rownames(IntegratedAdultsClassifier)<-new_col_names

temp1=cases_Adults_Relapse1_unmatched$AdultDiagnosis[!is.na(cases_Adults_Relapse1_unmatched$AdultDiagnosis)]
temp2=cases_Adults_Relapse1_unmatched$AdultRelapse[!is.na(cases_Adults_Relapse1_unmatched$AdultRelapse)]
cases_Adults_Relapse1_unmatched<-NULL

cases_Adults_Relapse1_unmatched$AdultDiagnosis=temp1
cases_Adults_Relapse1_unmatched$AdultRelapse=temp2

#cases_Adults_Relapse1_unmatched$AdultDiagnosis=str_replace_all(cases_Adults_Relapse1_unmatched$AdultDiagnosis , "_|-", ".")
#cases_Adults_Relapse1_unmatched$AdultRelapse=str_replace_all(cases_Adults_Relapse1_unmatched$AdultRelapse , "_|-", ".")


#Filter from the whole cohort only the cases of Diagnosis and Relapse 1 that are adults cases
filterCases=which(rownames(IntegratedAdultsClassifier) %in% append(append(levels(cases_Adults_Relapse1_unmatched$AdultDiagnosis),levels(cases_Adults_Relapse1_unmatched$AdultRelapse)),rownames(lietalcountDataBN) ))
IntegratedAdultsClassifier=IntegratedAdultsClassifier[filterCases,]
decisionIntegratedAdults=decisionIntegratedAdults[filterCases]


#because we normalized on the whole dataset so i decided not to remove until everything was normalized
remove_cols=nearZeroVar(IntegratedAdultsClassifier)
if(length(remove_cols)!=0)
  IntegratedAdultsClassifier=IntegratedAdultsClassifier[,-remove_cols]


#Annotate the whole decision table

IntegratedAdultsAnnotated=as.data.frame(annotateDecisionTable(IntegratedAdultsCorrected ,genes))


write.table(as.data.frame(IntegratedAdultsAnnotated),"Adults_Results/Integrated_Adults_Results/IntegratedAdultsNormalized.csv",sep = ",",row.names = T,col.names = T)

