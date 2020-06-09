#Read Lietal data without normalization
#Read Linda cohort2 data before normalization
#Integrate whole cohort 
#remove contorls and then run PCA and remove batches
#Basicly this rely on Linda_GE_Classifier2BN and lietalcountDataBN which is saved in the latest workspace

Linda_GE_Classifer2BN=Linda_GE_Classifier2
rownames(lietalcountData)<-gsub("(.*)\\.(.*)","\\1",rownames(lietalcountData))
lietalcountDataBN<-t(t(t(lietalcountData)))
Common_GenesAdults=colnames(Linda_GE_Classifer2BN)[which(colnames(Linda_GE_Classifer2BN) %in% colnames(lietalcountDataBN))]
IntegratedAdultsClassifierDF=rbind(Linda_GE_Classifer2BN[,Common_GenesAdults],lietalcountDataBN[,Common_GenesAdults])
decisionIntegratedAdultsDF=append(decision_Linda2,decisionLietal)

#Normalize the whole cohort first then 

IntegratedAdultsClassifierDF=normalizeGE(IntegratedAdultsClassifierDF,as.factor(decisionIntegratedAdultsDF),TRUE,FALSE,TRUE,FALSE)

#Plot PCA before removal of batch 
plotPCA(getwd(),IntegratedAdultsClassifierDF,append(rep("Linda",dim(Linda_GE_Classifer2BN)[1]),rep("leital",dim(lietalcountDataBN)[1])),"PCA_Adults_Integrated")
#Create grouping for batches 
batchesAdultsIntegrated=rep("batch1",dim(Linda_GE_Classifer2BN[,Common_GenesAdults])[1])
#remove batches after filtering
batchesAdultsIntegrated=rep("batch1",dim(Linda_GE_Classifer2BN)[1])
batchesAdultsIntegrated=append(batchesAdultsIntegrated,rep("batch2",dim(lietalcountDataBN[,Common_GenesAdults])[1]))
#Correcting batches
IntegratedAdultsCorrected=removeBatchEffect(IntegratedAdultsClassifierDF,as.factor(batchesAdultsIntegrated),batchesAdultsIntegrated,1)
IntegratedAdultsCorrected=t(IntegratedAdultsCorrected)
plotPCA(getwd(),IntegratedAdultsCorrected,decisionIntegratedAdultsDF,"PCA_Adults_Integrated")
IntegratedAdultsClassifier=IntegratedAdultsCorrected
#Split the cohort to only include R1 and D 

cases_Adults_Relapse1_unmatchedtrial<-read.xlsx("Linda_Cohort_Metadata/AdultsData/AdultsMetadata.xlsx", sheetName = "Relapse1")
#The coloumn names are marked with "." so we need to change to -
new_col_names=str_replace_all(rownames(IntegratedAdultsClassifier), "\\.","-")
rownames(IntegratedAdultsClassifier)<-new_col_names
new_col_names=str_replace_all(rownames(IntegratedAdultsClassifier), "Patient_AML","lietal")
rownames(IntegratedAdultsClassifier)<-new_col_names
new_col_names=str_replace_all(rownames(lietalcountDataBN), "Patient_AML","lietal")
#rownames(IntegratedAdultsClassifier)<-new_col_names
#rownames(IntegratedAdultsClassifier)<-new_col_names
rownames(lietalcountDataBN)<-new_col_names

temp1=cases_Adults_Relapse1_unmatched$AdultDiagnosis[!is.na(cases_Adults_Relapse1_unmatched$AdultDiagnosis)]
temp2=cases_Adults_Relapse1_unmatched$AdultRelapse[!is.na(cases_Adults_Relapse1_unmatched$AdultRelapse)]
cases_Adults_Relapse1_unmatched<-NULL

cases_Adults_Relapse1_unmatched$AdultDiagnosis=temp1
cases_Adults_Relapse1_unmatched$AdultRelapse=temp2

#cases_Adults_Relapse1_unmatched$AdultDiagnosis=str_replace_all(cases_Adults_Relapse1_unmatched$AdultDiagnosis , "_|-", ".")
#cases_Adults_Relapse1_unmatched$AdultRelapse=str_replace_all(cases_Adults_Relapse1_unmatched$AdultRelapse , "_|-", ".")


#Filter from the whole cohort only the cases of Diagnosis and Relapse 1 that are adults cases
filterCases=which(rownames(IntegratedAdultsClassifier) %in% append(append(levels(cases_Adults_Relapse1_unmatched$AdultDiagnosis),levels(cases_Adults_Relapse1_unmatched$AdultRelapse)),rownames(lietalcountDataBN) ))
filterCases=which(rownames(IntegratedAdultsClassifier) %in% append(append(levels(cases_Adults_Relapse1_unmatched$AdultDiagnosis),levels(cases_Adults_Relapse1_unmatched$AdultRelapse)),rownames(lietalcountDataBN) ))
filterCases=which(rownames(IntegratedAdultsClassifier) %in% append(append(levels(cases_Adults_Relapse1_unmatched$AdultDiagnosis),levels(cases_Adults_Relapse1_unmatched$AdultRelapse)),rownames(lietalcountDataBN) ))
#Filter only the unpaired latest pure relapse
filterCases=which(rownames(IntegratedAdultsClassifier) %in% append(sampleUsageAdults,rownames(lietalcountDataBN)))
IntegratedAdultsClassifier=IntegratedAdultsClassifier[filterCases,]
decisionIntegratedAdultsDF=decisionIntegratedAdultsDF[filterCases]


#because we normalized on the whole dataset so i decided not to remove until everything was normalized
remove_cols=nearZeroVar(IntegratedAdultsClassifier)
if(length(remove_cols)!=0)
  IntegratedAdultsClassifier=IntegratedAdultsClassifier[,-remove_cols]


#Annotate the whole decision table
#IntegratedAdultsCorrected is all the genes 
#Integrated AdultsClassifier comes from IntegratedAdultsClassifierDF which is filtered for the latest pure relapse

IntegratedAdultsAnnotated=as.data.frame(annotateDecisionTable( IntegratedAdultsClassifier ,genes))


write.table(as.data.frame(IntegratedAdultsAnnotated),"Adults_Results/Integrated_Adults_Results/IntegratedAdultsNormalized.csv",sep = ",",row.names = T,col.names = T)



IntegratedMetaData=as.data.frame(append(rep("Linda",55),rep("leital",dim(lietalcountData)[1])))

colnames(IntegratedMetaData)="cohort"

form_Integrated_Adults <- ~ (cohort) 


checkVariableEffects(list(as.data.frame(IntegratedAdultsClassifier)),list(as.data.frame(IntegratedMetaData)),form_Integrated_Adults,paste(getwd(),"/Adults_Results/Integrated_Adults_Results/",sep=""),"variableEffectsAdultsIntegratedBefore")


#----------------Box plots for anti and pro apoptotic---------------------------------------------

limit=1:length(sampleUsageAdults)

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BCL2"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BCL2","decision")

p1<-ggplot(df.m, aes(x=decision, y=BCL2)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BCL2A1"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BCL2A1","decision")

p2<-ggplot(df.m, aes(x=decision, y=BCL2A1)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BCL2L1"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BCL2L1","decision")
p3<-ggplot(df.m, aes(x=decision, y=BCL2L1)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BCL2L2"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BCL2L2","decision")
p4<-ggplot(df.m, aes(x=decision, y=BCL2L2)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"MCL1"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("MCL1","decision")
p5<-ggplot(df.m, aes(x=decision, y=MCL1)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

dev.off()
tiff(paste("Adults_Results/Integrated_Adults_Results/Anti-Apoptotic","-",Sys.Date(),".tiff"), units="in", width=10, height=5, res=400)
plot_grid( p1, p2,p3,p4,p5, labels="AUTO",label_size = 10)
dev.off()
#----------------Box plots for anti and pro apoptotic---------------------------------------------
df.m=data.frame(IntegratedAdultsAnnotated[limit,"BAD"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BAD","decision")

p1<-ggplot(df.m, aes(x=decision, y=BAD)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BAK1"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BAK1","decision")

p2<-ggplot(df.m, aes(x=decision, y=BAK1)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BAX"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BAX","decision")
p3<-ggplot(df.m, aes(x=decision, y=BAX)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BBC3"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BBC3","decision")
p4<-ggplot(df.m, aes(x=decision, y=BBC3)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BCL2L11"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BCL2L11","decision")
p5<-ggplot(df.m, aes(x=decision, y=BCL2L11)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

 
df.m=data.frame(IntegratedAdultsAnnotated[kimit,"BID"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BID","decision")
p6<-ggplot(df.m, aes(x=decision, y=BID)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"BMF"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("BMF","decision")
p7<-ggplot(df.m, aes(x=decision, y=BMF)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

df.m=data.frame(IntegratedAdultsAnnotated[limit,"PMAIP1"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("PMAIP1","decision")
p7<-ggplot(df.m, aes(x=decision, y=PMAIP1)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")
dev.off()

df.m=data.frame(IntegratedAdultsAnnotated[limit,"CD6"],decisionIntegratedAdultsDF[limit])
colnames(df.m)<-c("CD6","decision")
p8<-ggplot(df.m, aes(x=decision, y=CD6)) + 
  geom_boxplot()+geom_jitter(shape=16, position=position_jitter(0.2))+labs(x = "")

dev.off()

tiff(paste("Adults_Results/Integrated_Adults_Results/Pro-Apoptotic","-",Sys.Date(),".tiff"), units="in", width=10, height=5, res=400)
plot_grid( p1, p2,p3,p4,p5,p6,p7,p8, labels="AUTO",label_size = 10)
dev.off()


