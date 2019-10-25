#------------------------------Merge TARGET GE classifier and the Linda_GE_Classifier (Matched)-----------------------------------------------------------------------------
#Run the script for creating the "TARGET_GE_Classifier" and "Linda_GE_Classifier" raw values and then get Target_Genes and Linda_Genes from "PediatricAMLClassifier"
Common_Genes=Target_Genes[which(Target_Genes %in% Linda_Genes)]
Integrated_Matched=rbind(TARGET_GE_Classifier[,Common_Genes],Linda_GE_Classifier[,Common_Genes])
decision_Integrated_Matched=append(decision,decision_Linda)
#Create grouping for batches 
batches_integrated=rep("batch1",dim(TARGET_GE_Classifier)[1])
batches_integrated=append(batches_integrated,rep("batch2",dim(Linda_GE_Classifier)[1]))

Integrated_Matched=normalizeGE(Integrated_Matched,as.factor(decision_Integrated_Matched),TRUE)
Integrated_Matched_Corrected=removeBatchEffect(Integrated_Matched,as.factor(batches_integrated),batches_integrated,1)
Integrated_Matched_Corrected=t(Integrated_Matched_Corrected)




#Annotate the whole decision table

Integrated_Matched_Annotated=as.data.frame(annotateDecisionTable(Integrated_Matched_Corrected,genes))


write.table(as.data.frame(Integrated_Matched_Annotated),"Integrated_Cohort_Results/NormalizedMatched.csv",sep = ",",row.names = T,col.names = T)

#An attempt to save as xls files but i have limitations on the number of coloumns

#myfile="Integrated_Cohort_Results/NormalizedMatched.xls" #DO NOT USE SPACES, it will cause errors

#wb = createWorkbook() #create a workbook object

#mysheet = createSheet(wb, "Sheet1") #create a sheet and name it
#addDataFrame(as.data.frame(Integrated_Matched_Annotated), sheet=mysheet, startColumn=1, row.names=FALSE) #add data to Sheet1

#.... add as many sheets/data as you want ....

#saveWorkbook(wb, myfile) #write to a physical Excel file here



#Plot Before Correction
plotPCA(getwd(),Integrated_Matched,batches_integrated,"")
#Plot After Correction
plotPCA(getwd(),Integrated_Matched_Corrected,batches_integrated,"")

exploreSamples(list(as.data.frame(Linda_GE_Classifier)),list("Diagnosis/Relapse"),NULL,paste(getwd(),"/Linda_Cohort_Results/",sep=""),"ExploreSamples")

