
library(xlsx)
library("RTCGAToolbox")
library("scater")
library("edgeR")
library("caret")
library(xlsx)
library(Boruta)


options(stringsAsFactors = FALSE)
#---------Matched Diagnosis and Relapse---------------------
files<-read.xlsx("TARGET_AML_SampleMatrix_Discovery_20180118.xlsx", sheetName = "Filtered")
colnames(files)<-c("T","R")
files=files[1:dim(files)[1]-1,]
samples=gsub("(TARGET-20-(.)+)-([0-9]{2}A)-[0-9]{2}R","\\1",files[,"T"])
filepath=list.files("mRNA-seq/L3/expression/BCCA/geneQuantification/")
filepath="mRNA-seq/L3/expression/BCCA/geneQuantification/"
TumorMatrix=readQuantification(as.character(files[,1]),paste(samples,"-T",sep=""),filepath)
decision=rep("Diagnosis",dim(TumorMatrix)[1])
RelapseMatrix=readQuantification(as.character(files[,2]),paste(samples,"-R",sep=""),filepath)
decision=append(decision,rep("Relapse1",dim(RelapseMatrix)[1]))
TARGET_GE_Classifier=as.data.frame(rbind(TumorMatrix,RelapseMatrix))
TARGET_GE_Classifier=normalizeGE(TARGET_GE_Classifier,as.factor(decision),TRUE)

#Annotate the whole decision table

TARGET_GE_Classifier_Annotated=annotateDecisionTable(TARGET_GE_Classifier,genes)

write.csv(TARGET_GE_Classifier_Annotated,"ResultsRules/NormalizedMatched.csv",row.names = TRUE,col.names = TRUE,fileEncoding = "UTF-8")


#write.csv.raw(TARGET_GE_Classifier_Annotated, file = "ResultsRules/NormalizedMatched.csv", append = FALSE, sep = ",",  fileEncoding = "")

#logTARGET_GE_CLassifier=TARGET_GE_Classifier+0.00001
#logTARGET_GE_CLassifier=log2(logTARGET_GE_CLassifier)
#-----------------All Data even if they are unmatched------------------------------------------------
files<-read.xlsx("TARGET_AML_SampleMatrix_Discovery_20180118.xlsx",header = FALSE, sheetName = "Diagnosis")
files=as.character(files$X1)
samples=gsub("(TARGET-(20|21)-(.)+)-([0-9]{2}A)-[0-9]{2}R","\\1",files)
#filepath=list.files("mRNA-seq/L3/expression/BCCA/geneQuantification/")
filepath="mRNA-seq/L3/expression/BCCA/geneQuantification/"
TumorMatrixAll=readQuantification(files,paste(samples,"-T",sep=""),filepath)
decisionAll=rep("Diagnosis",dim(TumorMatrixAll)[1])
files<-read.xlsx("TARGET_AML_SampleMatrix_Discovery_20180118.xlsx",header = FALSE, sheetName = "Relapse")
files=as.character(files$X1)
samples=gsub("(TARGET-20-(.)+)-([0-9]{2}A)-[0-9]{2}R","\\1",files)
filepath="mRNA-seq/L3/expression/BCCA/geneQuantification/"
RelapseMatrixAll=readQuantification(files,paste(samples,"-R",sep=""),filepath)
decisionAll=append(decisionAll,rep("Relapse1",dim(RelapseMatrixAll)[1]))
TARGET_GE_ClassifierAll=as.data.frame(rbind(TumorMatrixAll,RelapseMatrixAll))
rownames(TARGET_GE_ClassifierAll)<-append(rownames(TumorMatrixAll),rownames(RelapseMatrixAll))
TARGET_GE_ClassifierAllNormalized=normalizeGE(TARGET_GE_ClassifierAll,as.factor(decisionAll),FALSE,TRUE,TRUE)

write.csv(t(Linda_GE_Classifier),"Linda_Cohort_Results/WholeCohort/Normalized-ENSMBLID.csv")

TARGET_GE_ClassifierAll_Annotated=annotateDecisionTable(TARGET_GE_ClassifierAll,genes)


write.csv(t(TARGET_GE_ClassifierAll),"/Users/saryo614/Desktop/Projects/TARGET_GE_ClassifierAll_Ensmbl.csv")

#-----------------------------------Explore Samples---------------------------------------------
metadata=read.xlsx("harmonized/TARGET_AML_Discovery_ClinicalData_20170525.xlsx",header = TRUE, sheetName = "Clinical Data")
metadata_matched=metadata[which(metadata[,"TARGET.USI"] %in% samples ),]


exploreSamples(list(as.data.frame(TARGET_GE_Classifier)),list("Diagnosis/Relapse"),NULL,getwd(),"ExploreSamples")
plotPCAmeta(list(as.data.frame(TARGET_GE_Classifier)),list(as.data.frame(metadata_matched)),"Gender",getwd(),"")
plotPCAmeta(list(as.data.frame(TARGET_GE_Classifier)),list(as.data.frame(metadata_matched)),"Race",getwd(),"")
plotPCAmeta(list(as.data.frame(TARGET_GE_Classifier)),list(as.data.frame(metadata_matched)),"Ethnicity",getwd(),"")
plotPCAmeta(list(as.data.frame(TARGET_GE_Classifier)),list(as.data.frame(metadata_matched)),"Age.at.Diagnosis.in.Days",getwd(),"cluster")
plotPCAmeta(list(as.data.frame(TARGET_GE_Classifier)),list(as.data.frame(metadata_matched)),"Protocol",getwd(),"")

#---------Finding sources of Variance---------------------------------------------------------
metadata_for_exploratory=rbind(metadata_matched,metadata_matched)
rownames_exploratory=append(paste(metadata_matched[,"TARGET.USI"],rep("-T",dim(TumorMatrix)[1]),sep=""),paste(metadata_matched[,"TARGET.USI"],rep("-R",dim(TumorMatrix)[1]),sep=""))
rownames(metadata_for_exploratory)<-rownames_exploratory
form <- ~ (1|Gender) + (1|Race) + (1|Ethnicity) +  Age.at.Diagnosis.in.Days 

checkVariableEffects(list(as.data.frame(TARGET_GE_Classifier)),list(as.data.frame(metadata_for_exploratory)),form,getwd(),"variableEffects")

form <- ~ (1|Gender) + (1|Race) + (1|Ethnicity) +  Age.at.Diagnosis.in.Days 
temp=normalizeGE(TumorMatrix,rep("Tumor",dim(TumorMatrix)[1]))
checkVariableEffects(list(temp),list(as.data.frame(metadata_matched)),form,getwd(),"variableEffects-Tumor")

form <- ~ (1|Gender) + (1|Race) + (1|Ethnicity) +  Age.at.Diagnosis.in.Days 
temp=normalizeGE(TumorMatrix,rep("Relapse",dim(RelapseMatrix)[1]))
checkVariableEffects(list(temp),list(as.data.frame(metadata_matched)),form,getwd(),"variableEffects-Relapse")



