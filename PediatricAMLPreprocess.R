#---------Matched Diagnosis and Relapse---------------------
#Matched
#files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "ChildrenMatched")
files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "ChildrenUnmatched")

#All Cohort
#files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "AllSamples")

files2_Linda<-read.xlsx("Linda_Cohort_Metadata/Svea/20190123_AML_Metadata.xls", sheetName = "Sheet1")
files2_Linda<-files2_Linda[,c("Sample.ID","Age.at.Sampling..Yrs.")]

ages=files2_Linda[which(!is.na(files2_Linda$Age.at.Sampling..Yrs.)),]
ages=ages[which(ages$Age.at.Sampling..Yrs.!="N/A"),]
ages$Age.at.Sampling..Yrs.=as.numeric(as.character(ages$Age.at.Sampling..Yrs.))


#Change to Matched or Unmatched
#files_Linda<-str_replace_all(files_Linda$Matched, "_|-", ".")
files_Linda<-str_replace_all(files_Linda$Unmatched, "_|-", ".")
#files_Linda<-str_replace_all(files_Linda$All, "_|-", ".")
#files=files[1:dim(files)[1]-1,]
#samples=gsub("(TARGET-20-(.)+)-([0-9]{2}A)-[0-9]{2}R","\\1",files[,"T"])
#filepath=list.files("Linda_Cohort_gene_counts/")
#filepath="mRNA-seq/L3/expression/BCCA/geneQuantification/"
#Linda_GE_Classifier=readQuantification(as.character(files_Linda$Matched),as.character(files_Linda$Matched),filepath)
Linda_GE_Classifier=fread("Linda_Cohort_gene_counts/lanes_combined_merged_gene_counts.txt")


#rownames(Linda_GE_Classifier)<-as.character(unlist(Linda_GE_Classifier[,1]))
rows=as.character(unlist(Linda_GE_Classifier[,1]))
Linda_GE_Classifier=subset(Linda_GE_Classifier, select=-1)
#rownames(Linda_GE_Classifier)<-rows
Linda_GE_Classifier=t(Linda_GE_Classifier)
names(Linda_GE_Classifier)<-rows
Linda_GE_Classifier=as.data.frame(Linda_GE_Classifier)
Linda_GE_Classifier=Linda_GE_Classifier[which(rownames(Linda_GE_Classifier) %in% files_Linda),]
#Linda_GE_Classifier=t(Linda_GE_Classifier)
colnames(Linda_GE_Classifier)<-rows

#Reading the data and preprocessing it for adults >18
Linda_GE_Classifier2 = setNames(data.frame(t(Linda_GE_Classifier2[,-1])), Linda_GE_Classifier2[,1])
Linda_GE_Classifier2=fread("Linda_Cohort_gene_counts/merged_gene_counts.txt")
colnames=Linda_GE_Classifier2$ENSEMBL_ID
Linda_GE_Classifier2=as.data.frame(Linda_GE_Classifier2)

#Linda_GE_Classifier2=Linda_GE_Classifier2[,which(grepl("merged\\.(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))_.*",colnames(Linda_GE_Classifier2)))]

#Use it when you want pediatric cohort only 
Remaining=files_Linda[which(!(files_Linda %in% rownames(Linda_GE_Classifier) ))]

#Use it when we are dealing with all the files
-----------------------------------------------
batches=gsub("(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))_(S[0-9]([0-9])?)_(L[0-9][0-9][0-9])_.*","\\8",colnames(Linda_GE_Classifier2))
batches=gsub("BM([0-9][0-9][0-9][0-9])_(S[0-9]([0-9])?)_(L[0-9][0-9][0-9])_.*","\\4",batches)
batches=gsub("merged.(L[0-9][0-9][0-9])","\\1",batches)

new_col_names1=gsub("(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))_.*","\\1",colnames(Linda_GE_Classifier2))
new_col_names1=gsub("(BM([0-9][0-9][0-9][0-9]))_(S[0-9]([0-9])?)_(L[0-9][0-9][0-9])_.*","\\1",new_col_names1)
new_col_names1=gsub("merged\\.(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))","\\1",new_col_names1)
new_col_names=new_col_names1
#--------------------------------------------
#Use it when dealing only with merged files
#-------------------------------------------
new_col_names=gsub("merged\\.(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))_.*","\\1",colnames(Linda_GE_Classifier2))
#---------------------------------------------
new_col_names=str_replace_all(new_col_names, "_|-", ".")
colnames(Linda_GE_Classifier2)<-new_col_names
#Take from Linda_GE_Classifier2 only the remaining which exist in the Remaining variable now
#Use it only when dealing with merged files
Linda_GE_Classifier2=Linda_GE_Classifier2[,which(new_col_names %in% Remaining )]
#-------------------------------------------------
colnames=Linda_GE_Classifier2$ENSEMBL.ID
colnames=colnames[2:length(colnames)]
Linda_GE_Classifier2=Linda_GE_Classifier2[,2:dim(Linda_GE_Classifier2)[2]]
x=Linda_GE_Classifier2[,2:dim(Linda_GE_Classifier2)[2]]
Linda_GE_Classifier2=x

#Run these two either you are working with matched or unmatched
Linda_GE_Classifier2=t(Linda_GE_Classifier2)
colnames(Linda_GE_Classifier2)<-colnames
#-----------------------------------------------------------
Linda_GE_Classifier=as.data.frame(rbind(Linda_GE_Classifier,Linda_GE_Classifier2))



decision_Linda<-rep(c("Diagnosis","Relapse1"),ceiling(length(rownames(Linda_GE_Classifier))/2))
#Linda_GE_Classifier=normalizeGE(Linda_GE_Classifier,as.factor(decision_Linda))

Linda_GE_Classifier=normalizeGE(Linda_GE_Classifier,as.factor(decision_Linda),TRUE,TRUE,TRUE)


#Annotate the whole decision table

Linda_GE_Classifier_Annotated=annotateDecisionTable(Linda_GE_Classifier,genes)

write.csv(t(Linda_GE_Classifier_Annotated),"Linda_Cohort_Results/NormalizedUnmatched.csv")

write.csv.raw(Linda_GE_Classifier_Annotated, file = "Linda_Cohort_Results/NormalizedUnmatched.csv", append = FALSE, sep = ",", nsep="\t",
              col.names = !is.null(colnames(Linda_GE_Classifier_Annotated)), fileEncoding = "")
#----------------------------------------------------------------------------------------

exploreSamples(list(as.data.frame(Linda_GE_Classifier)),list("Diagnosis/Relapse"),NULL,paste(getwd(),"/Linda_Cohort_Results/",sep=""),"ExploreSamples")
#--------------------------------------------------------------------------------------------------
metadata_Linda=read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "Sheet1")


metadata_Linda[,"Sample.ID"] <-str_replace_all(metadata_Linda[,"Sample.ID"] , "_|-", ".")
metadata_matched_Linda=metadata_Linda[which(metadata_Linda[,"Sample.ID"] %in% rownames(Linda_GE_Classifier) ),]
#In case we want to run it on Normalized data for the whole cohort
metadata_matched_Linda=metadata_Linda[which(metadata_Linda[,"Sample.ID"] %in% rownames(Linda_GE_Classifier[index,]) ),]


paste(getwd(),"/Linda_Cohort_Results/",sep="")

form_Linda <- ~ (1|Sex) + (1|Stage) + (1|Age.at.onset..Grouped.) +(1|Blast.count.....Grouped.)+(1|FAB.subtype)+(1|EFS..Grouped.)


checkVariableEffects(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),form_Linda,paste(getwd(),"/Linda_Cohort_Results/",sep=""),"variableEffects")
#In case we want to run it on Normalized data for the whole cohort
checkVariableEffects(list(as.data.frame(Linda_GE_Classifier[index,])),list(as.data.frame(metadata_matched_Linda)),form_Linda,paste(getwd(),"/Linda_Cohort_Results/",sep=""),"variableEffects")


#--------------------------Unmatched Diagnosis and Relapse--------------------------------------------
files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "ChildrenUnmatched")
#files2_Linda<-read.xlsx("Linda_Cohort_Metadata/Svea/20190123_AML_Metadata.xls", sheetName = "Sheet1")
#files2_Linda<-files2_Linda[,c("Sample.ID","Age.at.Sampling..Yrs.")]
#Repeat the same steps as before (Run the previous script) but decision_Linda will change here
decision_Linda=vector(length(rownames(Linda_GE_Classifier)),mode='list')
decision_Linda[which(grepl("*.D(.P)?",rownames(Linda_GE_Classifier)))]="Diagnosis"
decision_Linda[which(grepl("*.R1(.P)?",rownames(Linda_GE_Classifier)))]="Relapse1"
decision_Linda=unlist(decision_Linda)


#-------------------------Explore samples  based on the normalization of the whole cohort---------------------
#This script is based on running the NormalizationPediatric a priori
index=which(grepl("Diagnosis|Relapse1",decision_LindaAll))

exploreSamples(list(as.data.frame(Linda_GE_Classifier[index,])),list("Diagnosis/Relapse"),NULL,paste(getwd(),"/Linda_Cohort_Results/WholeCohort/",sep=""),"ExploreSamples")
