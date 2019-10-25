#Adults Preprocessing (Checking for confounders)


#Relapse unmatched
#Filtering Adults according to age If you want to check for the ages of the cases you include

cases_Adults_Relapse1_unmatched<-read.xlsx("Linda_Cohort_Metadata/AdultsData/AdultsMetadata.xlsx", sheetName = "Relapse1")
#The coloumn names are marked with "." so we need to change to -
new_col_names=str_replace_all(rownames(Linda_GE_Classifier2), "\\.","-")
rownames(Linda_GE_Classifier2)<-new_col_names

temp1=cases_Adults_Relapse1_unmatched$AdultDiagnosis[!is.na(cases_Adults_Relapse1_unmatched$AdultDiagnosis)]
temp2=cases_Adults_Relapse1_unmatched$AdultRelapse[!is.na(cases_Adults_Relapse1_unmatched$AdultRelapse)]
cases_Adults_Relapse1_unmatched<-NULL

cases_Adults_Relapse1_unmatched$AdultDiagnosis=temp1
cases_Adults_Relapse1_unmatched$AdultRelapse=temp2

#cases_Adults_Relapse1_unmatched$AdultDiagnosis=str_replace_all(cases_Adults_Relapse1_unmatched$AdultDiagnosis , "_|-", ".")
#cases_Adults_Relapse1_unmatched$AdultRelapse=str_replace_all(cases_Adults_Relapse1_unmatched$AdultRelapse , "_|-", ".")
 

#Filter from the whole cohort only the cases of Diagnosis and Relapse 1 that are adults cases
 filterCases=which(rownames(Linda_GE_Classifier2) %in% append(levels(cases_Adults_Relapse1_unmatched$AdultDiagnosis),levels(cases_Adults_Relapse1_unmatched$AdultRelapse) ))
 Linda_GE_Classifier2=Linda_GE_Classifier2[filterCases,]
 decision_Linda2=decision_Linda2[filterCases]
 
files2_Linda<-read.xlsx("Linda_Cohort_Metadata/Svea/20190123_AML_Metadata.xls", sheetName = "Sheet1")

#Just to check ages of Adult cases
#ages=files2_Linda[which(!is.na(files2_Linda$Age.at.Sampling..Yrs.)),]
#ages=ages[which(ages$Age.at.Sampling..Yrs.!="N/A"),]
#ages$Age.at.Sampling..Yrs.=as.numeric(as.character(ages$Age.at.Sampling..Yrs.))

#remove all the control cases
Linda_GE_Classifier2=cbind.data.frame(Linda_GE_Classifier2,decision_Linda2)

#I dont need this anymore coz i filtered above
#Linda_GE_Classifier2=Linda_GE_Classifier2[-which(grepl("BM.*",rownames (Linda_GE_Classifier2))),]
#if i run the normalize with (run pipeline==FALSE) to remove all the nearzero cols

#because we normalized on the whole dataset so i decided not to remove until everything was normalized
remove_cols=nearZeroVar(Linda_GE_Classifier2)
if(length(remove_cols)!=0)
  Linda_GE_Classifier2=Linda_GE_Classifier2[,-remove_cols]


#check confounders for Linda cohort

metadata_Linda=read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "Sheet1")

metadata_Adults_Linda<-read.xlsx("Linda_Cohort_Metadata/Svea/20190123_AML_Metadata.xls", sheetName = "Sheet1")


#metadata_Linda[,"Sample.ID"] <-str_replace_all(metadata_Linda[,"Sample.ID"] , "_|-", ".")
metadata_Adults_Linda=metadata_Adults_Linda[which(metadata_Adults_Linda[,"Sample.ID"] %in% rownames(Linda_GE_Classifier2) ),]

paste(getwd(),"/Linda_Adults_Cohort_Results/",sep="")

form_Linda_Adults <- ~ (1|Sex) + (1|Stage) + (1|Age.at.onset..Grouped.) +(1|Blast.count.....Grouped.)+(1|FAB.subtype)+(1|EFS..Grouped.)


checkVariableEffects(list(as.data.frame(Linda_GE_Classifier2[,1:dim(Linda_GE_Classifier2)[2]-1])),list(as.data.frame(metadata_Adults_Linda)),form_Linda,paste(getwd(),"/Linda_Adults_Cohort_Results/",sep=""),"variableEffectsAdults")


#check confounders for Lietal cohort


metadata_Lietal=fread("Linda_Cohort_Metadata/AdultsData/LietalMetadata.txt")
#Remove all normal bone marrow samples
metadata_Lietal=metadata_Lietal[1:(dim(metadata_Lietal)[1]-14),]
state=unlist(lapply(metadata_Lietal$Disease_stage,function(x) if(x=="Diagnosis"){return ("DX")}else return ("RX")))
metadata_Lietal$SUBJECT_ID=paste("Patient_",metadata_Lietal$SUBJECT_ID,"-",state,sep = "")

metadata_matched_Lietal=metadata_Lietal[which(metadata_Lietal$SUBJECT_ID %in% rownames(lietalcountData) ),]
rownames(metadata_matched_Lietal)=metadata_matched_Lietal$SUBJECT_ID

#Order according to the same order in LietalCountData
i2=match(rownames(lietalcountData),rownames(metadata_matched_Lietal))
metadata_matched_Lietal=metadata_matched_Lietal[i2,]


paste(getwd(),"/Linda_Adults_Cohort_Results/Lietal_Cohort_Results",sep="")

form_Lietal <- ~ (Age) + (1|sex) + (1|Cytogenetics) +(1|mutGenes)+(1|FAB)+(1|WBC)+(1|TTR)


checkVariableEffects(list(as.data.frame(lietalcountData[,1:dim(lietalcountData)[2]-1])),list(as.data.frame(metadata_matched_Lietal)),form_Lietal,paste(getwd(),"/Linda_Adults_Cohort_Results/Lietal_Cohort_Results",sep=""),"variableEffectsAdults")




