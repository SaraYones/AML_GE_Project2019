#---------Matched Diagnosis and Relapse---------------------
#Matched
#files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "ChildrenMatched")
#All Cohort
files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "AllSamples")

files2_Linda<-read.xlsx("Linda_Cohort_Metadata/Svea/20190123_AML_Metadata.xls", sheetName = "Sheet1")
files2_Linda<-files2_Linda[,c("Sample.ID","Age.at.Sampling..Yrs.")]
#Change to Matched or Unmatched
#files_Linda<-str_replace_all(files_Linda$Matched, "_|-", ".")
#files_Linda<-str_replace_all(files_Linda$Unmatched, "_|-", ".")
files_Linda<-str_replace_all(files_Linda$All, "_|-", ".")
#files=files[1:dim(files)[1]-1,]
#samples=gsub("(TARGET-20-(.)+)-([0-9]{2}A)-[0-9]{2}R","\\1",files[,"T"])
#filepath=list.files("Linda_Cohort_gene_counts/")
#filepath="mRNA-seq/L3/expression/BCCA/geneQuantification/"
#Linda_GE_Classifier=readQuantification(as.character(files_Linda$Matched),as.character(files_Linda$Matched),filepath)
Linda_GE_Classifier=fread("Linda_Cohort_gene_counts/lanes_combined_merged_gene_counts.txt")

#Data of lietal to compare adults cohort

lietalcountData=readRDS("lietal-countData")
decisionLietal=sapply(gsub("((Patient_AML_[0-9][0-9][0-9]-)(DX|RX))","\\3",colnames(lietalcountData)),function(x) if(x=="DX") {return("Diagnosis")} else return("Relapse1"))
names(decisionLietal)=NULL


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
Linda_GE_Classifier2=fread("Linda_Cohort_gene_counts/NewAdultsGeneCounts/merged_gene_counts.txt")
colnames=Linda_GE_Classifier2$ENSEMBL_ID
Linda_GE_Classifier2=as.data.frame(Linda_GE_Classifier2)

#Linda_GE_Classifier2=Linda_GE_Classifier2[,which(grepl("merged\\.(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))_.*",colnames(Linda_GE_Classifier2)))]

#Use it when you want pediatric cohort only 
Remaining=files_Linda[which(!(files_Linda %in% rownames(Linda_GE_Classifier) ))]

#Use it when we are dealing with all the files
-----------------------------------------------
x=Linda_GE_Classifier2[,2:dim(Linda_GE_Classifier2)[2]]
batches=gsub("(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))_(S[0-9]([0-9])?)_(L[0-9][0-9][0-9])_.*","\\8",colnames(x))
batches=gsub("BM([0-9][0-9][0-9][0-9])_(S[0-9]([0-9])?)_(L[0-9][0-9][0-9])_.*","\\4",batches)
batches=gsub("merged.(L[0-9][0-9][0-9])","\\1",batches)

new_col_names1=gsub("(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))_.*","\\1",colnames(x))
new_col_names1=gsub("(BM([0-9][0-9][0-9][0-9]))_(S[0-9]([0-9])?)_(L[0-9][0-9][0-9])_.*","\\1",new_col_names1)
new_col_names1=gsub("merged\\.(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))","\\1",new_col_names1)
new_col_names=new_col_names1

#--------------------------------------------
#Use it when dealing only with merged files
#-------------------------------------------
new_col_names=gsub("merged\\.(AML([0-9][0-9][0-9])-(([A-Z][0-9]?)|(Ref|P)))_.*","\\1",colnames(Linda_GE_Classifier2))
#---------------------------------------------
new_col_names=str_replace_all(new_col_names, "_|-", ".")
colnames=Linda_GE_Classifier2$ENSEMBL_ID
x=Linda_GE_Classifier2[,2:dim(Linda_GE_Classifier2)[2]]
Linda_GE_Classifier2=x
colnames(Linda_GE_Classifier2)<-new_col_names
#Take from Linda_GE_Classifier2 only the remaining which exist in the Remaining variable now
#Use it only when dealing with merged files
Linda_GE_Classifier2=Linda_GE_Classifier2[,which(new_col_names %in% Remaining )]
#-------------------------------------------------
Linda_GE_Classifier2=t(Linda_GE_Classifier2)
colnames(Linda_GE_Classifier2)<-colnames
#-----------------------------------------------------------
#Linda_GE_Classifier=normalizeGE(Linda_GE_Classifier,as.factor(decision_Linda),TRUE)
x=normalizeGE(Linda_GE_Classifier2,as.factor(decision_Linda2),FALSE,TRUE,TRUE)

#normalizing Lietal data
y=normalizeGE(t(lietalcountData),as.factor(decisionLietal),TRUE,TRUE,TRUE)

#logx=log2(x+0.00001)
#plotPCA(getwd(),logx,batches,"PCA_Adult")
#Linda_GE_Classifier2=removeBatchEffect(logx,as.factor(batches),batches,1)

#If you dont want to remove batch effects first and investigate the cofounderss
Linda_GE_Classifier2temp=Linda_GE_Classifier2
lietalcountDataTemp=lietalcountData

Linda_GE_Classifier2=x
lietalcountData=y

#if you want to get back to normal
Linda_GE_Classifier2=Linda_GE_Classifier2temp
lietalcountData=lietalcountDataTemp
  
Linda_GE_Classifier2=removeBatchEffect(x,as.factor(batches),batches,1)
Linda_GE_Classifier2=removeBatchEffect(Linda_GE_Classifier2,as.factor(batches),batches,1)
Linda_GE_Classifier2=t(Linda_GE_Classifier2)
pcaResult=plotPCA(getwd(),Linda_GE_Classifier2,decision_Linda2,"PCA_Adult")
#Found there is no batch effects from PCA
pcaResult=plotPCA(getwd(),lietalcountData,decisionLietal,"PCA_AdultLietal")


#Annotate the whole decision table
#With ENSMBLID
order_reference=read.xlsx("Linda_AdultCohort_Results/Sample_order_adult_AML.xlsx",sheetName = "Sheet1",header=FALSE)
order_reference=as.character(order_reference[,1])

write.csv(t(orderOnReference(Linda_GE_Classifier2,order_reference)),"Linda_AdultCohort_Results/Normalized-ENSMBLID.csv")

#Save the output file in this format
#Column 1: ENSMBILD
#Column 2: gene name
#Column 3: First AMLXXX
#Column 4-X: As Linda referred: the IDs in order (ascending); right now we have first the unmerged and then the merged
#Column XX: CD34+ BM samples

Linda_GE_Classifier_Annotated=annotateDecisionTable(Linda_GE_Classifier2,genes)

Linda_GE_Classifier_Annotated=orderOnReference(Linda_GE_Classifier_Annotated,order_reference)

Linda_GE_Classifier_Annotated=t(Linda_GE_Classifier_Annotated)

Linda_GE_Classifier_Annotated=cbind(colnames(Linda_GE_Classifier2),Linda_GE_Classifier_Annotated)

Linda_GE_Classifier_Annotated=as.data.frame(Linda_GE_Classifier_Annotated)
#Make rown names the first coloumn
setDT(Linda_GE_Classifier_Annotated, keep.rownames = TRUE)[]
#Set coloumn names again
colnames(Linda_GE_Classifier_Annotated)<-append(c("gene_name","ENSMBL_ID"),order_reference)

write.csv(Linda_GE_Classifier_Annotated,"Linda_AdultCohort_Results/Normalized-Annotated.csv")

write.csv.raw(Linda_GE_Classifier_Annotated, file = "Linda_Cohort_Results/NormalizedUnmatched.csv", append = FALSE, sep = ",", nsep="\t",
              col.names = !is.null(colnames(Linda_GE_Classifier_Annotated)), fileEncoding = "")
#----------------------------------------------------------------------------------------

exploreSamples(list(as.data.frame(Linda_GE_Classifier)),list("Diagnosis/Relapse"),NULL,paste(getwd(),"/Linda_Cohort_Results/",sep=""),"ExploreSamples")
#--------------------------------------------------------------------------------------------------
metadata_Linda=read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "Sheet1")
metadata_Lietal=fread("Linda_Cohort_Metadata/AdultsData/LietalMetadata.txt")
#Remove all normal bone marrow samples
metadata_Lietal=metadata_Lietal[1:(dim(metadata_Lietal)[1]-14),]
state=unlist(lapply(metadata_Lietal$Disease_stage,function(x) if(x=="Diagnosis"){return ("DX")}else return ("RX")))
metadata_Lietal$SUBJECT_ID=paste("Patient_",metadata_Lietal$SUBJECT_ID,"-",state,sep = "")

metadata_Linda[,"Sample.ID"] <-str_replace_all(metadata_Linda[,"Sample.ID"] , "_|-", ".")
metadata_matched_Linda2=metadata_Linda[which(metadata_Linda[,"Sample.ID"] %in% rownames(Linda_GE_Classifier2) ),]
metadata_matched_Lietal=metadata_Lietal[which(metadata_Lietal$SUBJECT_ID %in% rownames(lietalcountData) ),]
rownames(metadata_matched_Linda2)=metadata_matched_Linda2[,"Sample.ID"]
rownames(metadata_matched_Lietal)=metadata_matched_Lietal$SUBJECT_ID
#Order according to the same order in Linda_GE_Classifier2
i1=match(rownames(Linda_GE_Classifier2),rownames(metadata_matched_Linda2))
metadata_matched_Linda2=metadata_matched_Linda2[i1,]
metadata_matched_Linda2$Batches=batches


paste(getwd(),"/Linda_Cohort_Results/",sep="")

form_Linda <- ~ (1|Sex) + (1|Stage) + (1|Age.at.onset..Grouped.) +(1|Blast.count.....Grouped.)+(1|FAB.subtype)+(1|EFS..Grouped.)+(1|Source) +(1|Batches)


checkVariableEffects(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),form_Linda,paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"variableEffects")

plotPCAmeta(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),"Age.at.onset..Grouped.",paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"AfterBatchRemoval")

plotPCAmeta(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),"Sex",paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"AfterBatchRemoval")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),"FAB.subtype",paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"AfterBatchRemoval")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),"EFS..Grouped.",paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"AfterBatchRemoval")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),"Cell.viability..Grouped.",paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"AfterBatchRemoval")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),"Current.sample.resistant.",paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"AfterBatchRemoval")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),"Cell.viability..Grouped.",paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"AfterBatchRemoval")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier2)),list(as.data.frame(metadata_matched_Linda2)),"Source",paste(getwd(),"/Linda_AdultCohort_Results/",sep=""),"")

#--------------------------Unmatched Diagnosis and Relapse--------------------------------------------
files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "ChildrenUnmatched")
#files2_Linda<-read.xlsx("Linda_Cohort_Metadata/Svea/20190123_AML_Metadata.xls", sheetName = "Sheet1")
#files2_Linda<-files2_Linda[,c("Sample.ID","Age.at.Sampling..Yrs.")]
#Repeat the same steps as before (Run the previous script) but decision_Linda will change here
decision_Linda2=vector(length(rownames(Linda_GE_Classifier2)),mode='list')
decision_Linda2[which(grepl("*(.Ref)?(.D)?(.P)?",rownames(Linda_GE_Classifier2)))]="Diagnosis"
decision_Linda2[which(grepl("*.R1(.P)?",rownames(Linda_GE_Classifier2)))]="Relapse1"
decision_Linda2[which(grepl("*.R2(.P)?",rownames(Linda_GE_Classifier2)))]="Relapse2"
decision_Linda2[which(grepl("*.R3(.P)?",rownames(Linda_GE_Classifier2)))]="Relapse3"
decision_Linda2[which(grepl("BM.*",rownames (Linda_GE_Classifier2)))]="Control"
#Relabeled Cases
decision_Linda2[which(grepl("AML056.Ref",rownames (Linda_GE_Classifier2)))]="Relapse1"
decision_Linda2[which(grepl("AML017.R2",rownames (Linda_GE_Classifier2)))]="Relapse1"
rownames(Linda_GE_Classifier2)[which(grepl("AML056.Ref",rownames (Linda_GE_Classifier2)))]="AML056.R1"
rownames(Linda_GE_Classifier2)[which(grepl("AML017.R2",rownames (Linda_GE_Classifier2)))]="AML017.R1"
#---------------
decision_Linda2=unlist(decision_Linda2)



#library(pvclust)

#clustering=pvclust(Linda_GE_Classifier2, method.dist="cor",method.hclust="average", nboot=10)
#result <- pvclust(Countries, method.dist="cor", 
 #                 method.hclust="average", nboot=10)

#Sepecifying the number of clusters
#https://www.r-bloggers.com/pca-and-k-means-clustering-of-delta-aircraft/
#
wss <- (nrow(Linda_GE_Classifier2)-1)*sum(apply(Linda_GE_Classifier2,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(Linda_GE_Classifier2,
                                     centers=i,nstart=25,iter.max = 1000)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")


pc <- prcomp(Linda_GE_Classifier2)

# First for principal components
comp <- data.frame(pc$x[,1:3])
k <- kmeans(comp, 3, nstart=25, iter.max=1000)

palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(comp, col=k$clust, pch=16)

Linda_GE_Classifier2PCA=Linda_GE_Classifier2
Linda_GE_Classifier2PCA$clusters=k$cluster

autoplot(prcomp(Linda_GE_Classifier2PCA), data = Linda_GE_Classifier2, colour = 'clusters', shape = FALSE, label.size = 3)

batches=Linda_GE_Classifier2PCA$clusters
