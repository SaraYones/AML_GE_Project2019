#---------Matched Diagnosis and Relapse---------------------
#Matched
#files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "ChildrenMatched")
#All Cohort
files_Linda<-read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "AllSamples")
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
#-----------------------------------------------------------
#Linda_GE_Classifier=normalizeGE(Linda_GE_Classifier,as.factor(decision_Linda),TRUE)
x=normalizeGE(Linda_GE_Classifier,as.factor(decision_LindaAll),TRUE,TRUE,TRUE)
#logx=log2(x+0.00001)
#plotPCA(getwd(),logx,batches,"PCA_Adult")
#Linda_GE_Classifier2=removeBatchEffect(logx,as.factor(batches),batches,1)

#If you dont want to remove batch effects first and investigate the cofounderss
Linda_GE_Classifiertemp=Linda_GE_Classifier

Linda_GE_Classifier=x
#if you want to get back to normal
Linda_GE_Classifier2=Linda_GE_Classifier2temp

Linda_GE_Classifier2=removeBatchEffect(x,as.factor(batches),batches,1)
Linda_GE_Classifier2=t(Linda_GE_Classifier2)
PCAresult=plotPCA(getwd(),Linda_GE_Classifier,as.factor(decision_LindaAll),"PCA_Adult")


#Annotate the whole decision table
#With ENSMBLID
write.csv(t(Linda_GE_Classifier),"Linda_Cohort_Results/WholeCohort/Normalized-ENSMBLID.csv")

Linda_GE_Classifier_Annotated=annotateDecisionTable(Linda_GE_Classifier,genes)


write.csv(t(Linda_GE_Classifier_Annotated),"Linda_Cohort_Results/WholeCohort/Normalized-Annotated.csv")

write.csv.raw(Linda_GE_Classifier_Annotated, file = "Linda_Cohort_Results/NormalizedUnmatched.csv", append = FALSE, sep = ",", nsep="\t",
              col.names = !is.null(colnames(Linda_GE_Classifier_Annotated)), fileEncoding = "")
#----------------------------------------------------------------------------------------

exploreSamples(list(as.data.frame(Linda_GE_Classifier)),list("Diagnosis/Relapse"),NULL,paste(getwd(),"/Linda_Cohort_Results/",sep=""),"ExploreSamples")
#--------------------------------------------------------------------------------------------------
metadata_Linda=read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "Sheet1")


metadata_Linda[,"Sample.ID"] <-str_replace_all(metadata_Linda[,"Sample.ID"] , "_|-", ".")
metadata_matched_Linda=metadata_Linda[which(metadata_Linda[,"Sample.ID"] %in% rownames(Linda_GE_Classifier) ),]
rownames(metadata_matched_Linda)=metadata_matched_Linda[,"Sample.ID"]
#Order according to the same order in Linda_GE_Classifier2
i1=match(rownames(Linda_GE_Classifier),rownames(metadata_matched_Linda))
metadata_matched_Linda=metadata_matched_Linda[i1,]



paste(getwd(),"/Linda_Cohort_Results/",sep="")

form_Linda <- ~ (1|Sex) + (1|Stage) + (1|Age.at.onset..Grouped.) +(1|Blast.count.....Grouped.)+(1|FAB.subtype)+(1|EFS..Grouped.)+(1|Source)


checkVariableEffects(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),form_Linda,paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep=""),"variableEffects"," ")

plotPCAmeta(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),"Age.at.onset..Grouped.",paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep="")," ")

plotPCAmeta(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),"Sex",paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep="")," ")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),"FAB.subtype",paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep="")," ")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),"EFS..Grouped.",paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep="")," ")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),"Cell.viability..Grouped.",paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep="")," ")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),"Current.sample.resistant.",paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep="")," ")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),"Cell.viability..Grouped.",paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep="")," ")
plotPCAmeta(list(as.data.frame(Linda_GE_Classifier)),list(as.data.frame(metadata_matched_Linda)),"Source",paste(getwd(),"/Linda_Cohort_Results/WholeCohort",sep="")," ")

#Creating the decision Variable for the whole cohort

#Repeat the same steps as before (Run the previous script) but decision_Linda will change here
decision_LindaAll=vector(length(rownames(Linda_GE_Classifier)),mode='list')
decision_LindaAll[which(grepl("*(.Ref)?(.D)?(.P)?",rownames(Linda_GE_Classifier)))]="Diagnosis"
decision_LindaAll[which(grepl("*.R1(.P)?",rownames(Linda_GE_Classifier)))]="Relapse1"
decision_LindaAll[which(grepl("*.R2(.P)?",rownames(Linda_GE_Classifier)))]="Relapse2"
decision_LindaAll[which(grepl("*.R3(.P)?",rownames(Linda_GE_Classifier)))]="Relapse3"
decision_LindaAll[which(grepl("BM.*",rownames (Linda_GE_Classifier)))]="Control"

decision_LindaAll=unlist(decision_LindaAll)



library(pvclust)

clustering=pvclust(Linda_GE_Classifier2, method.dist="cor", 
                   method.hclust="average", nboot=10)
result <- pvclust(Countries, method.dist="cor", 
                  method.hclust="average", nboot=10)

#Sepecifying the number of clusters
#https://www.r-bloggers.com/pca-and-k-means-clustering-of-delta-aircraft/
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

