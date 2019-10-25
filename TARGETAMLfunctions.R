#---------------Preprocess----------------------------------
readQuantification =function(files,samples,filepath){
  print(files[1])
  tempfile=paste(as.character(files[1]),as.character(".gene.quantification.txt"),sep="")
  print(tempfile)
  tempfile=read.table(file=paste(filepath,as.character(tempfile),sep=""),header = T,sep="\t")
  
  # names(FilteredGeneList)
  quantificationMatrix=matrix(0,nrow=length(samples),ncol=dim(tempfile)[1])
  genes=tempfile[,"gene"]
  colnames(quantificationMatrix)=genes
  rownames(quantificationMatrix)=samples
  for (i in  1:length(samples))
  {
    tempfile=paste(as.character(files[i]),as.character(".gene.quantification.txt"),sep="")
    tempfile=read.table(file=paste(filepath,as.character(tempfile),sep=""),header = T,sep="\t")
    # names(FilteredGeneList)
    raw_counts=tempfile[,"raw_counts"]
    
    quantificationMatrix[as.character(samples[i]),]=raw_counts
  }
  
  return( quantificationMatrix)
}


normalizeGE=function(GeneExpression,group)
{#GeneExpression=TARGET_GE_Classifier
#LAMLGE=t(GeneExpression)
remove_cols=nearZeroVar(GeneExpression)
#group=decision
#cpmd=cpm(t((GeneExpression[,-remove_cols])))
cpmd=t((GeneExpression[,-remove_cols]))
dge <- DGEList(counts=cpmd,group=group)
keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")
dge <- estimateCommonDisp(dge)
dge<-cpm(dge$pseudo.counts,log=TRUE)

#dge <- DGEList(counts=cpmd,group=group)
LAMLGEClassifier=as.data.frame(t(dge))
View(LAMLGEClassifier)
dim(LAMLGEClassifier)
return(LAMLGEClassifier)
}
#--------------Annotation-------------------------------
readGeneDB=function(file)
{
  #genes2 <- fread("gencode.v19.annotation.gtf")
  genes <- fread(file)
  setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
  genes$gene_name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))
  genes$gene_id <- unlist(lapply(genes$attributes, extract_attributes, "gene_id"))
  genes$gene_type<- unlist(lapply(genes$attributes, extract_attributes, "gene_type"))
  # [optional] focus, for example, only on entries of type "gene", 
  # which will drastically reduce the file size
  genes <- genes[type == "gene"]
  return(genes)
  
}

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

getAnnotatedGenes <- function(genes,reference)
{

  #reference=MCFSFeatures[1:20]
  filtered=genes[which(genes$geneID %in% reference),]
  filtered=filtered[order(sapply(filtered$geneID, function(x) which(x == reference))), ]
  View(filtered[order(sapply(filtered$geneID, function(x) which(x == reference))), ])
  
  nonAnnotated=reference[which(!(reference %in% genes$geneID))]
  return(list(filtered,nonAnnotated))
  
}
annotateInOrder=function(filteredMCFS,Features)
{
  
  annotation=NULL
  for(i in 1:length(Features))
  {
    if(Features[i] %in% filteredMCFS[[1]]$geneID)
      annotation=append(annotation,filteredMCFS[[1]]$gene_name[which(filteredMCFS[[1]]$geneID==Features[i])])
    else{
      annotation=append(annotation,Features[i])
    }
    
  }
  return(annotation)
  
}


#--------Plotting Functions--------------------
plotDistributionCategory=function(GE,titles,variable,filepath,discretizemethod)
{  
  
  print(length(GE))
  
  my.plots <- vector(length(GE), mode='list')
 
  
  for(i in 1:length(GE)){
    op <- par(mar = c(16,7,7,4) + 0.1)
    end_point = 0.5 + length(GE[[i]]) #+ length(GE[[i]])-1 #this is the line which does the trick (together with barplot "space = 1" parameter)
    
    plt=barplot(table(GE[[i]]),xlab = "",xaxt = "n",space=1,cex.axis = 0.5 ,cex.names = 0.5,cex.main=0.8,main=titles[i])
    par(op)
    text( plt, par("usr")[3]-0.25, 
         srt = 60, adj= 1, xpd = TRUE,
         labels = names(table(GE[[i]])), cex=0.65)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
    
  }
  
  graphics.off()
 # update_geom_default("point", list(size=1))
 # theme_set(theme_grey(base_size=6))
  sc <- c(0.5,0.75,1)
  
  pdf(width=7*sc[1],height=7*sc[1],pointsize=1,paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}

