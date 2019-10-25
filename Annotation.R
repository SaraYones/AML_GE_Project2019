readGeneDB=function(file)
{
  #genes2 <- fread("gencode.v19.annotation.gtf")
  genes <- fread(file)
  setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
  genes$gene_name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))
  genes$gene_id <- unlist(lapply(genes$attributes, extract_attributes, "gene_id"))
  genes$gene_type<- unlist(lapply(genes$attributes, extract_attributes, "gene_type"))
  genes$geneID<-gsub("((.*)(\\..*))","\\2",genes$gene_id)
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
  
  #reference=Linda_MCFSFeatures[1:50]
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
    {
      annotation=append(annotation,filteredMCFS[[1]]$gene_name[which(filteredMCFS[[1]]$geneID==Features[i])])
    print("I am here")
      }
      else{
      annotation=append(annotation,Features[i])
    }
    
  }
  return(annotation)
  
}

annotateDecisionTable=function(decisiontable,genes)
{
  
  AllGenes=getAnnotatedGenes(genes,colnames(decisiontable))
  
  print(dim(AllGenes))
  
  AllGenes_ENSMBL=append(AllGenes[[1]]$geneID,AllGenes[[2]])
  
  colnames=annotateInOrder(AllGenes,colnames(decisiontable))
  
#  print(colnames)
  
 decisiontable_Annotated=decisiontable
  colnames(decisiontable_Annotated)<-colnames
  
  return(decisiontable_Annotated)
  
}



