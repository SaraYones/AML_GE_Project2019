discretizeAge=function(values,method)
{
  
  #if using SLE
  # age=as.numeric(gsub("(age: (([0-9]*[.])?[0-9]+))","\\2",values))
  #if using AML data
  age=values
  print(age)
  discretizedAgeFrequency=discretize(age,method = method, breaks = 4,  labels = c("young", "youth","adult","senior"))
  return(discretizedAgeFrequency)
}


removeBatchEffect= function(logGeneExpressionMatrixBatch,pheno,requiredPheno,flag){
  
  #Ask Mateusz how to deal with this because I dont want an exact name for the factor i want to substitute the value of requiredPheno
  #pheno2=as.data.frame(pheno) 
  
  #Creating a model for sva 
  # mod = model.matrix(~cyto_risk_group, data=pheno)
  # mod0 = model.matrix(~1,as.data.frame(pheno))
  #   n.sv = num.sv(t(logGeneExpressionMatrixBatch),mod,method="leek")
  #  svobj = sva(t(logGeneExpressionMatrixBatch),mod,mod0,n.sv=n.sv)
  #   modSv = cbind(mod,svobj$sv)
  #  mod0Sv = cbind(mod0,svobj$sv)
  #fit = lmFit(t(logGeneExpressionMatrixBatch),modSv)
  #Batches known
  if(flag==1)
    #Creating a model for COMBAT with known batches
  {
    modcombat = model.matrix(~1, data=as.data.frame(seq(1,length(pheno))))
    print(modcombat)
    combat_edata = ComBat(dat=as.matrix(t(logGeneExpressionMatrixBatch)), batch=pheno, mod=modcombat)
  }
  
  
  
  return(combat_edata)
}

pheno = pData(bladderEset)

removeColsWithZeroVar=function(GeneExpression)
{
  for(i in 1:length(GeneExpression))
  {
    remove_cols=nearZeroVar(as.data.frame(GeneExpression[i]))
    temp=as.data.frame(GeneExpression[i])
    
    if(length(remove_cols)!=0)
    {
      GeneExpression[[i]]=as.data.frame(temp[,-remove_cols])
    }else{
      GeneExpression[[i]]=as.data.frame(temp)
      
    }
    
  }
  return(GeneExpression)
}

checkVariableEffects=function(GeneExpr,metadata,form,filepath,variable,discretizemethod="")
{ my.plots <- vector(length(GeneExpr), mode='list')
myplots=NULL
for(i in 1:length(GeneExpr)){
  varPart <- fitExtractVarPartModel( t(as.data.frame(GeneExpr[[i]])), form, as.data.frame(metadata[[i]]))
  vp <- sortCols(varPart)
  print(plotVarPart(vp))
  myplots=recordPlot()
  my.plots[[i]]=myplots
}
graphics.off()

pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
for (my.plot in my.plots){
  replayPlot(my.plot)
}
graphics.off()

}


coorelationVariables=function(metadata,decisionSLE,form,filepath,variable,discretizemethod="")
  
{
  
  
  my.plots <- vector(length(metadata), mode='list')
  myplots=NULL
  for(i in 1:length(metadata)){
    print("hello")
    metadata[[i]]=as.data.frame(cbind(metadata[[i]],decisionSLE[[i]]))
    names(metadata[[i]])[names(metadata[[i]])=="decisionSLE[[i]]"]="decisionSLE"
    View(metadata[[i]])
    C = canCorPairs( form, metadata[[i]])
    print( plotCorrMatrix( C ))
    myplots=recordPlot()
    my.plots[[i]]=myplots
  }
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots){
    replayPlot(my.plot)
  }
  graphics.off()
}

exploreSamples=function(GeneExpression,titles,form,filepath,variable,discretizemethod="")
{
  
  par(mar=c(7,5,1,1))
  my.plots <- vector(length(GeneExpression), mode='list')
  myplots=NULL
  for(i in 1:length(GeneExpression)){
    #print("hello")
    
    sample <- rownames(as.data.frame(GeneExpression[[i]]))
    d.f <- data.frame(sample,as.data.frame(GeneExpression[[i]]))
    d.f2 <- melt(d.f, id.vars = "sample")
    
    
    
    boxplot(value~sample,
            data=d.f2,
            main=titles[[i]],
            xlab="samples",
            ylab="gene expr",
            col="orange",
            border="brown",las=2)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
  }
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots){
    replayPlot(my.plot)
  }
  graphics.off()
  
  
}

#For Gene Expression data

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


normalizeGE=function(GeneExpression,group,runpipeline,remove,logflag,secondpassFilter)
{
  #LAMLGE=t(GeneExpression)
  if (secondpassFilter==FALSE)
  {
  if(runpipeline==TRUE)
  {
    remove_cols=nearZeroVar(GeneExpression)
    cpmd=(GeneExpression[,-remove_cols])
    print(dim(cpmd))
    cpmd=t(cpmd)
    print("")
  }
  else{
    print("here i am ")
    cpmd=t((GeneExpression))
  }
  

  #group=decision
  #cpmd=cpm(t((GeneExpression[,-remove_cols])))
  
  
  if(remove==TRUE)
  {
print("I am here")
  dge <- DGEList(counts=cpmd,group=group)
  keep <- rowSums(cpm(dge)>1) >= 5 
  print(keep)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  print("hello")
  dge <- calcNormFactors(dge, method="TMM")
  print("hi")
  dge <- estimateCommonDisp(dge)
  
  }else{
   #cpmd=t(GeneExpression)
    dge <- DGEList(counts=cpmd, group=group)
    dge <- calcNormFactors(dge, method="TMM")
    dge <- estimateCommonDisp(dge)
      #dge<-cpm(dge$pseudo.counts)
    }

    if(logflag==TRUE)
    {
      print("hi")
      dge<-cpm(dge$pseudo.counts,log=TRUE)
    }else{
      dge<-cpm(dge$pseudo.counts) 
      }

  #dge <- DGEList(counts=cpmd,group=group)
  LAMLGEClassifier=as.data.frame(t(dge))
  View(LAMLGEClassifier)
  dim(LAMLGEClassifier)
  return(LAMLGEClassifier)
  }
  else{
    
    remove_cols=nearZeroVar(GeneExpression)
    if(length(remove_cols)!=0)
    {
      cpmd=(GeneExpression[,-remove_cols])
    }else{
      cpmd=GeneExpression
    }
    dge <- DGEList(counts=t(cpmd),group=group)
    keep <- rowSums(cpm(dge$counts)>1)>= 5
    keep=which(keep==TRUE)
    keep=names(keep)
    return(list(remove_genes=remove_cols,keep_genes=keep))
  }
}

orderOnReference=function(matrix,reference)
{
  ordered=matrix[order(sapply(rownames(matrix), function(x) which(x == reference))), ]
  return(ordered)
}
