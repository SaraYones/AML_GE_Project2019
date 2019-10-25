#-----------------Functions Classifier-------------------------------------------------
extractFeaturesBoruta=function(boruta.train)
{
  confirmed=rownames(as.data.frame(boruta.train$finalDecision[which(boruta.train$finalDecision=="Confirmed")]))
  tentative=rownames(as.data.frame(boruta.train$finalDecision[which(boruta.train$finalDecision=="Tentative")]))
  BorutaFeatures=append(confirmed,tentative)
  return(BorutaFeatures)
}

FilterFeatures=function(file,numberOfFeatures)
{
  
  features=fread(file)
  features=features[seq(1,numberOfFeatures),"attribute"]
  return(features$attribute)
  
}
compareAccuracies=function(Decisiontable,limit,features)
{
  Accuracies=NULL
  
  for(i in seq(10, limit, by = 10)){
    print(i)
    
    # resultRosettaWithUSMod7=rosetta(DecisiontableBinary,classifier="StandardVoter",discrete=TRUE,underSample = TRUE)
    #For SLE 
    #resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discreteMask=TRUE,discrete = TRUE,ruleFiltration=TRUE,ruleFiltrSupport=c(1,3))
    #FOR AML
 #Johnson
  resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",discreteParam=3)

    #resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5, discreteParam=3)
 #Genetic
   # resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=10,reducer="Genetic",ruleFiltration=TRUE,ruleFiltrSupport=c(1,3) , discreteParam=3)
  # print(resultRosetta)
    Accuracies=append(Accuracies,resultRosetta$quality$accuracyMean)
   # print(Accuracies)
  }
  plot(Accuracies,type='l',xaxt="n",main="Accuracies")
  axis(1,at=seq(1,as.numeric(limit/10), by = 1),labels=as.character(seq(10, limit, by = 10)))
  return(Accuracies)
  
}
#For Gene expression Decision tables
prepareDT=function(DT,values)
{
  temp=NULL;
  temp=DT[,1:dim(DT)[2]-1]
  temp=as.matrix(temp)
  temp[temp=="up"]<-values[3]
  temp[temp=="down"]<-values[1]
  temp[temp=="normal"]<-values[2]
  temp=apply(temp,2,as.numeric)
  temp=as.data.frame(temp)
  decisiontemp=unlist(lapply(DT$decisionSLE,as.character))
  temp=cbind.data.frame(temp,decisiontemp)
  #temp=as.data.frame(apply(temp,2,as.factor))
  names(temp)[names(temp)=="decisiontemp"]="decisionSLE"
  
  #temp$decision=as.factor(temp$decision)
  #temp$decisionSLE=as.integer(temp$decisionSLE)
  rownames(temp)<-rownames(DT)
  return (temp)
}
clusterRules=function(result,objects)
{
  #result=filterResultRosetta13
  #objects=rownames(DescretizedDF13)
  resultMatrix=matrix(0,nrow=dim(result)[1],ncol=length(objects))
  colnames(resultMatrix)<-objects
  for(i in 1:dim(result)[1])
  {
    SUPP_SET_LHS=unlist(as.list(strsplit(as.character(result$supportSetLHS[[i]]), ",")))
    #SUPP_SET_RHS=unlist(as.list(strsplit(as.character(result$SUPP_SET_RHS[[i]]), ",")))
    #SUPP_SET=intersect(SUPP_SET_LHS,SUPP_SET_RHS)
    resultMatrix[i,which(colnames(resultMatrix) %in% SUPP_SET_LHS)]=1
    
  }
  resultMatrix=resultMatrix[,which(colSums(resultMatrix)!=0)]
  #resultMatrix=resultMatrix[,which(colSums(resultMatrix)!=-dim(resultMatrix)[1])]
  return(resultMatrix)
}
