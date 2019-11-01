Classifier <- setRefClass("Classifier",
                          #          fields = list(classifier = "data.frame",MCFSFeatures = "list", Accuracies = "list",classifierMCFSgenes="data.frame",
                          #             path="string",Accuracies_Genetic="list",classifierMCFSgenes="list",Classifier_MCFS="as.data.frame",resultRosettaMCFSGenetic="as.data.frame"
                          # ),
                          fields = list(classifier = "data.frame",flagAccuracy="character",MCFSFeatures = "character", Accuracies = "numeric",AccuraciesGenetic="numeric", path="character",classifierMCFSgenes="character",Classifier_MCFS="data.frame",annotation="character",resultRosettaMCFSGenetic="list",resultRosettaMCFSJohnson="list",resultRosettaMCFSMerged="list",recalculatedResultRosettaMCFSJohnson="data.frame",recalculatedResultRosettaMCFSGenetic="data.frame"),
                          methods =  list(findAccuracies = function() {
                            classifier <<- classifier
                            MCFSFeatures<<-MCFSFeatures
                            flagAccuracy<<-flagAccuracy
                            
                            if(flagAccuracy=="Genetic")
                            {
                              
                              AccuraciesGenetic<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy)
                            }else if(flagAccuracy=="Johnson")
                            {
                              Sys.sleep(4)
                              Accuracies<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy)
                              
                              
                              
                            }
                            
                            
                            
                          },
                          createAnnotatedClassifier = function() {
                            
                            Accuracies<<-Accuracies
                            AccuraciesGenetic<<-AccuraciesGenetic
                            classifierMCFSgenes<<-classifierMCFSgenes
                            Classifier_MCFS<<-Classifier_MCFS
                            resultRosettaMCFSJohnson<<- resultRosettaMCFSJohnson
                            resultRosettaMCFSGenetic<<- resultRosettaMCFSGenetic
                            recalculatedResultRosettaMCFSJohnson<<-recalculatedResultRosettaMCFSJohnson
                            recalculatedResultRosettaMCFSGenetic<<-recalculatedResultRosettaMCFSGenetic
                            #   recalculatedResultRosettaMCFS<<-recalculatedResultRosettaMCFS
                            classifier<<-classifier
                            flagAccuracy<<-flagAccuracy
                            resultRosettaMCFSMerged<<-resultRosettaMCFSMerged
                            genes<<-genes
                            if(flagAccuracy=="Johnson")
                            {
                              maxAccuracy<-which(Accuracies==max(Accuracies))*10
                            }else if(flagAccuracy=="Genetic") {
                              
                              maxAccuracy<-which(AccuraciesGenetic==max(AccuraciesGenetic))*10
                              
                            }
                            
                            filtered_MCFS<<-getAnnotatedGenes(genes,MCFSFeatures[1:maxAccuracy])
                            #All Genes MCFS
                            
                            classifierMCFSgenes<<-append(filtered_MCFS[[1]]$geneID,filtered_MCFS[[2]])
                            
                            
                            #Ordering according to MCFS Features
                            classifierMCFSgenes<<-classifierMCFSgenes[!is.na(match(MCFSFeatures, classifierMCFSgenes))]
                            
                            
                            annotation<<-annotateInOrder(filtered_MCFS,classifierMCFSgenes)
                            
                            annotation<<-make.unique(annotation)
                            
                            
                            Classifier_MCFS<<-classifier[,append(classifierMCFSgenes,"decision")]
                            
                            
                            colnames(Classifier_MCFS)<<-append(annotation,"decision")
                            
                            
                            resultRosettaMCFSGenetic<<-rosetta(Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=10,reducer="Genetic",ruleFiltration=TRUE,ruleFiltrSupport=c(1,3) , discreteParam=3)
                            Sys.sleep(10)
                            # resultRosettaMCFSJohnson<<-NULL
                            resultRosettaMCFSJohnson<<-rosetta(Classifier_MCFS)
                            resultRosettaMCFSMerged<<-list("Johnson"=resultRosettaMCFSJohnson,"Genetic"=resultRosettaMCFSGenetic)
                            
                            recalculatedResultRosettaMCFSJohnson<<-recalculateRules(Classifier_MCFS,resultRosettaMCFSJohnson$main)
                            Sys.sleep(10)
                            # recalculatedResultRosettaMCFSGenetic<<-NULL
                            recalculatedResultRosettaMCFSGenetic<<-recalculateRules(Classifier_MCFS,resultRosettaMCFSGenetic$main)
                            resultRosettaMCFSMerged<<-list("JohnsonResult"=resultRosettaMCFSJohnson,"GeneticResult"=resultRosettaMCFSGenetic,"recalculatedResultJohnson"=recalculatedResultRosettaMCFSJohnson,"recalculatedResultGenetic"=recalculatedResultRosettaMCFSGenetic)
                            plotDistributionCategory(list(filtered_MCFS[[1]]$gene_type),list(c("MCFS")),"Distribution_Classes_FS",paste(getwd(),"/Adults_Results/Linda_Adults_Cohort_Results",sep=""),"")
                            
                          },
                          
                          
                          clusterRulesandWriteoutput = function() {
                            #recalculatedResultRosettaMCFS<<-recalculatedResultRosettaMCFS
                            resultRosettaMCFSGenetic<<-resultRosettaMCFSGenetic
                            recalculatedResultRosettaMCFSGenetic<<-recalculatedResultRosettaMCFSGenetic
                            resultRosettaMCFSJohnson<<-resultRosettaMCFSJohnson
                            recalculatedResultRosettaMCFSJohnson<<-recalculatedResultRosettaMCFSJohnson
                            path<<-path
                            
                            clusteredRulesMCFS<<-clusterRules(recalculatedResultRosettaMCFSJohnson,rownames(Classifier_MCFS))
                            
                            clusteredRulesMCFSGenetic<<-clusterRules(recalculatedResultRosettaMCFSGenetic,rownames(Classifier_MCFS))
                            writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Genetic",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFSGenetic,recalculatedResultRosettaMCFSGenetic,resultRosettaMCFSGenetic$main,enrichment = NULL,"enrichment",FALSE)
                            
                            writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Johnson",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFS, recalculatedResultRosettaMCFSJohnson,resultRosettaMCFSJohnson$main,enrichment = NULL,"enrichment",FALSE)
                            
                            
                          }
                          
                          
                          )
                          
)








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
compareAccuracies=function(Decisiontable,limit,features,flag)
{
  Accuracies=NULL
  resultRosetta=NULL
  
  for(i in seq(10, limit, by = 10)){
    print(i)
    
    # resultRosettaWithUSMod7=rosetta(DecisiontableBinary,classifier="StandardVoter",discrete=TRUE,underSample = TRUE)
    #For SLE 
    #resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discreteMask=TRUE,discrete = TRUE,ruleFiltration=TRUE,ruleFiltrSupport=c(1,3))
    #FOR AML
 #Johnson
   # View(Decisiontable)
    #print(class(Decisiontable))
    #print(limit)
    #print(features)
    if(flag =="Johnson")
  { resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter", discreteMethod="EqualFrequency",discreteParam=3)
    } else if(flag =="Genetic")
    {
  #resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5, discreteParam=3)
 # print("hiii")
 #    print(str(Decisiontable[,append(features[1:i],"decision")]))
 #        #Genetic
 #      View(Decisiontable[,append(features[1:i],"decision")])
 #      print(str(Decisiontable[,append(features[1:i],"decision")]))
    resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],cvNum = 10,reducer="Genetic",ruleFiltration=TRUE)
    
  
  }# print(resultRosetta)
    Accuracies=append(Accuracies,resultRosetta$quality$accuracyMean)
    resultRosetta=NULL
    Sys.sleep(4)
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
