
Classifier <- setRefClass("Classifier",
                          #          fields = list(classifier = "data.frame",MCFSFeatures = "list", Accuracies = "list",classifierMCFSgenes="data.frame",
                          #             path="string",Accuracies_Genetic="list",classifierMCFSgenes="list",Classifier_MCFS="as.data.frame",resultRosettaMCFSGenetic="as.data.frame"
                          # ),
                          fields = list(classifier = "data.frame",flagAccuracy="character",MCFSFeatures = "character",rank="numeric", Accuracies = "numeric",AccuraciesGenetic="numeric", path="character",classifierMCFSgenes="character",Classifier_MCFS="data.frame",annotation="character",resultRosettaMCFSGenetic="list",resultRosettaMCFSJohnson="list",resultRosettaMCFSMerged="list",recalculatedResultRosettaMCFSJohnson="data.frame",recalculatedResultRosettaMCFSGenetic="data.frame",enrichment="list",kegg="gseaResult",ontology="character",numberOfFeatures="numeric",keyType="character",flagEnrichment="logical",underSample="logical"),
                          methods =  list(findAccuracies = function() {
                            classifier <<- classifier
                            MCFSFeatures<<-MCFSFeatures
                            flagAccuracy<<-flagAccuracy
                            underSample<<-underSample
                           #ontology<<--ontology
                            numberOfFeatures<<--numberOfFeatures
                          #  keyType<<--keyType
                           # kegg<<--kegg
                           # enrichment<<--enrichment
                            
                            if(flagAccuracy=="Genetic")
                            {
                              if(underSample==TRUE)
                                AccuraciesGenetic<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy,underSample = TRUE)
                              else
                                AccuraciesGenetic<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy,underSample = FALSE)
                              
                            }else if(flagAccuracy=="Johnson")
                            {
                              Sys.sleep(4)
                              if(underSample==TRUE)
                                Accuracies<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy,TRUE)
                              else
                                Accuracies<<-compareAccuracies(classifier,200,as.character(MCFSFeatures),flagAccuracy,FALSE)
                              
                              
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
                            underSample<<-underSample
                            #   recalculatedResultRosettaMCFS<<-recalculatedResultRosettaMCFS
                            classifier<<-classifier
                            flagAccuracy<<-flagAccuracy
                            resultRosettaMCFSMerged<<-resultRosettaMCFSMerged
                            genes<<-genes
                            if(flagAccuracy=="Johnson")
                            {
                              if(length(which(Accuracies==max(Accuracies))>1))
                                
                              {
                                temp=which(Accuracies==max(Accuracies))
                                maxAccuracy<-temp[length(temp)]*10
                              }else{
                                maxAccuracy<-which(Accuracies==max(Accuracies))*10
                              }
                            }else if(flagAccuracy=="Genetic") {
                              if(length(which(AccuraciesGenetic==max(AccuraciesGenetic))>1))
                                
                              {
                                temp=which(AccuraciesGenetic==max(AccuraciesGenetic))
                                maxAccuracy<-temp[length(temp)]*10
                              }else{
                                maxAccuracy<-which(AccuraciesGenetic==max(AccuraciesGenetic))*10
                              }
                              
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
                            
                            if(underSample==TRUE)
                              resultRosettaMCFSGenetic<<-rosetta(Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5,reducer="Genetic",underSample =T,ruleFiltration=TRUE,ruleFiltrSupport=c(1,3) , discreteParam=3)
                            else
                              resultRosettaMCFSGenetic<<-rosetta(Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5,reducer="Genetic",underSample =F,ruleFiltration=TRUE,ruleFiltrSupport=c(1,3) , discreteParam=3)
                            
                            Sys.sleep(10)
                            # resultRosettaMCFSJohnson<<-NULL
                            if(underSample==TRUE)
                              resultRosettaMCFSJohnson<<-rosetta(Classifier_MCFS,cvNum=5,underSample=T)
                            else
                              resultRosettaMCFSJohnson<<-rosetta(Classifier_MCFS,underSample=F)
                            resultRosettaMCFSMerged<<-list("Johnson"=resultRosettaMCFSJohnson,"Genetic"=resultRosettaMCFSGenetic)
                            
                            recalculatedResultRosettaMCFSJohnson<<-recalculateRules(Classifier_MCFS,resultRosettaMCFSJohnson$main)
                            Sys.sleep(10)
                            # recalculatedResultRosettaMCFSGenetic<<-NULL
                            recalculatedResultRosettaMCFSGenetic<<-recalculateRules(Classifier_MCFS,resultRosettaMCFSGenetic$main)
                            resultRosettaMCFSMerged<<-list("JohnsonResult"=resultRosettaMCFSJohnson,"GeneticResult"=resultRosettaMCFSGenetic,"recalculatedResultJohnson"=recalculatedResultRosettaMCFSJohnson,"recalculatedResultGenetic"=recalculatedResultRosettaMCFSGenetic)
                            plotDistributionCategory(list(filtered_MCFS[[1]]$gene_type),list(c("MCFS")),paste("Distribution_Classes_FS","_",flagAccuracy,sep=""),paste(getwd(),path,sep=""),"")
                            
                          },
                          
                          
                          clusterRulesandWriteoutput = function() {
                            #recalculatedResultRosettaMCFS<<-recalculatedResultRosettaMCFS
                            resultRosettaMCFSGenetic<<-resultRosettaMCFSGenetic
                            recalculatedResultRosettaMCFSGenetic<<-recalculatedResultRosettaMCFSGenetic
                            resultRosettaMCFSJohnson<<-resultRosettaMCFSJohnson
                            recalculatedResultRosettaMCFSJohnson<<-recalculatedResultRosettaMCFSJohnson
                            path<<-path
                            enrichment<<-enrichment
                            flagEnrichment<<-flagEnrichment
                            if(flagEnrichment==FALSE)
                            {
                            tempJohnson=recalculatedResultRosettaMCFSJohnson[(which(!(recalculatedResultRosettaMCFSJohnson$supportLHS==0))),]
                            tempGenetic=recalculatedResultRosettaMCFSGenetic[(which(!(recalculatedResultRosettaMCFSGenetic$supportLHS==0))),]
                            enrichment<<-enrichment
                              
                            clusteredRulesMCFS<<-clusterRules(tempJohnson,rownames(Classifier_MCFS))
                            
                            clusteredRulesMCFSGenetic<<-clusterRules(tempGenetic,rownames(Classifier_MCFS))
                            enrichment<<-computeEnrichment()
                             
                            writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Genetic",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFSGenetic,tempGenetic,resultRosettaMCFSGenetic$main,enrichment,"enrichment",FALSE,FALSE)
                         
                            writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Johnson",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFS, tempJohnson,resultRosettaMCFSJohnson$main,enrichment ,"enrichment",FALSE,FALSE)
                            }
                            else
                            {
                              print(length(enrichment))
                              writeOutput(paste(getwd(),path,sep = ""),Sys.Date(),"AllGenes","Genetic",append(annotation,filtered_MCFS[[2]]),clusteredRulesMCFSGenetic,tempGenetic,resultRosettaMCFSGenetic$main,enrichment,"enrichment",TRUE,FALSE)
                              
                              
                              
                            }
                            
                          }, computeEnrichment = function() {
                            
                            
                            MCFSFeatures<<-MCFSFeatures
                            classifier<<-classifier
                            enrichment<<-enrichment
                            ontology<<-ontology
                            numberOfFeatures<<-numberOfFeatures
                            keyType<<-keyType
                            classifierMCFSgenes<<-classifierMCFSgenes
                            optimalNumber<<-length(classifierMCFSgenes)
                            enrichment<<-enrichment
                            
                            enrichment1=GOenrichment(MCFSFeatures[1:numberOfFeatures],colnames(classifier)[1:length(classifier)-1],keyType,ontology)
                            enrichment2=GOenrichment(MCFSFeatures[1:optimalNumber],colnames(classifier)[1:length(classifier)-1],keyType,ontology)
                            enrichment3=GOenrichment(MCFSFeatures,colnames(classifier)[1:length(classifier)-1],keyType,ontology)
                            enrichment=list(enrichment1=enrichment1,enrichment2=enrichment2,enrichment3=enrichment3)
                           # plotEnrichment(enrichment,"Enrichment  All Genes",paste(paste(getwd(),path,sep = ""),Sys.Date(),"/","enrichment.csv",sep = ""))
                            
                          },
                          computeKEGG = function() {
                            kegg<<-kegg
                            keyType<<-keyType
                            MCFSFeatures<<-MCFSFeatures
                            numberOfFeatures<<-numberOfFeatures
                            rank<<-rank
                            geneList=NULL
                            geneList=1:numberOfFeatures
                            names(geneList)=MCFSFeatures[1:numberOfFeatures]
                            
                            eg = bitr(MCFSFeatures[1:numberOfFeatures], fromType=keyType,toType="ENTREZID", OrgDb="org.Hs.eg.db")
                            # geneList[match(eg$ENSEMBL, names(geneList))]
                            # temp=geneList[match(eg$ENSEMBL, names(geneList))]
                            if(keyType=="ENSMBL")
                            {
                              temp=eg[match(names(geneList),eg$ENSEMBL),]
                              temp=temp[which(!is.na(temp)),]
                              temp=temp[which(!is.na(temp$ENSEMBL)),]
                              temp=unique(temp)
                              
                              geneList=geneList[names(geneList)%in% temp$ENSEMBL]
                              View(geneList)
                            } else if(keyType=="SYMBOL")
                              
                            {
                              temp=eg[match(names(geneList),eg$SYMBOL),]
                              temp=temp[which(!is.na(temp)),]
                              temp=temp[which(!is.na(temp$SYMBOL)),]
                              temp=unique(temp)
                              
                              geneList=geneList[names(geneList)%in% temp$SYMBOL]
                              View(geneList)
                              
                              
                            }
                            
                            names(geneList)=temp$ENTREZID
                            
                            kegg <- gseKEGG(geneList     = sort(geneList, decreasing = TRUE),
                                            organism     = 'hsa',
                                            nPerm        = 1000,
                                            minGSSize    = 1000,
                                            pvalueCutoff = 0.8,
                                            verbose      = FALSE)
                            
                            
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
  features=features[seq(1,numberOfFeatures),c("attribute","RI_norm")]
  return(list(features$attribute,features$RI_norm))
  
}
compareAccuracies=function(Decisiontable,limit,features,flag,underSample)
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
  { if(underSample==TRUE)
      resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter", discreteMethod="EqualFrequency",cvNum = 5,underSample=T,discreteParam=3)
  else
  resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter", discreteMethod="EqualFrequency",underSample=F,discreteParam=3)
  
    } else if(flag =="Genetic")
    {
  #resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=5, discreteParam=3)
 # print("hiii")
 #    print(str(Decisiontable[,append(features[1:i],"decision")]))
 #        #Genetic
 #      View(Decisiontable[,append(features[1:i],"decision")])
 #      print(str(Decisiontable[,append(features[1:i],"decision")]))
      if(underSample==TRUE)
        resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],cvNum = 5,reducer="Genetic",underSample=T,ruleFiltration=TRUE)
    else{
      resultRosetta=rosetta(Decisiontable[,append(features[1:i],"decision")],cvNum = 10,reducer="Genetic",underSample=F,ruleFiltration=TRUE)
      
    }
  
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
