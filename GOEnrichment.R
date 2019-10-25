GOenrichment=function(features,background,flag)
{
  
  if (flag=="ENSEMBL")
  {
    
    ego <- enrichGO(gene          = features,
                    universe      = background,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.01)
    
    bp2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    
  }
  return(bp2)
  
  
  
}

findReleventBP=function(enrichment,features,decisiontable,pval)
{
#  enrichment=enrichmentTARGET
 # features=append(MCFSFeatures[1:20],"decision")
#  decisiontable=TARGET_GE_Classifier
  result=rosetta(decisiontable[,features],classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",reducer="Genetic",ruleFiltration=TRUE, discreteParam=3)
  result=recalculateRules(decisiontable[,features],result$main)
  features=getFeatures(decisiontable[,features],result,1.0)
  #based that we have only 2 decisions can be extended for multi class
  features=unique(append(levels(features[[1]]),levels(features[[2]])))
  
  print(class(features))
  print(features)
  BP=vector(length(features), mode='list')
  print(length(BP))
  print(length(BP))
  print(length(features))
  for(i in 1:length(features))
  {
    print(features[i])
    print("Iam here")
   for(j in 1:dim(enrichment)[1])
  {
    
    temp=strsplit(enrichment[i]$geneID, "/")[[1]]
    if(features[i] %in% temp)
      {BP[i]=enrichment[i]$Description
      print('yes')
     }else{BP[i]="None"}
    
   }
  }
  
   return(as.data.frame(cbind(features,BP)))
}
 
    
#findReleventBP(enrichmentTARGET,append(MCFSFeatures[1:20],"decision"),TARGET_GE_Classifier)

