 library(sva)
 library(bladderbatch)
 data(bladderdata)
 library(pamr)
 library(limma)

pheno = as.data.frame(decision_Linda2)

edata = t(Linda_GE_Classifier2)


mod = model.matrix(~as.factor(decision_Linda2), data=pheno)

mod0 = model.matrix(~1,data=pheno)

n.sv = num.sv(edata,mod,method="leek")

svobj = sva(edata,mod,mod0,n.sv=n.sv)

pValues = f.pvalue(edata,mod,mod0)

modSv = cbind(mod,svobj$sv)

mod0Sv = cbind(mod0,svobj$sv)

fit = lmFit(edata,modSv)

contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
> fitContrasts = contrasts.fit(fit,contrast.matrix)


#Clean the data after computing surrogate variables (remove the effect of variables computed by Surragate variables)

corrected=cleanY(t(Linda_GE_Classifier2),mod,svobj$sv)
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}


pcaResult=plotPCA(getwd(),t(corrected),decision_Linda2,"PCA_Adult")


