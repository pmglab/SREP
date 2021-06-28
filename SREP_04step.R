# One:Unicox primary select signature relative with prognosis #
# install.packages("Matrix")
Gene_pair_Signature_fun <- function(path=null,inputfile=null,outputfile=null,p=null){
  library(survival)

  setwd(path)                                                                            #set the directory by yourself
  pFilter=p                                                                              #setup the filtration paremeter, default value 0.05
  rtTrain=read.table(inputfile,header=T,sep="\t",check.names=F,row.names=1)              #read the inputfile
  
  outTab=data.frame()
  sigGenes=c("os","censor")
  for(gene in colnames(rtTrain[,3:ncol(rtTrain)])){
    cox=coxph(Surv(os, censor) ~ rtTrain[,gene], data = rtTrain)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    
    if(coxP<pFilter){
      diff=survdiff(Surv(os, censor) ~rtTrain[,gene],data = rtTrain)
      pValue=1-pchisq(diff$chisq,df=1)
      if(pValue<pFilter){
        sigGenes=c(sigGenes,gene)
        outTab=rbind(outTab,
                     cbind(gene=gene,
                           KM=pValue,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           coxPvalue=coxP) )
      }
    }
  }
  write.table(outTab,file="cggaPair_Train_UniCox.txt",sep="\t",row.names=F,quote=F)     #The outputfile of gene and p-value
  surSigExp=rtTrain[,sigGenes]
  surSigExp=cbind(id=row.names(surSigExp),surSigExp)
  write.table(surSigExp,file="cggaPairtTrain_UniSigExp.txt",sep="\t",row.names=F,quote=F)
  
  # Two: Lasso-cox method #
  
  library(survminer)
  library(glmnet)
  
  rtUni=read.table("cggaPairtTrain_UniSigExp.txt",header=T,sep="\t",row.names=1,check.names=F)      
  rtUni$os[rtUni$os<=0]=0.003
  set.seed(1234)
  
  #construct Lasso model
  x=as.matrix(rtUni[,c(3:ncol(rtUni))])
  y=data.matrix(Surv(rtUni$os,rtUni$censor))
  fit=glmnet(x, y, family = "cox", maxit = 1000)
  pdf("cggaPairTrain_lambda.pdf")
  plot(fit, xvar = "lambda", label = TRUE)
  dev.off()
  
  cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)
  pdf("cggaTrain_cvfit.pdf")
  plot(cvfit)
  abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
  dev.off()
  
  
  #output model fomula
  coef=coef(fit, s = cvfit$lambda.min)
  index=which(coef != 0)
  actCoef=coef[index]
  lassoGene=row.names(coef)[index]
  geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
  write.table(geneCoef,file="cggaTrain_geneCoef.txt",sep="\t",quote=F,row.names=F)
  write.table(lassoGene,file="cggaTrain_lassoGene.txt",sep="\t",quote=F,row.names=F,col.names=F)
  
  # Three: Cox-model #
  library(survival)
  rt <- read.table("cggaTrain_lassoGene.txt",header= F,sep="\t",check.names=F)
  
  rt1 <- read.table("cggaPairtTrain_UniSigExp.txt",header = T,sep = "\t",check.names = F)
  
  colnum<-which(colnames(rt1) %in% rt[,1])
  
  rt2<-cbind(rt1[,c(1:3)],rt1[,colnum])
  
  write.table(rt2, file = "R_multiInput.txt",sep="\t",quote=F,row.names=F)
  
  
  library(survival)
  library(survminer)
  set.seed(1234)
  
  rt3=read.table("R_multiInput.txt",header=T,sep="\t",check.names=F,row.names=1)
  
  multiCox=coxph(Surv(os, censor) ~ ., data = rt3)
  multiCox=step(multiCox,direction = "both")
  multiCoxSum=summary(multiCox)
  multiCoxSum
  
  outTab=data.frame(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outTab=data.frame(id=row.names(outTab),outTab)
  
  library(stringr)
  library(magrittr)
  outTab$id = outTab$id%>%str_remove_all('`')
  write.table(outTab, file = 'multicox.txt', sep = '\t', row.names = F, quote = F)
  write.table(outTab$id, file = 'Multicox_gene.txt', sep = '\t', row.names = F, quote = F)
  
  riskScore=predict(multiCox,type="risk",newdata=rt3)
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("os","censor",coxGene)
  risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
  write.table(cbind(id=rownames(cbind(rt3[,outCol],riskScore,risk)),cbind(rt3[,outCol],riskScore,risk)),
              file=outputfile,
              sep="\t",
              quote=F,
              row.names=F)
}  

Gene_pair_Signature_fun("D:/GBM",inputfile = "cggatrain_event.txt", outputfile = "cggatrain_risk.txt",0.05)