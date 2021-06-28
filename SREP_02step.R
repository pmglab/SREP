###This part you need to creat a workdiretory, then put the gene-expression profiles in this diretory according to your study. ###
###Matrixfiles have no numberic definition, but you need to konw that the initial matrixfiles must be sampleid as rownames while genenames as columnnames.###
###To sum up, what you want to do is to: 1.create a workdiretory; 2.put your candidategene file and bulk-cell gene-expression proflies in this diretory ###
Relative_expression_genepair_construction_fun <- function(workdir, candidategenefile = null){
	library(limma)
	setwd(workdir)
	v <- list.files("./")
	index <- grep("_matrix.txt", v, ignore.case = F)
	c <- v[index]
	rt <- read.table(candidategenefile,header= F,sep="\t",check.names=F)

	for (i in 1:length(v[index])){
	  rt1 <- read.table(c[i],header = T,sep = "\t",check.names = F)
	  rt1[,c(1)]
	  strname=strsplit(c[i],split = "_")[[1]][1]
	  filename = paste0(strname,'_diffgene.txt')
	  
	  column <- which(colnames(rt1) %in% rt[,1])   #intersection
	  rt2 <- rt1[,c(1,column)]
	  write.table(rt2, file = filename,sep = '\t',quote = F,row.names = F)
	}

	#Matrix Construction-----------------------------------------------------------------------------------------#
	v1 <- list.files("./")
	index <- grep("_diffgene.txt", v1, ignore.case = F)
	c1 <- v1[index]

	for (i in 1:length(v1[index])){
	  rt=read.table(c1[i],header = F,sep = "\t",row.names=1,check.names = F)
	  rt=as.matrix(rt)
	  strname=strsplit(c1[i],split = "_")[[1]][1]
	  filename = paste0(strname,'_trans.txt')
	  rt0=t(rt)
	  write.table(rt0,file = filename,sep = '\t',quote = F,row.names = F)
	}

	#----------------------------------------------------------------------------------------------------#
	v2 <- list.files("./")
	index <- grep("_trans.txt", v2, ignore.case = F)
	c2 <- v2[index]
	matlist <- list()
	genelist <- list()
	for (i in 1:length(v2[index])){
	  if (grepl("TCGA", c2[i], ignore.case = F)){
		rt=read.table(c2[i],header=T,sep="\t",check.names=F)
		rt=as.matrix(rt)
		rownames(rt)=rt[,1]
		exp=rt[,2:ncol(rt)]
		dimnames=list(rownames(exp),colnames(exp))
		data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
		data=avereps(data)
		#Delete the normal samples
		group=sapply(strsplit(colnames(data),"\\-"),"[",4)
		group=sapply(strsplit(group,""),"[",1)
		group=gsub("2","1",group)
		data=data[,group==0]
		data=data[apply(data,1,mad)>0.5,]     #according to article standard,MAD>0.5
	  }else{
		rt=read.table(c2[i],header=T,sep="\t",check.names=F)
		#strname=strsplit(c[i],split = "_")[[1]][1]
		rt=as.matrix(rt)
		rownames(rt)=rt[,1]
		exp=rt[,2:ncol(rt)]
		dimnames=list(rownames(exp),colnames(exp))
		data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
		data=avereps(data)
		data=data[apply(data,1,mad)>0.5,]     #according to article standard,MAD>0.5
	  }
	  matlist[[i]] <- data
	  genelist[[i]] <- row.names(data)
	}

	sameGene=Reduce(intersect,genelist)

	for(i in 1:length(c2)){
	  
	  mat <- matlist[[i]]
	  mat <- mat[sameGene,]
	  mat=rbind(ID=colnames(mat),mat)
	  
	  mat_name = strsplit(c2[i],"_")[[1]][1]
	  mat_name = paste0(mat_name,'_expshare.txt')
	  write.table(mat, file=mat_name, sep="\t", quote=F, col.names=F)
	  
	}

	#Please ignore this section, it doesn't contribute to this study-----------------------------------------------------------------------------------#
	# v3 <- list.files("./")
	# index <- grep("_expshare.txt", v3, ignore.case = F)
	# c3 <- v3[index]
	# for (i in 1:length(v3[index])){
	#   rt=read.table(c3[i],header = F,sep = "\t",row.names=1,check.names = F)
	#   rt=as.matrix(rt)
	#   strname=strsplit(c3[i],split = "_")[[1]][1]
	#   filename = paste0(strname,'_sharetrans.txt')
	#   rt0=t(rt)
	#   write.table(rt0,file = filename,sep = '\t',quote = F,row.names = F)
	# }

	#------------------------------------------------------------------------------------------------#
	v4 <- list.files("./")
	index <- grep("_expshare.txt", v4, ignore.case = F)
	c4 <- v4[index]

	matlist0 <- list()
	genelist0 <- list()
	for (k in 1:length(v4[index])){
	  strname <- strsplit(c4[k], split = "_")[[1]][1]
	  newdata <- paste0(strname,"Pair")
	  newdata <- data.frame()
	  rt = read.table(c4[k],header=T,sep="\t",check.names=F,row.names=1)
	  sampleNum=ncol(rt)
	  for(i in 1:(nrow(rt)-1)){
		for(j in (i+1):nrow(rt)){
		  pair=ifelse(rt[i,]>rt[j,],1,0)
		  pairRatio=sum(pair)/sampleNum
		  if((pairRatio>0.2)&(pairRatio<0.8)){
			rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
			newdata=rbind(newdata,pair)
		  }
		}
	  }
	  matlist0[[k]] <- newdata
	  genelist0[[k]] <- row.names(newdata)
	}

	sameGene2=Reduce(intersect,genelist0)

	for(i in 1:length(c4)){
	  
	  mat <- matlist0[[i]]
	  mat <- mat[sameGene2,]
	  mat=rbind(ID=colnames(mat),mat)
	  
	  mat_name = strsplit(c4[i],split = "_")[[1]][1]
	  mat_name = paste0(mat_name,'_genePair.txt')
	  write.table(mat, file=mat_name, sep="\t", quote=F, col.names=F)
	  
	}
}

Relative_expression_genepair_construction_fun("D:/GBM","RUNNER_gene.txt")