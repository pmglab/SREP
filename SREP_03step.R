###This part mainly provide code to help split the trainingset and testingset, and combine overall-survival time and status into the final matrix.###
Trainingset_split_combine_oscensor_fun <- function(workdir = null, Inputfile = null, os_censor_file = null, Outputfile = null){

	setwd(workdir) # set the workdir by yourself
	options(stringsAsFactors = F)
	rm(list = ls())

	rt <- read.table(Inputfile,header = F,sep = "\t",row.names=1,check.names = F)    #set the inputfile which you consider it as the trainingset
	rt0=t(rt)
	write.table(rt0,file = "Dataset01Pair_trans.txt",sep = '\t',quote = F,row.names = F)
	
	rt1 <- read.table("Dataset01Pair_trans.txt",sep = '\t',header = T,check.names = F)
	str(rt1)
	set.seed(2)
	ind <- sample(2, nrow(rt1),replace = T, prob=c(0.7,0.3))
	trainset=rt1[ind==1,]
	testset=rt1[ind==2,]
	dim(trainset)
	dim(testset)
	write.table(trainset,file = 'Dataset01train.txt',sep = '\t',quote = F,row.names = F)
	write.table(testset,file = 'Dataset01test.txt',sep = '\t',quote = F,row.names = F)

	rt <- read.table("Dataset01train.txt",header = F,sep = "\t",row.names=1,check.names = F)
	rt0=t(rt)
	write.table(rt0,file = "Dataset01train_trans.txt",sep = '\t',quote = F,row.names = F)

	exp <- read.table("Dataset01train_trans.txt",header = T,sep = "\t")  
	time <- read.table(os_censor_file,header = T,sep = "\t")    #set the file which involved the sampleid, overall_survival and survival_status. 
	colnames(time)[1] <- colnames(exp)[1]

	exp <- as.data.frame(t(exp)) 
	colnames(exp)<- exp[1,] 
	exp<- exp[-1,] 
	result <- cbind(time,exp[match(time$ID,row.names(exp)),])
	result <- na.omit(result)    #exclude those non-avaliable values
	write.table(result,Outputfile,sep = "\t",quote = F,row.names = F)
}
Trainingset_split_combine_oscensor_fun("D:/GBM", "cggaPair.txt", "cgga_os_cens.txt", "cggatrain_event.txt")
