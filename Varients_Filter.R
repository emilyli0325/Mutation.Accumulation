########################################################################################################################
########################### diploid filter methods adjust ##############################################################
########################################################################################################################
setwd("C:/Users/xli/Desktop/05032018_MA/Result")
cpSNP<- read.delim("merge.SNP.table", sep="\t",header=T)
> dim(cpSNP)
[1] 1071  126

# Extract genotype
genoSNP <- cpSNP[,c(1:6,seq(7,126,3))]


# Extract depth for each locus
dpSNP <- cpSNP[,c(1:6,seq(9,126,3))]
DP.geno <-  cpSNP[,seq(9,126,3)]
DP.geno[is.na(DP.geno)] <- 0.01

# Extract max allele depth for each loci
AdpSNP <- cpSNP[,c(1:6,seq(8,126,3))]

> dim(dpSNP)
[1] 1071   46

# individual consensus
maxA <- function(X){max(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
DPmaxA <- matrix(ncol=length(7:46),nrow=1071)
consIndi <- matrix(ncol=length(7:46),nrow=1071)
for(j in 7:46){
	DPmaxA[,c(j-6)]<-sapply(AdpSNP[,j],maxA,simplify="array")
	consIndi[,c(j-6)]<- DPmaxA[,c(j-6)]/dpSNP[,j]
	}
consIndi.fail.percen <- matrix(ncol=1,nrow=1071)
for(i in 1:1071){
consIndi.fail.percen[i,] <- length(which((consIndi[i,])<0.9))/length(which(is.na(consIndi[i,])== F)) 
}

genoSNP.filt <- genoSNP[,7:46]
genoSNP.filt[consIndi<0.9|DPmaxA<3]<- NA  ### filter consIndi and DP

genoSNP.filt.Addcons <- cbind(cpSNP[,1:6], genoSNP.filt, consIndi.fail.percen)
write.csv(genoSNP.filt.Addcons, "genoSNP.filt.Addcons3.csv")




genoSNP.filt.Addcons <- cbind(cpSNP[,1:6], genoSNP.filt, consIndi.fail.percen)
dim(genoSNP.filt.Addcons)
[1] 5605   81


### generate consREF
GB.REF <- GB.cl.genoNO[,8:47]
> dim(GB.REF)
[1] 5306   40

Mode <- function(x) {x<-x[complete.cases(x)]; ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]}
consensus.strgeno<-apply(GB.REF,1,Mode)
GB.cl.genoNO.REF <- cbind(GB.cl.genoNO, consensus.strgeno)
colnames(GB.cl.genoNO.REF)[7] <- "consFail"


### caculate concordance
concordance <- function(X,Y){a1<-length(which(X==Y)); a2<-length(which(X!=Y)); a3<-cbind(a1,a2); return(a3)}
comp.two.runs <- cbind(GB.cl.genoNO.REF[,1:7], GB.cl.genoNO.REF[,48:49], GB.10.c)
comp.two.runs <- comp.two.runs[comp.two.runs$consFail<0.5,]
comp.two.runs <- comp.two.runs[comp.two.runs$genoNO > 24,] #genotype rate
C.afterrmHet <- matrix(ncol=2,nrow=34)
for(j in seq(10,77,2)){
C.afterrmHet[(j-8)/2,] <- concordance(as.character(comp.two.runs[,j]), as.character(comp.two.runs[,j+1]))
}



### count each sample SNP geno No. and mutant No.
SNP.info.REF <- GB.cl.genoNO.REF[GB.cl.genoNO.REF$consFail<0.5 & GB.cl.genoNO.REF$genoNO > 24,] #genotype rate
MutantNO <- matrix(ncol=1,nrow=40)
for(i in 1:40){
 MutantNO[i,1] <- length(which(SNP.info.REF[,i+7]!= as.character(SNP.info.REF$consensus.strgeno)))
}

Sample.genoNO <- matrix(ncol=1,nrow=40)
for(i in 8:47){
 Sample.genoNO[i-7,1] <- length(which(!is.na(SNP.info.REF[,i])))
}
sum(MutantNO[,1])
sum(C.afterrmHet[,2])
write.csv(SNP.info.REF, "SNP.info.REF.csv")
write.csv(comp.two.runs, "SNP_comp.two.runs.csv")

