########################################################################################################################
####################################################################################################################
########################################################################################################################

setwd("C:/Users/xli/Desktop/table")

### caculate concordance
concordance <- function(X,Y){a1<-length(which(X==Y)); a2<-length(which(X!=Y)); a3<-cbind(a1,a2); return(a3)}

C.afterrmHet <- matrix(ncol=2,nrow=42)
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT[,j]), as.character(GT[,j+40]))
}

### FreeBayes
Bayes.basic <- read.delim("MA.freeBayes.basicfilter.table", sep="\t",header=T)
dim(Bayes.basic)
[1] 13129   249
GT <- Bayes.basic[,seq(7,249,3)]
GT[GT== "."]<-NA

Bayes.hard <- read.delim("MA.freeBayes.hardfilter2.table", sep="\t",header=T)
GT <- Bayes.hard[,seq(7,249,3)]
GT[GT== "."]<-NA


### GATK INDEL
GATK.basicfilter <- read.delim("MA.GATK.basicfilter.INDEL.table", sep="\t",header=T)
GT <- GATK.basicfilter[,seq(7,249,3)]
GT[GT== "."]<-NA

GATK.VQSR <- read.delim("MA.GATK.basicfilter.INDEL2.VQSR.table", sep="\t",header=T)
GT <- GATK.VQSR[,seq(7,249,3)]
GT[GT== "."]<-NA

GATK.hard <- read.delim("MA.GATK.hard.filter2.INDEL.table", sep="\t",header=T)
GT <- GATK.hard[,seq(7,249,3)]
GT[GT== "."]<-NA


### GATK SNP
GATK.basic <- read.delim("MA.GATK.basicfilter.SNP.table", sep="\t",header=T)
GT <- GATK.basic[,seq(7,249,3)]
GT[GT== "."]<-NA

GATK.hard <- read.delim("MA.GATK.hard.filter2.SNP.table", sep="\t",header=T)
GT <- GATK.hard[,seq(7,249,3)]
GT[GT== "."]<-NA

### HipSTR

HipSTR <- read.delim("MA.HipSTR.basicfilter.INDEL.table", sep="\t",header=T)
dim(HipSTR)
[1] 19504   168
GT <- HipSTR[,seq(7,168,2)]
GT[GT== "."]<-NA

########################################################################################################################
####################################################################################################################
########################################################################################################################
# cons filter

### FreeBayes
Bayes.basic <- read.delim("MA.freeBayes.basicfilter.table", sep="\t",header=T)
dim(Bayes.basic)
[1] 13129   249

GT <- Bayes.basic[,seq(7,249,3)]
GT[GT== "."]<-NA

AD <- Bayes.basic[,seq(8,249,3)]
> dim(AD)
[1] 13129    81

# individual consensus
maxA <- function(X){max(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
sumA <- function(X){sum(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
DPmaxA <- matrix(ncol=81,nrow=13129)
DPsumA <- matrix(ncol=81,nrow=13129)


consIndi <- matrix(ncol=81,nrow=13129)
for(j in 1:81){
	DPmaxA[,j]<-sapply(AD[,j],maxA,simplify="array")
	DPsumA[,j]<-sapply(AD[,j],sumA,simplify="array")
	consIndi[,j]<- DPmaxA[,j]/DPsumA[,j]
	}
	
	
consIndi.fail.percen <- matrix(ncol=1,nrow=13129)
for(i in 1:13129){
consIndi.fail.percen[i,] <- length(which((consIndi[i,])<0.8))/length(which(is.na(consIndi[i,])== F)) 
}


GT[consIndi<0.8]<-NA
GT[DPmaxA<3]<-NA


GT.filter <- GT[consIndi.fail.percen<0.5,]
> dim(GT.filter)
[1] 12165    81

C.afterrmHet[1,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,20]))
C.afterrmHet[2,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,60]))
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT.filter[,j]), as.character(GT.filter[,j+40]))
}

Bayes.basic.filter <- cbind(Bayes.basic[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)
write.csv(Bayes.basic.filter, "Bayes.basic.filter.csv")

### GATK INDEL
GATK.Indel <- read.delim("MA.GATK.basicfilter.INDEL.table", sep="\t",header=T)

GT <- GATK.Indel[,seq(7,249,3)]
GT[GT== "."]<-NA

AD <- GATK.Indel[,seq(8,249,3)]
> dim(AD)
[1] 7392   81

# individual consensus
maxA <- function(X){max(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
sumA <- function(X){sum(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
DPmaxA <- matrix(ncol=81,nrow=7392)
DPsumA <- matrix(ncol=81,nrow=7392)


consIndi <- matrix(ncol=81,nrow=7392)
for(j in 1:81){
	DPmaxA[,j]<-sapply(AD[,j],maxA,simplify="array")
	DPsumA[,j]<-sapply(AD[,j],sumA,simplify="array")
	consIndi[,j]<- DPmaxA[,j]/DPsumA[,j]
	}
	
	
consIndi.fail.percen <- matrix(ncol=1,nrow=7392)
for(i in 1:7392){
consIndi.fail.percen[i,] <- length(which((consIndi[i,])<0.8))/length(which(is.na(consIndi[i,])== F)) 
}


GT[consIndi<0.8]<-NA
GT[DPmaxA<3]<-NA


GT.filter <- GT[consIndi.fail.percen<0.5,]
> dim(GT.filter)
[1] 6953   81

C.afterrmHet[1,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,20]))
C.afterrmHet[2,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,60]))
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT.filter[,j]), as.character(GT.filter[,j+40]))
}

GATK.Indel.filter <- cbind(GATK.Indel[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)

### GATK SNP

GATK.SNP <- read.delim("MA.GATK.basicfilter.SNP.table", sep="\t",header=T)

GT <- GATK.SNP[,seq(7,249,3)]
GT[GT== "."]<-NA

AD <- GATK.SNP[,seq(8,249,3)]
> dim(AD)
[1] 6131   81

# individual consensus
maxA <- function(X){max(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
sumA <- function(X){sum(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
DPmaxA <- matrix(ncol=81,nrow=6131)
DPsumA <- matrix(ncol=81,nrow=6131)


consIndi <- matrix(ncol=81,nrow=6131)
for(j in 1:81){
	DPmaxA[,j]<-sapply(AD[,j],maxA,simplify="array")
	DPsumA[,j]<-sapply(AD[,j],sumA,simplify="array")
	consIndi[,j]<- DPmaxA[,j]/DPsumA[,j]
	}
	
	
consIndi.fail.percen <- matrix(ncol=1,nrow=6131)
for(i in 1:6131){
consIndi.fail.percen[i,] <- length(which((consIndi[i,])<0.8))/length(which(is.na(consIndi[i,])== F)) 
}


GT[consIndi<0.8]<-NA
GT[DPmaxA<3]<-NA


GT.filter <- GT[consIndi.fail.percen<0.5,]
> dim(GT.filter)
[1] 6953   81

C.afterrmHet[1,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,20]))
C.afterrmHet[2,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,60]))
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT.filter[,j]), as.character(GT.filter[,j+40]))
}

GATK.SNP.filter <- cbind(GATK.SNP[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)


### HipSTR
HipSTR <- read.delim("MA.HipSTR.basicfilter.INDEL.table", sep="\t",header=T)
dim(HipSTR)
[1] 19504   168
GT <- HipSTR[,seq(7,168,2)]
GT[GT== "."]<-NA

AD <- HipSTR[,seq(8,168,2)]
> dim(AD)
[1] 19504    81


# individual consensus
maxA <- function(X){max(as.numeric(unlist(strsplit(strsplit(as.character(X),";")[[1]],"\\|"))[!c(TRUE,FALSE)]))}
sumA <- function(X){sum(as.numeric(unlist(strsplit(strsplit(as.character(X),";")[[1]],"\\|"))[!c(TRUE,FALSE)]))}

DPmaxA <- matrix(ncol=81,nrow=19504)
DPsumA <- matrix(ncol=81,nrow=19504)


consIndi <- matrix(ncol=81,nrow=19504)
for(j in 1:81){
	DPmaxA[,j]<-sapply(AD[,j],maxA,simplify="array")
	DPsumA[,j]<-sapply(AD[,j],sumA,simplify="array")
	consIndi[,j]<- DPmaxA[,j]/DPsumA[,j]
	}
	
	
consIndi.fail.percen <- matrix(ncol=1,nrow=19504)
for(i in 1:19504){
consIndi.fail.percen[i,] <- length(which((consIndi[i,])<0.8))/length(which(is.na(consIndi[i,])== F)) 
}


GT[consIndi<0.8]<-NA
GT[DPmaxA<3]<-NA


GT.filter <- GT[consIndi.fail.percen<0.5,]
> dim(GT.filter)
[1] 15411    81

C.afterrmHet[1,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,20]))
C.afterrmHet[2,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,60]))
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT.filter[,j]), as.character(GT.filter[,j+40]))
}

HipSTR.filter <- cbind(HipSTR[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)

write.csv(HipSTR.filter, "HipSTR.filter.csv")



########################################################################################################################
####################################################################################################################
########################################################################################################################
# merge two runs
# FreeBayes cons merge 
#Bayes.basic.filter <- cbind(Bayes.basic[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)
dim(Bayes.basic.filter)
[1] 12165    88

collapse_SNP<-function(X,Y){
X<-as.character(X); Y<-as.character(Y);
Z<-NULL; 
Z <- ifelse(is.na(X), Y, X)
Z[which(X!=Y & is.na(X)==F & is.na(Y)==F)]<-NA; 
return(Z)
}

GB.cl<-NULL

for(i in c(9:48)){
    GB.cl<-cbind(GB.cl,collapse_SNP(Bayes.basic.filter[,i],Bayes.basic.filter[,c(i+40)]))
}
> dim(GB.cl)
[1] 12165    40

# caculate genotype number
genoNO <- matrix(ncol=1,nrow=12165)
for(i in 1:12165){
genoNO[i,1] <- length(which(!is.na(GB.cl[i,])))
}

# generate consensus REF
Mode <- function(x) {x<-x[complete.cases(x)]; ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]}
cons.REF <-apply(GB.cl,1,Mode)


# caculate mutant number
MutantNO <- matrix(ncol=1,nrow=12165)
for(i in 1:12165){
 MutantNO[i,1] <- length(which(GB.cl[i,]!= as.character(cons.REF[i])))
}

Bayes.cons.final <- cbind(Bayes.basic.filter[,1:7],cons.REF,genoNO,MutantNO, GB.cl)
Bayes.cons.final2 <- subset(Bayes.cons.final, genoNO>3 & MutantNO>0)
write.csv(Bayes.cons.final2, "Bayes.cons.final2.csv")

# GATK Indel cons merge
# GATK.Indel.filter <- cbind(GATK.Indel[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)
dim(GATK.Indel.filter)
[1] 6953   88
GB.cl<-NULL
for(i in c(9:48)){
    GB.cl<-cbind(GB.cl,collapse_SNP(GATK.Indel.filter[,i],GATK.Indel.filter[,c(i+40)]))
}
# caculate genotype number
genoNO <- matrix(ncol=1,nrow=6953)
for(i in 1:6953){
genoNO[i,1] <- length(which(!is.na(GB.cl[i,])))
}

# generate consensus REF
Mode <- function(x) {x<-x[complete.cases(x)]; ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]}
cons.REF <-apply(GB.cl,1,Mode)


# caculate mutant number
MutantNO <- matrix(ncol=1,nrow=6953)
for(i in 1:6953){
 MutantNO[i,1] <- length(which(GB.cl[i,]!= as.character(cons.REF[i])))
}

GATK.Indel.final <- cbind(GATK.Indel.filter[,1:7],cons.REF,genoNO,MutantNO, GB.cl)
GATK.Indel.final2 <- subset(GATK.Indel.final, genoNO>3 & MutantNO>0)
write.csv(GATK.Indel.final2, "GATK.Indel.final2.csv")


# GATK SNP cons merge
# GATK.SNP.filter <- cbind(GATK.SNP[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)
dim(GATK.SNP.filter)
[1] 4479   88
GB.cl<-NULL
for(i in c(9:48)){
    GB.cl<-cbind(GB.cl,collapse_SNP(GATK.SNP.filter[,i],GATK.SNP.filter[,c(i+40)]))
}
# caculate genotype number
genoNO <- matrix(ncol=1,nrow=4479)
for(i in 1:4479){
genoNO[i,1] <- length(which(!is.na(GB.cl[i,])))
}

# generate consensus REF
Mode <- function(x) {x<-x[complete.cases(x)]; ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]}
cons.REF <-apply(GB.cl,1,Mode)


# caculate mutant number
MutantNO <- matrix(ncol=1,nrow=4479)
for(i in 1:4479){
 MutantNO[i,1] <- length(which(GB.cl[i,]!= as.character(cons.REF[i])))
}

GATK.SNP.final <- cbind(GATK.SNP.filter[,1:7],cons.REF,genoNO,MutantNO, GB.cl)
GATK.SNP.final2 <- subset(GATK.SNP.final, genoNO>3 & MutantNO>0)
write.csv(GATK.SNP.final2, "GATK.SNP.final2.csv")

# HipSTR cons merge
# HipSTR.filter <- cbind(HipSTR[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)

dim(HipSTR.filter)
[1] 15411    88
GB.cl<-NULL
for(i in c(9:48)){
    GB.cl<-cbind(GB.cl,collapse_SNP(HipSTR.filter[,i],HipSTR.filter[,c(i+40)]))
}
# caculate genotype number
genoNO <- matrix(ncol=1,nrow=15411)
for(i in 1:15411){
genoNO[i,1] <- length(which(!is.na(GB.cl[i,])))
}

# generate consensus REF
Mode <- function(x) {x<-x[complete.cases(x)]; ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]}
cons.REF <-apply(GB.cl,1,Mode)


# caculate mutant number
MutantNO <- matrix(ncol=1,nrow=15411)
for(i in 1:15411){
 MutantNO[i,1] <- length(which(GB.cl[i,]!= as.character(cons.REF[i])))
}

HipSTR.final <- cbind(HipSTR.filter[,1:7],cons.REF,genoNO,MutantNO, GB.cl)
HipSTR.final2 <- subset(HipSTR.final, genoNO>3 & MutantNO>0)
write.csv(HipSTR.final2, "HipSTR.final2.csv")
