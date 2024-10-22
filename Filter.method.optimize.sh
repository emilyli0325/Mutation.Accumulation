# filter method optimize
########################################################################################################################
setwd("C:/Users/xli/Desktop/06062018/filterMethodsOptimize")


Input <- read.delim("GATK.Indel.hardfilter.07182018.table", sep="\t",header=T)
dim(Input)
[1] 11409   248

GT <- Input[,seq(6,248,3)]
GT[GT== "."]<-NA
AD <- Input[,seq(7,248,3)]
PL <- Input[,seq(8,248,3)]

> dim(GT)
[1] 11409    81

# individual consensus
maxA <- function(X){max(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
sumA <- function(X){sum(as.numeric(unlist(strsplit(as.character(X),",")[[1]])))}
DPmaxA <- matrix(ncol=81,nrow=11409)
DPsumA <- matrix(ncol=81,nrow=11409)

#purity
consIndi <- matrix(ncol=81,nrow=11409)
for(j in 1:81){
	DPmaxA[,j]<-sapply(AD[,j],maxA,simplify="array")
	DPsumA[,j]<-sapply(AD[,j],sumA,simplify="array")
	consIndi[,j]<- DPmaxA[,j]/DPsumA[,j]
	}

#PLdiff	
PLIndi <- matrix(ncol=81,nrow=11409)
for(j in 1:81){
	PLIndi[,j]<- sapply(PL[,j],maxA,simplify="array")
	}

GT <- GT[,2:81]
consIndi <- consIndi[,2:81]
PLIndi <- PLIndi[,2:81]

	
loci.true <-function(X,Y){
X<-as.character(X); Y<-as.character(Y);
Z<-NULL; 
Z <- ifelse(is.na(X), X, Y)
Z[which(X!=Y & is.na(Y)==F)]<-NA; 
return(Z)
}
	
loci.false <-function(X,Y){
X<-as.character(X); Y<-as.character(Y);
Z<-NULL; 
Z <- ifelse(is.na(X), X, Y)
Z[which(X==Y & is.na(Y)==F)]<-NA; 
return(Z)
}

GT.true <- NULL	
for(i in c(1:40)){
    GT.true <- cbind(GT.true,loci.true(GT[,i],GT[,c(i+40)]))
}	
GT.true<-cbind(GT.true,GT.true)

GT.false <- NULL	
for(i in c(1:40)){
    GT.false <- cbind(GT.false,loci.false(GT[,i],GT[,c(i+40)]))
}	
GT.false<-cbind(GT.false,GT.false)

###########
#true
consIndi.true <- consIndi
consIndi.true[is.na(GT.true)==T] <- NA
PLIndi.true <- PLIndi
PLIndi.true[is.na(GT.true)==T] <- NA
cons.PL.true <- NULL
cons.PL.true <- cbind(as.numeric(consIndi.true),as.numeric(PLIndi.true))
cons.PL.true.1 <- cons.PL.true[rowSums(!is.na(cons.PL.true)) > 0, ]
colnames(cons.PL.true.1)<- c("cons", "PL")
head(cons.PL.true.1)

cons.PL.true.table <- count(as.data.frame(cons.PL.true.1))
cons.PL.true.table <- cons.PL.true.table[order(cons.PL.true.table$freq),]
> dim(cons.PL.true.table)
[1] 87598     3
plot(cons~PL, data=cons.PL.true.table[1:86000,], xlim=range(0:5000), pch=16, cex = freq/20, col= "lightgrey")

#false
PLIndi.false <- PLIndi
PLIndi.false[is.na(GT.false)==T] <- NA

consIndi.false <- consIndi
consIndi.false[is.na(GT.false)==T] <- NA

cons.PL.false <- NULL
cons.PL.false <- cbind(as.numeric(consIndi.false[,1:40]),as.numeric(PLIndi.false[,1:40]),as.numeric(consIndi.false[,41:80]),as.numeric(PLIndi.false[,41:80]))
cons.PL.false.1 <- cons.PL.false[rowSums(!is.na(cons.PL.false)) > 0, ]
colnames(cons.PL.false.1)<- c("cons1", "PL1", "cons2", "PL2")

cons.PL.false.1 <- as.data.frame(cons.PL.false.1)
cons.PL.false.2 <- cons.PL.false.1[with(cons.PL.false.1,which(PL1 <= PL2)), ][,1:2]
cons.PL.false.3 <- cons.PL.false.1[with(cons.PL.false.1,which(PL1 > PL2)), ][,3:4]
colnames(cons.PL.false.2)<- c("cons", "PL")
colnames(cons.PL.false.3)<- c("cons", "PL")
cons.PL.false.4 <- rbind(cons.PL.false.2,cons.PL.false.3)

lapply(as.data.frame(cons.PL.false.4), class)
library(dplyr)
cons.PL.false.table <- count(cons.PL.false.4)
lapply(cons.PL.false.table, class)
cons.PL.false.table <- cons.PL.false.table[rev(order(cons.PL.false.table$freq)),]

points(cons~PL, data=cons.PL.false.table, xlim=range(0:5000), cex = freq/20, pch=16, col= "red")
abline(v=200,col="black", lty=2, lwd=2)



pdf('test.pdf', width=5, height=5)
plot(cons~PL, data=cons.PL.true.table[1:86000,], xlim=range(0:5000), pch=16, cex = freq/20, col= "lightgrey")
points(cons~PL, data=cons.PL.false.table, xlim=range(0:5000), cex = freq/20, pch=16, col= "red")
abline(v=250,col="black", lty=2, lwd=1.5)
abline(h=0.8,col="black", lty=2, lwd=1.5)
dev.off()



########################################################################################################################


cons.PL.true.table <- rename(count(as.data.frame(cons.PL.true.1), cons, PL), Freq = n)



test <- cons.PL.true.table[c(100000:100010),]





library(plyr)
counts <- ddply(data.frame(cons.PL.true.1), .(data.frame(cons.PL.true.1)$PL, data.frame(cons.PL.true.1)$cons), nrow)
names(counts) <- c("PL", "cons", "Freq")



plot(cons~PL, data=cons.PL.true.table, cex = as.numeric(Freq))




loci.false <-function(X,Y){
X<-as.character(X); Y<-as.character(Y);
Z<-NULL; 
Z <- ifelse(is.na(X), X, Y)
Z[which(X==Y & is.na(Y)==F)]<-NA; 
return(Z)
}



GB.cl<-NULL
for(i in c(9:48)){
    GB.cl<-cbind(GB.cl,collapse_SNP(filter[,i],filter[,c(i+40)]))
}

# caculate genotype number
genoNO <- matrix(ncol=1,nrow=10947)
for(i in 1:10947){
genoNO[i,1] <- length(which(!is.na(GB.cl[i,])))
}

# generate consensus REF
Mode <- function(x) {x<-x[complete.cases(x)]; ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]}
cons.REF <-apply(GB.cl,1,Mode)


# caculate mutant number
MutantNO <- matrix(ncol=1,nrow=10947)
for(i in 1:10947){
 MutantNO[i,1] <- length(which(GB.cl[i,]!= as.character(cons.REF[i])))
}

final <- subset(cbind(filter[,1:7],cons.REF,genoNO,MutantNO, GB.cl),genoNO>3 & MutantNO>0)

GT.filter <- GT.filter[genoNO>3 & MutantNO>0,]
C.afterrmHet <- matrix(ncol=2,nrow=42)
C.afterrmHet[1,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,20]))
C.afterrmHet[2,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,60]))
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT.filter[,j]), as.character(GT.filter[,j+40]))
}

write.csv(final, "GATK.Indel.final.withoutVQSR.csv")
