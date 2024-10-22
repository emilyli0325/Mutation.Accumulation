#/master/xli/software/mreps/mreps -exp 3 -minp 1 -maxp 9 -minsize 10 -maxsize 1000 -fasta PlasmoDB-32_Pfalciparum3D7_Genome.fasta > PlasmoDB-32_exp3_minp1_maxp9_minsize10_maxsize1000.txt


#1 haploid seperate
/master/xli/software/HipSTR/HipSTR --bams 3D7.recal.bam,MA3D7.recal.bam,MA12.recal.bam,MA18.recal.bam,MA19.recal.bam,MA20.recal.bam,MA21.recal.bam,MA22.recal.bam,MA23.recal.bam,MA24.recal.bam,MA26.recal.bam,MA27.recal.bam,MA29.recal.bam,MA30.recal.bam,MA32.recal.bam,MA35.recal.bam,MA36.recal.bam,MA37.recal.bam,MA38.recal.bam,MA39.recal.bam,MA40.recal.bam,MA41.recal.bam,MA42.recal.bam,MA43.recal.bam,MA44.recal.bam,MA45.recal.bam,MA46.recal.bam,MA47.recal.bam,MA49.recal.bam,MA51.recal.bam,MA52.recal.bam,MA53.recal.bam,MA54.recal.bam,MA55.recal.bam,MA56.recal.bam,MA57.recal.bam,MA58.recal.bam,MA6.recal.bam,MA60.recal.bam,MA61.recal.bam,MA62.recal.bam,Sample_12.recal.bam,Sample_18.recal.bam,Sample_19.recal.bam,Sample_20.recal.bam,Sample_21.recal.bam,Sample_22.recal.bam,Sample_23.recal.bam,Sample_24.recal.bam,Sample_26.recal.bam,Sample_27.recal.bam,Sample_29.recal.bam,Sample_30.recal.bam,Sample_32.recal.bam,Sample_35.recal.bam,Sample_36.recal.bam,Sample_37.recal.bam,Sample_38.recal.bam,Sample_39.recal.bam,Sample_3D7.recal.bam,Sample_40.recal.bam,Sample_41.recal.bam,Sample_42.recal.bam,Sample_43.recal.bam,Sample_44.recal.bam,Sample_45.recal.bam,Sample_46.recal.bam,Sample_47.recal.bam,Sample_49.recal.bam,Sample_51.recal.bam,Sample_52.recal.bam,Sample_53.recal.bam,Sample_54.recal.bam,Sample_55.recal.bam,Sample_56.recal.bam,Sample_57.recal.bam,Sample_58.recal.bam,Sample_6.recal.bam,Sample_60.recal.bam,Sample_61.recal.bam,Sample_62.recal.bam --fasta /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta --regions /master/xli/Index/HipSTR/Pf3D7_1-9_10-70.bed --str-vcf RUN_seperate_haploid.vcf.gz  --haploid-chrs Pf3D7_01_v3,Pf3D7_02_v3,Pf3D7_03_v3,Pf3D7_04_v3,Pf3D7_05_v3,Pf3D7_06_v3,Pf3D7_07_v3,Pf3D7_08_v3,Pf3D7_09_v3,Pf3D7_10_v3,Pf3D7_11_v3,Pf3D7_12_v3,Pf3D7_13_v3,Pf3D7_14_v3

#2 haploid merge
/master/xli/software/HipSTR/HipSTR --bams merge12.recal.bam,merge18.recal.bam,merge19.recal.bam,merge20.recal.bam,merge21.recal.bam,merge22.recal.bam,merge23.recal.bam,merge24.recal.bam,merge26.recal.bam,merge27.recal.bam,merge29.recal.bam,merge30.recal.bam,merge32.recal.bam,merge35.recal.bam,merge36.recal.bam,merge37.recal.bam,merge38.recal.bam,merge39.recal.bam,merge3D7.recal.bam,merge40.recal.bam,merge41.recal.bam,merge42.recal.bam,merge43.recal.bam,merge44.recal.bam,merge45.recal.bam,merge46.recal.bam,merge47.recal.bam,merge49.recal.bam,merge51.recal.bam,merge52.recal.bam,merge53.recal.bam,merge54.recal.bam,merge55.recal.bam,merge56.recal.bam,merge57.recal.bam,merge58.recal.bam,merge6.recal.bam,merge60.recal.bam,merge61.recal.bam,merge62.recal.bam  --fasta /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta --regions /master/xli/Index/HipSTR/Pf3D7_1-9_10-70.bed --str-vcf RUN_merge_haploid.vcf.gz  --haploid-chrs Pf3D7_01_v3,Pf3D7_02_v3,Pf3D7_03_v3,Pf3D7_04_v3,Pf3D7_05_v3,Pf3D7_06_v3,Pf3D7_07_v3,Pf3D7_08_v3,Pf3D7_09_v3,Pf3D7_10_v3,Pf3D7_11_v3,Pf3D7_12_v3,Pf3D7_13_v3,Pf3D7_14_v3

#3 diploid seperate
/master/xli/software/HipSTR/HipSTR --bams 3D7.recal.bam,MA3D7.recal.bam,MA12.recal.bam,MA18.recal.bam,MA19.recal.bam,MA20.recal.bam,MA21.recal.bam,MA22.recal.bam,MA23.recal.bam,MA24.recal.bam,MA26.recal.bam,MA27.recal.bam,MA29.recal.bam,MA30.recal.bam,MA32.recal.bam,MA35.recal.bam,MA36.recal.bam,MA37.recal.bam,MA38.recal.bam,MA39.recal.bam,MA40.recal.bam,MA41.recal.bam,MA42.recal.bam,MA43.recal.bam,MA44.recal.bam,MA45.recal.bam,MA46.recal.bam,MA47.recal.bam,MA49.recal.bam,MA51.recal.bam,MA52.recal.bam,MA53.recal.bam,MA54.recal.bam,MA55.recal.bam,MA56.recal.bam,MA57.recal.bam,MA58.recal.bam,MA6.recal.bam,MA60.recal.bam,MA61.recal.bam,MA62.recal.bam,Sample_12.recal.bam,Sample_18.recal.bam,Sample_19.recal.bam,Sample_20.recal.bam,Sample_21.recal.bam,Sample_22.recal.bam,Sample_23.recal.bam,Sample_24.recal.bam,Sample_26.recal.bam,Sample_27.recal.bam,Sample_29.recal.bam,Sample_30.recal.bam,Sample_32.recal.bam,Sample_35.recal.bam,Sample_36.recal.bam,Sample_37.recal.bam,Sample_38.recal.bam,Sample_39.recal.bam,Sample_3D7.recal.bam,Sample_40.recal.bam,Sample_41.recal.bam,Sample_42.recal.bam,Sample_43.recal.bam,Sample_44.recal.bam,Sample_45.recal.bam,Sample_46.recal.bam,Sample_47.recal.bam,Sample_49.recal.bam,Sample_51.recal.bam,Sample_52.recal.bam,Sample_53.recal.bam,Sample_54.recal.bam,Sample_55.recal.bam,Sample_56.recal.bam,Sample_57.recal.bam,Sample_58.recal.bam,Sample_6.recal.bam,Sample_60.recal.bam,Sample_61.recal.bam,Sample_62.recal.bam --fasta /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta --regions /master/xli/Index/HipSTR/Pf3D7_1-9_10-70.bed --str-vcf RUN.seperate.diploid.vcf.gz


#4 diploid merge
/master/xli/software/HipSTR/HipSTR --bams merge12.recal.bam,merge18.recal.bam,merge19.recal.bam,merge20.recal.bam,merge21.recal.bam,merge22.recal.bam,merge23.recal.bam,merge24.recal.bam,merge26.recal.bam,merge27.recal.bam,merge29.recal.bam,merge30.recal.bam,merge32.recal.bam,merge35.recal.bam,merge36.recal.bam,merge37.recal.bam,merge38.recal.bam,merge39.recal.bam,merge3D7.recal.bam,merge40.recal.bam,merge41.recal.bam,merge42.recal.bam,merge43.recal.bam,merge44.recal.bam,merge45.recal.bam,merge46.recal.bam,merge47.recal.bam,merge49.recal.bam,merge51.recal.bam,merge52.recal.bam,merge53.recal.bam,merge54.recal.bam,merge55.recal.bam,merge56.recal.bam,merge57.recal.bam,merge58.recal.bam,merge6.recal.bam,merge60.recal.bam,merge61.recal.bam,merge62.recal.bam  --fasta /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta --regions /master/xli/Index/HipSTR/Pf3D7_1-9_10-70.bed --str-vcf RUN.merge.diploid.vcf.gz


python /master/xli/software/HipSTR/scripts/filter_vcf.py --vcf RUN.seperate.diploid.vcf.gz --min-call-qual 0.8 --max-call-flank-indel 0.15 --max-call-stutter 0.15 --min-call-spanning-depth 2 --min-loc-calls 3 >RUN.seperate.diploid.filter.vcf

python /master/xli/software/HipSTR/scripts/filter_vcf.py --vcf RUN.merge.diploid.vcf.gz --min-call-qual 0.8 --max-call-flank-indel 0.15 --max-call-stutter 0.15 --min-call-spanning-depth 2 --min-loc-calls 3 >RUN.merge.diploid.filter.vcf

python /master/xli/software/HipSTR/scripts/filter_haploid_vcf.py --vcf RUN_seperate_haploid.vcf.gz --max-call-flank-indel 0.15 --max-call-stutter 0.15 --min-call-qual 0.9 --no-spanning >RUN_seperate_haploid_filter.vcf

python /master/xli/software/HipSTR/scripts/filter_haploid_vcf.py --vcf RUN_merge_haploid.vcf.gz --max-call-flank-indel 0.15 --max-call-stutter 0.15 --min-call-qual 0.9 --no-spanning >RUN_merge_haploid_filter.vcf


java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V RUN_merge_haploid_filter.vcf \
-F CHROM -F POS -F ID -F REF -F ALT -GF GT -GF GB -GF Q -GF MALLREADS \
-o RUN_merge_haploid_filter.table


/master/xli/software/HipSTR/HipSTR --bams 3D7.recal.bam,MA3D7.recal.bam,MA12.recal.bam,MA18.recal.bam,MA19.recal.bam,MA20.recal.bam,MA21.recal.bam,MA22.recal.bam,MA23.recal.bam,MA24.recal.bam,MA26.recal.bam,MA27.recal.bam,MA29.recal.bam,MA30.recal.bam,MA32.recal.bam,MA35.recal.bam,MA36.recal.bam,MA37.recal.bam,MA38.recal.bam,MA39.recal.bam,MA40.recal.bam,MA41.recal.bam,MA42.recal.bam,MA43.recal.bam,MA44.recal.bam,MA45.recal.bam,MA46.recal.bam,MA47.recal.bam,MA49.recal.bam,MA51.recal.bam,MA52.recal.bam,MA53.recal.bam,MA54.recal.bam,MA55.recal.bam,MA56.recal.bam,MA57.recal.bam,MA58.recal.bam,MA6.recal.bam,MA60.recal.bam,MA61.recal.bam,MA62.recal.bam,Sample_12.recal.bam,Sample_18.recal.bam,Sample_19.recal.bam,Sample_20.recal.bam,Sample_21.recal.bam,Sample_22.recal.bam,Sample_23.recal.bam,Sample_24.recal.bam,Sample_26.recal.bam,Sample_27.recal.bam,Sample_29.recal.bam,Sample_30.recal.bam,Sample_32.recal.bam,Sample_35.recal.bam,Sample_36.recal.bam,Sample_37.recal.bam,Sample_38.recal.bam,Sample_39.recal.bam,Sample_3D7.recal.bam,Sample_40.recal.bam,Sample_41.recal.bam,Sample_42.recal.bam,Sample_43.recal.bam,Sample_44.recal.bam,Sample_45.recal.bam,Sample_46.recal.bam,Sample_47.recal.bam,Sample_49.recal.bam,Sample_51.recal.bam,Sample_52.recal.bam,Sample_53.recal.bam,Sample_54.recal.bam,Sample_55.recal.bam,Sample_56.recal.bam,Sample_57.recal.bam,Sample_58.recal.bam,Sample_6.recal.bam,Sample_60.recal.bam,Sample_61.recal.bam,Sample_62.recal.bam --fasta /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta --regions HipSTR.vis.bed --str-vcf HipSTR.vcf.gz  --haploid-chrs Pf3D7_01_v3,Pf3D7_02_v3,Pf3D7_03_v3,Pf3D7_04_v3,Pf3D7_05_v3,Pf3D7_06_v3,Pf3D7_07_v3,Pf3D7_08_v3,Pf3D7_09_v3,Pf3D7_10_v3,Pf3D7_11_v3,Pf3D7_12_v3,Pf3D7_13_v3,Pf3D7_14_v3 --viz-out HipSTR.vis.vcf.gz

tabix -p HipSTR.vis.bed HipSTR.vis.vcf.gz

/master/xli/software/HipSTR/VizAlnPdf HipSTR.vis.vcf.gz Pf3D7_03_v3 613596 MA21/  test.pdf 1




GT:      ALLREADS:DFLANKINDEL:DP:DSTUTTER:GB:  GLDIFF:     MALLREADS:Q
 0:           0|4:          0:13:       0: 0:    3.99:           0|4:1.0
 0: -2|1;0|18;2|2:          0:45:       3: 0:   17.11: -2|1;0|15;2|2:1.0
 0:-4|2;-2|5;0|16:          1:88:       9: 0:    3.99:     -2|4;0|13:1.0
 
GT:GB:  Q:DP:DSTUTTER:DFLANKINDEL:GLDIFF:      ALLREADS :MALLREADS
 0: 0:1.0:88:       9:          1:  3.99:-4|2;-2|5;0|16 :-2|4;0|13
 0: 0:1.0:85:       1:          3:  4.51:-50|1;-2|1;0|12:-2|1;0|11
 
 

python /master/xli/software/HipSTR/scripts/filter_haploid_vcf.py --vcf RUN_seperate_haploid.vcf.gz --max-call-flank-indel 0.15 --max-call-stutter 0.15 --min-call-qual 0.9 --no-spanning > HipSTR_basicfilter.vcf

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V MA.HipSTR.basicfilter.vcf \
-F ID -F CHROM -F REF -F ALT -F POS -F BPDIFFS -GF GT -GF MALLREADS -GF GLDIFF \
-o MA.HipSTR.basicfilter.table2

awk '$6 != 0 {print $0}' MA.HipSTR.basicfilter.table2 > MA.HipSTR.basicfilter.table3


python /master/xli/software/HipSTR/scripts/filter_haploid_vcf.py --vcf RUN_seperate_haploid.vcf.gz --max-call-flank-indel 0.10 --max-call-stutter 0.10 --min-call-qual 0.9 --no-spanning > HipSTR_basicfilter_2.vcf

awk '$5 != "." {print $0}' HipSTR_basicfilter_2.vcf > HipSTR_basicfilter_3.vcf

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V HipSTR_basicfilter_3.vcf \
-F ID -F CHROM -F REF -F ALT -F POS -F BPDIFFS -GF GT -GF MALLREADS -GF GLDIFF \
-o HipSTR_basicfilter_3.table

awk '$6 != 0 {print $0}' HipSTR_basicfilter_3.table > HipSTR_basicfilter_4.table


awk '$5 != "." {print $0}' RUN_seperate_haploid.vcf > RUN_seperate_haploid_withALT.vcf
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V RUN_seperate_haploid_withALT.vcf \
-F ID -F CHROM -F REF -F ALT -F POS -F BPDIFFS -GF GT -GF MALLREADS -GF GLDIFF \
-o RUN_seperate_haploid_withALT.table
awk '$6 != 0 {print $0}' RUN_seperate_haploid_withALT.table > HipSTR_raw.table




########################################################################################################################
setwd("C:/Users/xli/Desktop/Table3")
HipSTR <- read.delim("HipSTR_basicfilter_4.table", sep="\t",header=T)
dim(HipSTR)
[1] 17649   249

GT <- HipSTR[,seq(7,249,3)]
GT[GT== "."]<-NA

concordance <- function(X,Y){a1<-length(which(X==Y)); a2<-length(which(X!=Y)); a3<-cbind(a1,a2); return(a3)}
C.afterrmHet <- matrix(ncol=2,nrow=42)
C.afterrmHet[1,]<-concordance(as.character(GT[,1]),as.character(GT[,20]))
C.afterrmHet[2,]<-concordance(as.character(GT[,1]),as.character(GT[,60]))
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT[,j]), as.character(GT[,j+40]))
}


AD <- HipSTR[,seq(8,249,3)]
GLdiff <- HipSTR[,seq(9,249,3)]

> dim(AD)
[1] 17649    81


# individual consensus
maxA <- function(X){max(as.numeric(unlist(strsplit(strsplit(as.character(X),";")[[1]],"\\|"))[!c(TRUE,FALSE)]))}
sumA <- function(X){sum(as.numeric(unlist(strsplit(strsplit(as.character(X),";")[[1]],"\\|"))[!c(TRUE,FALSE)]))}

DPmaxA <- matrix(ncol=81,nrow=17649)
DPsumA <- matrix(ncol=81,nrow=17649)


consIndi <- matrix(ncol=81,nrow=17649)
for(j in 1:81){
	DPmaxA[,j]<-sapply(AD[,j],maxA,simplify="array")
	DPsumA[,j]<-sapply(AD[,j],sumA,simplify="array")
	consIndi[,j]<- DPmaxA[,j]/DPsumA[,j]
	}
	
	
consIndi.fail.percen <- matrix(ncol=1,nrow=17649)
for(i in 1:17649){
consIndi.fail.percen[i,] <- length(which((consIndi[i,])<0.8))/length(which(is.na(consIndi[i,])== F)) 
}


GT[consIndi<0.8]<-NA
GT[DPmaxA<3]<-NA
GT[GLdiff < 5 & consIndi == 1]<-NA
GT[GLdiff < 10 & consIndi < 1]<-NA

GT.filter <- GT[consIndi.fail.percen<0.5,]
> dim(GT.filter)
[1] 15975    81


C.afterrmHet <- matrix(ncol=2,nrow=42)

C.afterrmHet[1,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,20]))
C.afterrmHet[2,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,60]))
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT.filter[,j]), as.character(GT.filter[,j+40]))
}

HipSTR.filter <- cbind(HipSTR[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)

write.csv(HipSTR.filter, "HipSTR.filter.csv")
# HipSTR cons merge
# HipSTR.filter <- cbind(HipSTR[consIndi.fail.percen<0.5,][,1:6], consIndi.fail.percen[consIndi.fail.percen<0.5,], GT.filter)

dim(HipSTR.filter)
[1] 15975    88

collapse_SNP<-function(X,Y){
X<-as.character(X); Y<-as.character(Y);
Z<-NULL; 
Z <- ifelse(is.na(X), Y, X)
Z[which(X!=Y & is.na(X)==F & is.na(Y)==F)]<-NA; 
return(Z)
}


GB.cl<-NULL
for(i in c(9:48)){
    GB.cl<-cbind(GB.cl,collapse_SNP(HipSTR.filter[,i],HipSTR.filter[,c(i+40)]))
}
# caculate genotype number
genoNO <- matrix(ncol=1,nrow=15975)
for(i in 1:15975){
genoNO[i,1] <- length(which(!is.na(GB.cl[i,])))
}

# generate consensus REF
Mode <- function(x) {x<-x[complete.cases(x)]; ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]}
cons.REF <-apply(GB.cl,1,Mode)


# caculate mutant number
MutantNO <- matrix(ncol=1,nrow=15975)
for(i in 1:15975){
 MutantNO[i,1] <- length(which(GB.cl[i,]!= as.character(cons.REF[i])))
}

HipSTR.final <- cbind(HipSTR.filter[,1:7],cons.REF,genoNO,MutantNO, GB.cl)
HipSTR.final2 <- subset(HipSTR.final, genoNO>3 & MutantNO>0)
write.csv(HipSTR.final2, "HipSTR.final_2.csv")


GT.filter <- GT[consIndi.fail.percen<0.5,]
GT.filter <- GT.filter[genoNO>3 & MutantNO>0,]
C.afterrmHet <- matrix(ncol=2,nrow=42)
C.afterrmHet[1,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,20]))
C.afterrmHet[2,]<-concordance(as.character(GT.filter[,1]),as.character(GT.filter[,60]))
for(j in seq(2,41,1)){
C.afterrmHet[j+1,] <- concordance(as.character(GT.filter[,j]), as.character(GT.filter[,j+40]))
}




/master/xli/software/HipSTR/HipSTR --bams 3D7.recal.bam,MA3D7.recal.bam,MA12.recal.bam,MA18.recal.bam,MA19.recal.bam,MA20.recal.bam,MA21.recal.bam,MA22.recal.bam,MA23.recal.bam,MA24.recal.bam,MA26.recal.bam,MA27.recal.bam,MA29.recal.bam,MA30.recal.bam,MA32.recal.bam,MA35.recal.bam,MA36.recal.bam,MA37.recal.bam,MA38.recal.bam,MA39.recal.bam,MA40.recal.bam,MA41.recal.bam,MA42.recal.bam,MA43.recal.bam,MA44.recal.bam,MA45.recal.bam,MA46.recal.bam,MA47.recal.bam,MA49.recal.bam,MA51.recal.bam,MA52.recal.bam,MA53.recal.bam,MA54.recal.bam,MA55.recal.bam,MA56.recal.bam,MA57.recal.bam,MA58.recal.bam,MA6.recal.bam,MA60.recal.bam,MA61.recal.bam,MA62.recal.bam,Sample_12.recal.bam,Sample_18.recal.bam,Sample_19.recal.bam,Sample_20.recal.bam,Sample_21.recal.bam,Sample_22.recal.bam,Sample_23.recal.bam,Sample_24.recal.bam,Sample_26.recal.bam,Sample_27.recal.bam,Sample_29.recal.bam,Sample_30.recal.bam,Sample_32.recal.bam,Sample_35.recal.bam,Sample_36.recal.bam,Sample_37.recal.bam,Sample_38.recal.bam,Sample_39.recal.bam,Sample_3D7.recal.bam,Sample_40.recal.bam,Sample_41.recal.bam,Sample_42.recal.bam,Sample_43.recal.bam,Sample_44.recal.bam,Sample_45.recal.bam,Sample_46.recal.bam,Sample_47.recal.bam,Sample_49.recal.bam,Sample_51.recal.bam,Sample_52.recal.bam,Sample_53.recal.bam,Sample_54.recal.bam,Sample_55.recal.bam,Sample_56.recal.bam,Sample_57.recal.bam,Sample_58.recal.bam,Sample_6.recal.bam,Sample_60.recal.bam,Sample_61.recal.bam,Sample_62.recal.bam --fasta /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta --regions HipSTR.vis.bed --str-vcf HipSTR.vcf.gz --haploid-chrs Pf3D7_01_v3,Pf3D7_02_v3,Pf3D7_03_v3,Pf3D7_04_v3,Pf3D7_05_v3,Pf3D7_06_v3,Pf3D7_07_v3,Pf3D7_08_v3,Pf3D7_09_v3,Pf3D7_10_v3,Pf3D7_11_v3,Pf3D7_12_v3,Pf3D7_13_v3,Pf3D7_14_v3 --viz-out HipSTR.vis.vcf.gz

tabix HipSTR.vis.vcf.gz

/master/xli/software/HipSTR/VizAlnPdf HipSTR.vis.vcf.gz Pf3D7_03_v3 613596 MA21/ test.pdf 1






