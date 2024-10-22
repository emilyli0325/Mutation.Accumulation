### 06112019 MA review False negative ####

### Synthetic genome with 2000 SNPs and 3000 microsatellite mutations ####
# mutations are from P01 parents genotypes: /data/infectious/malaria_XUE/cross/Parents/GT.P

# SNPs: Biallele, called in all parents, core genome, not inside microsatellite
# Microsatellite: (1) Biallele, called in all parents, core genome, non-SNP, simple mutation 
#                 (2) 1-9bp in motif, 10-70bp in length, indel distribution same as Figure 7B
# Samples: 31 run1, 26 run2, 3 3D7

export PATH=$PATH:/master/xli/software/bedtools2/bin
export PATH=$PATH:/master/xli/software/vcftools/bin	
export PATH=$PATH:/master/xli/software/tabix-0.2.6
export PERL5LIB=$PERL5LIB:/master/xli/software/vcftools/src/perl
export PATH=$PATH:/master/xli/software/samtools-1.3/

############################################################################################################
## 2000 SNPs
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
-V Parent.SNP.VQSR.vcf \
--excludeNonVariants \
--maxNOCALLnumber 0 \
--removeUnusedAlternates \
--restrictAllelesTo BIALLELIC \
-o MA.dummy.ref/Parent.core.bi.SNP.VQSR.vcf

/master/xli/software/FreeBayes/vcflib/bin/vcffilter -f "FILTER = PASS" Parent.core.bi.SNP.VQSR.vcf > MA.dummy.SNP.candi.vcf	

bedtools intersect -a MA.dummy.SNP.candi.vcf -b /master/xli/Index/Genome_region/FINAL/Mocro.1-9.10-1000.wholeGenome.bed -wa -v | shuf -n 2000 >MA.dummy.SNP.2000.vcf
cat MA.dummy.SNP.candi.vcf | grep '^#' >header.txt
cat header.txt MA.dummy.SNP.2000.vcf >MA.dummy.SNP.2000.header.vcf
vcf-sort MA.dummy.SNP.2000.header.vcf > MA.dummy.SNP.2000.header.sort.vcf

## 3000 microsatellites
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
-V Parent.INDEL.VQSR.vcf \
--excludeNonVariants \
--maxNOCALLnumber 0 \
--removeUnusedAlternates \
--restrictAllelesTo BIALLELIC \
-o MA.dummy.ref/Parent.core.bi.INDEL.VQSR.vcf

/master/xli/software/FreeBayes/vcflib/bin/vcffilter -f "FILTER = PASS" Parent.core.bi.INDEL.VQSR.vcf > MA.dummy.INDEL.candi.vcf
bedtools intersect -a MA.dummy.INDEL.candi.vcf -b /master/xli/Index/Genome_region/FINAL/Mocro.1-9.10-1000.coreGenome.start-1.bed -wa | sort | uniq -u >MA.dummy.INDEL.all.vcf
cat MA.dummy.INDEL.candi.vcf | grep '^#' >header.Indel.txt
cat header.Indel.txt MA.dummy.INDEL.all.vcf >MA.dummy.INDEL.all.header.vcf
vcf-sort MA.dummy.INDEL.all.header.vcf > MA.dummy.INDEL.all.header.sort.vcf

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V MA.dummy.INDEL.all.header.sort.vcf \
-F CHROM -F POS -F REF -F ALT \
-o MA.dummy.INDEL.all.header.sort.table
		
### Use Access to get only indels at the start of a Microsatellite
8350 candidates 

### Generate random mutations ## in R
setwd("C:/Users/xli/Desktop/2019_paper_submission/MAforGBE_0320/review1/FalseNegative")
candi <- read.csv("Dummy.Candi.csv",header=T)
candi <- cbind(candi, floor(candi$exp))
colnames(candi)[19] <- "expZ"
> dim(candi)
[1] 8350   19

Indel <- matrix(ncol=1,nrow=8350)
for(j in 1:8350){
    n <- candi$expZ[j]
	i <- n*2
	m <- rbinom(1,i,.5)
	while (m == n) { m <- rbinom(1,i,.5) }		   
    Indel[j,1] <- m     
	}

candi <- cbind(candi[,1:19], Indel)
IndelSize <- candi$Indel - candi$expZ
IndelSize2 <- abs(IndelSize)
candi <- cbind(candi[,1:19], Indel, IndelSize,IndelSize2)

IndelSeq <- matrix(ncol=1,nrow=8350)
for(j in 1:8350){
      Motif <- candi$RepeatUnit[j]
	  i <- candi$IndelSize2[j]
      IndelSeq[j,1] <- paste(rep(Motif, times = i),collapse="") 
}

candi <- cbind(candi[,1:19], Indel, IndelSize,IndelSize2, IndelSeq)
write.csv(candi,"candi.csv")

pdf('IndelSize.pdf', width=6, height=6)
hist(IndelSize, pch=20, breaks=25, prob=T, main="",col="orange", xlim=c(-15,15), xlab="Repeat number change", ylab="Frequency")
dev.off()

#### make dummy REF and ALT, remove loci that too close to each other 
candi.filter <- read.csv("candi-filter.csv",header=T)

bp1 <- candi.filter[candi.filter$per == 1, ]
bp2 <- candi.filter[candi.filter$per == 2, ]
bp3 <- candi.filter[candi.filter$per == 3, ]
bp4 <- candi.filter[candi.filter$per == 4, ]
bp5 <- candi.filter[candi.filter$per == 5, ]
bp6 <- candi.filter[candi.filter$per == 6, ]
bp7 <- candi.filter[candi.filter$per == 7, ]
bp8 <- candi.filter[candi.filter$per == 8, ]
bp9 <- candi.filter[candi.filter$per == 9, ]

candi.filter2 <- rbind(bp1[sample(nrow(bp1), 838),], bp2[sample(nrow(bp2), 773),],bp3[sample(nrow(bp3), 392),],bp4[sample(nrow(bp4), 350),],bp5[sample(nrow(bp5), 200),],bp6[sample(nrow(bp6), 100),],bp7[sample(nrow(bp7), 50),],bp8[sample(nrow(bp8), 50),],bp9[sample(nrow(bp9), 50),]) 

write.csv(candi.filter2,"candi-filter2.csv")

#### merge SNP and Indel

java -cp /master/xli/software/GATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V MA.dummy.SNP.2000.header.sort.vcf \
-V MA.dummy.INDEL.3000.filter.vcf \
-out MA.dummy.vcf

vcf-sort MA.dummy.vcf > MA.dummy.sort.vcf

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V MA.dummy.sort.vcf \
-F CHROM -F POS -F REF -F ALT \
-o MA.dummy.sort.table


#### generate dummy reference
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T FastaAlternateReferenceMaker \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-o MA.dummy.fasta \
-V MA.dummy.sort.vcf

# change header in MA.dummy.fasta

# build Index
cp MA.dummy.fasta /data/infectious/malaria_XUE/MA/Genome_sequence_merge_two_run/FalseNegative/dummyRef
samtools faidx MA.dummy.fasta
java -jar /master/xli/software/picard/picard.jar CreateSequenceDictionary R=MA.dummy.fasta O=MA.dummy.dict
/master/xli/software/bwa3/bwa-0.7.15/bwa index MA.dummy.fasta

# mapping and genotyping as usual

# GATK done
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs -R /data/infectious/malaria_XUE/MA/Genome_sequence_merge_two_run/FalseNegative/dummyRef/MA.dummy.fasta \
-V 3D7.g.vcf -V AB-BSA-220-NF54.g.vcf -V AB-BSA-222-NHP4026.g.vcf -V MA12.g.vcf -V MA18.g.vcf -V MA19.g.vcf -V MA21.g.vcf -V MA22.g.vcf -V MA24.g.vcf -V MA26.g.vcf -V MA29.g.vcf -V MA30.g.vcf -V MA35.g.vcf -V MA36.g.vcf -V MA37.g.vcf -V MA38.g.vcf -V MA39.g.vcf -V MA3D7.g.vcf -V MA42.g.vcf -V MA43.g.vcf -V MA44.g.vcf -V MA46.g.vcf -V MA47.g.vcf -V MA49.g.vcf -V MA51.g.vcf -V MA53.g.vcf -V MA54.g.vcf -V MA55.g.vcf -V MA57.g.vcf -V MA58.g.vcf -V MA60.g.vcf -V MA61.g.vcf -V MA62.g.vcf -V MA6.g.vcf -V MKK2835.2G.g.vcf -V NHP1337.12C.g.vcf -V NHP4026-2G.g.vcf -V Sample.12.g.vcf -V Sample.18.g.vcf -V Sample.19.g.vcf -V Sample.21.g.vcf -V Sample.24.g.vcf -V Sample.26.g.vcf -V Sample.27.g.vcf -V Sample.29.g.vcf -V Sample.30.g.vcf -V Sample.35.g.vcf -V Sample.36.g.vcf -V Sample.37.g.vcf -V Sample.38.g.vcf -V Sample.39.g.vcf -V Sample.3D7.g.vcf -V Sample.42.g.vcf -V Sample.43.g.vcf -V Sample.44.g.vcf -V Sample.46.g.vcf -V Sample.49.g.vcf -V Sample.51.g.vcf -V Sample.53.g.vcf -V Sample.54.g.vcf -V Sample.55.g.vcf -V Sample.60.g.vcf -V Sample.61.g.vcf -V Sample.62.g.vcf -V Sample.6.g.vcf \
--useNewAFCalculator \
--sample_ploidy 1 \
-o FN.GATK.vcf

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /data/infectious/malaria_XUE/MA/Genome_sequence_merge_two_run/FalseNegative/dummyRef/MA.dummy.fasta \
-V FN.GATK.vcf \
-L TargetAll.intervals \
-o FN.GATK.target.vcf



# freeBayes
	export PATH=$PATH:/master/xli/software/FreeBayes/freebayes/bin
	export PATH=$PATH:/master/xli/software/FreeBayes/vcflib/bin
	export PATH=$PATH:/master/xli/software/vt

	
freebayes -f /data/infectious/malaria_XUE/MA/Genome_sequence_merge_two_run/FalseNegative/dummyRef/MA.dummy.fasta \
		  -L bam.list \
		  -t /master/xli/Index/Known_sites/Core.bed \
		  -v FN.freebayes.vcf \
		  -p 1 \
		  --min-repeat-size 3 \
		  --no-partial-observations \
		  -m 30 \
		  -q 20 

freebayes -f /data/infectious/malaria_XUE/MA/Genome_sequence_merge_two_run/FalseNegative/dummyRef/MA.dummy.fasta \
		  -L bam.list \
		  -t TargetAll-simple.bed \
		  -v FN.freebayes-target.vcf \
		  -p 1 \
		  --min-repeat-size 3 \
		  --no-partial-observations \
		  -m 30 \
		  -q 20 
		  		  
vcfallelicprimitives -kg FN.freebayes.vcf | /master/xli/software/vt/vt normalize -r /data/infectious/malaria_XUE/MA/Genome_sequence_merge_two_run/FalseNegative/dummyRef/MA.dummy.fasta -o FN.freebayes.norm.vcf -

# HipSTR
/master/xli/software/HipSTR/HipSTR \
--bams 3D7.dedup.sorted.bam,AB-BSA-220-NF54.dedup.sorted.bam,AB-BSA-222-NHP4026.dedup.sorted.bam,MA12.dedup.sorted.bam,MA18.dedup.sorted.bam,MA19.dedup.sorted.bam,MA21.dedup.sorted.bam,MA22.dedup.sorted.bam,MA24.dedup.sorted.bam,MA26.dedup.sorted.bam,MA29.dedup.sorted.bam,MA30.dedup.sorted.bam,MA35.dedup.sorted.bam,MA36.dedup.sorted.bam,MA37.dedup.sorted.bam,MA38.dedup.sorted.bam,MA39.dedup.sorted.bam,MA3D7.dedup.sorted.bam,MA42.dedup.sorted.bam,MA43.dedup.sorted.bam,MA44.dedup.sorted.bam,MA46.dedup.sorted.bam,MA47.dedup.sorted.bam,MA49.dedup.sorted.bam,MA51.dedup.sorted.bam,MA53.dedup.sorted.bam,MA54.dedup.sorted.bam,MA55.dedup.sorted.bam,MA57.dedup.sorted.bam,MA58.dedup.sorted.bam,MA60.dedup.sorted.bam,MA61.dedup.sorted.bam,MA62.dedup.sorted.bam,MA6.dedup.sorted.bam,MKK2835.2G.dedup.sorted.bam,NHP1337.12C.dedup.sorted.bam,NHP4026-2G.dedup.sorted.bam,Sample.12.dedup.sorted.bam,Sample.18.dedup.sorted.bam,Sample.19.dedup.sorted.bam,Sample.21.dedup.sorted.bam,Sample.24.dedup.sorted.bam,Sample.26.dedup.sorted.bam,Sample.27.dedup.sorted.bam,Sample.29.dedup.sorted.bam,Sample.30.dedup.sorted.bam,Sample.35.dedup.sorted.bam,Sample.36.dedup.sorted.bam,Sample.37.dedup.sorted.bam,Sample.38.dedup.sorted.bam,Sample.39.dedup.sorted.bam,Sample.3D7.dedup.sorted.bam,Sample.42.dedup.sorted.bam,Sample.43.dedup.sorted.bam,Sample.44.dedup.sorted.bam,Sample.46.dedup.sorted.bam,Sample.49.dedup.sorted.bam,Sample.51.dedup.sorted.bam,Sample.53.dedup.sorted.bam,Sample.54.dedup.sorted.bam,Sample.55.dedup.sorted.bam,Sample.60.dedup.sorted.bam,Sample.61.dedup.sorted.bam,Sample.62.dedup.sorted.bam,Sample.6.dedup.sorted.bam \
--fasta /data/infectious/malaria_XUE/MA/Genome_sequence_merge_two_run/FalseNegative/dummyRef/MA.dummy.fasta \
--regions TargetMicro.bed \
--str-vcf FN.HipSTR.vcf.gz \
--haploid-chrs Pf3D7_01_v3,Pf3D7_02_v3,Pf3D7_03_v3,Pf3D7_04_v3,Pf3D7_05_v3,Pf3D7_06_v3,Pf3D7_07_v3,Pf3D7_08_v3,Pf3D7_09_v3,Pf3D7_10_v3,Pf3D7_11_v3,Pf3D7_12_v3,Pf3D7_13_v3,Pf3D7_14_v3

python /master/xli/software/HipSTR/scripts/filter_haploid_vcf.py --vcf FN.HipSTR.vcf.gz --max-call-flank-indel 0.15 --max-call-stutter 0.15 --min-call-qual 0.9 --no-spanning > FN.HipSTR.filter.vcf

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V FN.HipSTR.filter.vcf \
-F CHROM -F POS -F ID -F REF -F ALT -GF GT -GF GB -GF Q -GF MALLREADS \
-o FN.HipSTR.filter.table


## with "--min-reads 15 --def-stutter-model": NO big difference
#/master/xli/software/HipSTR/HipSTR \
#--bams 3D7.dedup.sorted.bam,AB-BSA-220-NF54.dedup.sorted.bam,AB-BSA-222-NHP4026.dedup.sorted.bam,MA12.dedup.sorted.bam,MA18.dedup.sorted.bam,MA19.dedup.sorted.bam,MA21.dedup.sorted.bam,MA22.dedup.sorted.bam,MA24.dedup.sorted.bam,MA26.dedup.sorted.bam,MA29.dedup.sorted.bam,MA30.dedup.sorted.bam,MA35.dedup.sorted.bam,MA36.dedup.sorted.bam,MA37.dedup.sorted.bam,MA38.dedup.sorted.bam,MA39.dedup.sorted.bam,MA3D7.dedup.sorted.bam,MA42.dedup.sorted.bam,MA43.dedup.sorted.bam,MA44.dedup.sorted.bam,MA46.dedup.sorted.bam,MA47.dedup.sorted.bam,MA49.dedup.sorted.bam,MA51.dedup.sorted.bam,MA53.dedup.sorted.bam,MA54.dedup.sorted.bam,MA55.dedup.sorted.bam,MA57.dedup.sorted.bam,MA58.dedup.sorted.bam,MA60.dedup.sorted.bam,MA61.dedup.sorted.bam,MA62.dedup.sorted.bam,MA6.dedup.sorted.bam,MKK2835.2G.dedup.sorted.bam,NHP1337.12C.dedup.sorted.bam,NHP4026-2G.dedup.sorted.bam,Sample.12.dedup.sorted.bam,Sample.18.dedup.sorted.bam,Sample.19.dedup.sorted.bam,Sample.21.dedup.sorted.bam,Sample.24.dedup.sorted.bam,Sample.26.dedup.sorted.bam,Sample.27.dedup.sorted.bam,Sample.29.dedup.sorted.bam,Sample.30.dedup.sorted.bam,Sample.35.dedup.sorted.bam,Sample.36.dedup.sorted.bam,Sample.37.dedup.sorted.bam,Sample.38.dedup.sorted.bam,Sample.39.dedup.sorted.bam,Sample.3D7.dedup.sorted.bam,Sample.42.dedup.sorted.bam,Sample.43.dedup.sorted.bam,Sample.44.dedup.sorted.bam,Sample.46.dedup.sorted.bam,Sample.49.dedup.sorted.bam,Sample.51.dedup.sorted.bam,Sample.53.dedup.sorted.bam,Sample.54.dedup.sorted.bam,Sample.55.dedup.sorted.bam,Sample.60.dedup.sorted.bam,Sample.61.dedup.sorted.bam,Sample.62.dedup.sorted.bam,Sample.6.dedup.sorted.bam \
#--fasta /data/infectious/malaria_XUE/MA/Genome_sequence_merge_two_run/FalseNegative/dummyRef/MA.dummy.fasta \
#--regions TargetMicro.bed \
#--str-vcf FN.HipSTR-minReads15.vcf.gz --min-reads 15 --def-stutter-model \
#--haploid-chrs Pf3D7_01_v3,Pf3D7_02_v3,Pf3D7_03_v3,Pf3D7_04_v3,Pf3D7_05_v3,Pf3D7_06_v3,Pf3D7_07_v3,Pf3D7_08_v3,Pf3D7_09_v3,Pf3D7_10_v3,Pf3D7_11_v3,Pf3D7_12_v3,Pf3D7_13_v3,Pf3D7_14_v3

#python /master/xli/software/HipSTR/scripts/filter_haploid_vcf.py --vcf FN.HipSTR-minReads15.vcf.gz --max-call-flank-indel 0.15 --max-call-stutter 0.15 --min-call-qual 0.9 --no-spanning > FN.HipSTR-minReads15.filter.vcf










############################################################################################################

bedtools intersect -a MA.dummy.INDEL.candi.vcf -b /master/xli/Index/Genome_region/FINAL/Mocro.1-9.10-1000.coreGenome.start-1.bed -wa | sort | uniq -u | shuf -n 3000 >MA.dummy.INDEL.3000.vcf
cat MA.dummy.INDEL.candi.vcf | grep '^#' >header.Indel.txt
cat header.Indel.txt MA.dummy.INDEL.3000.vcf >MA.dummy.INDEL.3000.header.vcf
vcf-sort MA.dummy.INDEL.3000.header.vcf > MA.dummy.INDEL.3000.header.sort.vcf



