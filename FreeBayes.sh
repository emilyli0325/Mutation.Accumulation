####################################
############# From Fred ############
####################################
# Calling
freebayes -f ~/data/sm_genome/sma_v5.0.chr.fa \
    -r Schisto_mansoni.Chr_6:1518000-1526000  \
    -b $(find data -name *.bam | tr "\n" " ") \
    --population pop \
    -q 20            \
    -m 30            \
    -! 4             \
    -@ sm_dbSNP_Smp_089320.vcf.gz > Sm.SN.NE.TZ.OM_HR9_Sr.pChr6.fb.raw.snp_indel.vcf

# Normalize VCF
vt normalize -o Sm.SN.NE.TZ.OM_HR9_Sr.pChr6.fb.raw.snp_indel.norm.vcf \
    -r ~/data/sm_genome/sma_v5.0.chr.fa \
    Sm.SN.NE.TZ.OM_HR9_Sr.pChr6.fb.raw.snp_indel.vcf
	

####### For MA ######
export PATH=$PATH:/master/xli/software/FreeBayes/freebayes/bin
freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
		  -L bam.indi.txt \
		  -t /master/xli/Index/Known_sites/Core.bed \
		  -v MA.indi.freebayes.vcf \
		  -p 1 \
		  --min-repeat-size 3 \
		  --no-partial-observations \
		  -m 30 \
		  -q 20 

freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
		  -L bam.merge.txt \
		  -t /master/xli/Index/Known_sites/Core.bed \
		  -v MA.merge.freebayes.vcf \
		  -p 1 \
		  --min-repeat-size 3 \
		  --no-partial-observations \
		  -m 30 \
		  -q 20 

f1
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_09_v3:79101-1242137 -v MA.indi.chr9-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20 
f2
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_09_v3:1244484-1473560 -v MA.indi.chr9-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20 
f3
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_10_v3:68971-1571815 -v MA.indi.chr10.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f4
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_11_v3:110001-831968 -v MA.indi.chr11-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20

f5
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_11_v3:834246-2003320 -v MA.indi.chr11-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20

f6
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_12_v3:60301-766654 -v MA.indi.chr12-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20

f7
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_12_v3:780451-1282773 -v MA.indi.chr12-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f8
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_12_v3:1285068-1688600 -v MA.indi.chr12-3.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f9
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_12_v3:1745531-2163700 -v MA.indi.chr12-4.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f10
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_13_v3:74414-1168127 -v MA.indi.chr13-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f11
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_13_v3:1170426-2791900 -v MA.indi.chr13-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f12
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_14_v3:35775-1071523 -v MA.indi.chr14-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f13
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.indi.txt -r Pf3D7_14_v3:1075090-3255710 -v MA.indi.chr14-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20

f14
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_11_v3:110001-831968 -v MA.merge.chr11-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20

f15
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_11_v3:834246-2003320 -v MA.merge.chr11-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20

f16
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_12_v3:60301-766654 -v MA.merge.chr12-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20

f17
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_12_v3:780451-1282773 -v MA.merge.chr12-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f18
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_12_v3:1285068-1688600 -v MA.merge.chr12-3.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f19
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_12_v3:1745531-2163700 -v MA.merge.chr12-4.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f20
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_13_v3:74414-1168127 -v MA.merge.chr13-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f21
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_13_v3:1170426-2791900 -v MA.merge.chr13-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f22
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_14_v3:35775-1071523 -v MA.merge.chr14-1.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20
f23
/master/xli/software/FreeBayes/freebayes/bin/freebayes -f /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -L bam.merge.txt -r Pf3D7_14_v3:1075090-3255710 -v MA.merge.chr14-2.freebayes.vcf -p 1 --min-repeat-size 3 --no-partial-observations -m 30 -q 20


### cat vcf
java -cp /master/xli/software/GATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V MA.indi.chr1-chr8.freebayes.vcf \
-V MA.indi.chr9-1.freebayes.vcf \
-V MA.indi.chr9-2.freebayes.vcf \
-V MA.indi.chr10.freebayes.vcf \
-V MA.indi.chr11-1.freebayes.vcf \
-V MA.indi.chr11-2.freebayes.vcf \
-V MA.indi.chr12-1.freebayes.vcf \
-V MA.indi.chr12-2.freebayes.vcf \
-V MA.indi.chr12-3.freebayes.vcf \
-V MA.indi.chr12-4.freebayes.vcf \
-V MA.indi.chr13-1.freebayes.vcf \
-V MA.indi.chr13-2.freebayes.vcf \
-V MA.indi.chr14-1.freebayes.vcf \
-V MA.indi.chr14-2.freebayes.vcf \
-out MA.indi.vcf \
-assumeSorted

java -cp /master/xli/software/GATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V MA.merge.chr1-chr10.freebayes.vcf \
-V MA.merge.chr11-1.freebayes.vcf \
-V MA.merge.chr11-2.freebayes.vcf \
-V MA.merge.chr12-1.freebayes.vcf \
-V MA.merge.chr12-2.freebayes.vcf \
-V MA.merge.chr12-3.freebayes.vcf \
-V MA.merge.chr12-4.freebayes.vcf \
-V MA.merge.chr13-1.freebayes.vcf \
-V MA.merge.chr13-2.freebayes.vcf \
-V MA.merge.chr14-1.freebayes.vcf \
-V MA.merge.chr14-2.freebayes.vcf \
-out MA.merge.vcf \
-assumeSorted

vcfallelicprimitives -kg MA.indi.vcf | /master/xli/software/vt/vt normalize -r /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -o MA.indi.norm.vcf -
vcfallelicprimitives -kg MA.merge.vcf | /master/xli/software/vt/vt normalize -r /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta -o MA.merge.norm.vcf -

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant MA.merge.norm.vcf \
-o MA.merge.norm.rm50-59.vcf \
-xl_sn  merge59/ \
-xl_sn merge50/ 



java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant MA.merge.norm.rm50-59.vcf \
-o MA.merge.norm.rm50-59.rm-non-variant.vcf \
-env



java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant MA.indi.vcf \
-o MA.indi.rm50-59.vcf \
-xl_sn MA50/ \
-xl_sn MA59/ \
-xl_sn Sample_50/ \
-xl_sn Sample_59/ 



java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant MA.indi.rm50-59.vcf \
-o MA.indi.rm50-59.rm-non-variant.vcf \
-env


java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant RUN_seperate_haploid_filter.vcf \
-o MA.HipSTR.basicfilter.vcf \
-env





