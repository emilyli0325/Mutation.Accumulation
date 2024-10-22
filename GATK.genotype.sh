#### genotype ####
# run1
MA12.g.vcf  MA21.g.vcf  MA26.g.vcf  MA32.g.vcf  MA38.g.vcf   MA41.g.vcf  MA45.g.vcf  MA50.g.vcf  MA54.g.vcf  MA58.g.vcf  MA62.g.vcf
MA18.g.vcf  MA22.g.vcf  MA27.g.vcf  MA35.g.vcf  MA39.g.vcf   MA42.g.vcf  MA46.g.vcf  MA51.g.vcf  MA55.g.vcf  MA59.g.vcf  MA6.g.vcf
MA19.g.vcf  MA23.g.vcf  MA29.g.vcf  MA36.g.vcf  MA3D7.g.vcf  MA43.g.vcf  MA47.g.vcf  MA52.g.vcf  MA56.g.vcf  MA60.g.vcf
MA20.g.vcf  MA24.g.vcf  MA30.g.vcf  MA37.g.vcf  MA40.g.vcf   MA44.g.vcf  MA49.g.vcf  MA53.g.vcf  MA57.g.vcf  MA61.g.vcf

# run2
Sample_12.g.vcf  Sample_23.g.vcf  Sample_32.g.vcf  Sample_3D7.g.vcf  Sample_45.g.vcf  Sample_52.g.vcf  Sample_58.g.vcf
Sample_18.g.vcf  Sample_24.g.vcf  Sample_35.g.vcf  Sample_40.g.vcf   Sample_46.g.vcf  Sample_53.g.vcf  Sample_59.g.vcf
Sample_19.g.vcf  Sample_26.g.vcf  Sample_36.g.vcf  Sample_41.g.vcf   Sample_47.g.vcf  Sample_54.g.vcf  Sample_60.g.vcf
Sample_20.g.vcf  Sample_27.g.vcf  Sample_37.g.vcf  Sample_42.g.vcf   Sample_49.g.vcf  Sample_55.g.vcf  Sample_61.g.vcf
Sample_21.g.vcf  Sample_29.g.vcf  Sample_38.g.vcf  Sample_43.g.vcf   Sample_50.g.vcf  Sample_56.g.vcf  Sample_62.g.vcf
Sample_22.g.vcf  Sample_30.g.vcf  Sample_39.g.vcf  Sample_44.g.vcf   Sample_51.g.vcf  Sample_57.g.vcf  Sample_6.g.vcf

# run3
3D7.g.vcf

# merge
merge12.g.vcf  merge22.g.vcf  merge29.g.vcf  merge37.g.vcf   merge41.g.vcf  merge46.g.vcf  merge52.g.vcf  merge57.g.vcf  merge62.g.vcf
merge18.g.vcf  merge23.g.vcf  merge30.g.vcf  merge38.g.vcf   merge42.g.vcf  merge47.g.vcf  merge53.g.vcf  merge58.g.vcf  merge6.g.vcf
merge19.g.vcf  merge24.g.vcf  merge32.g.vcf  merge39.g.vcf   merge43.g.vcf  merge49.g.vcf  merge54.g.vcf  merge59.g.vcf
merge20.g.vcf  merge26.g.vcf  merge35.g.vcf  merge3D7.g.vcf  merge44.g.vcf  merge50.g.vcf  merge55.g.vcf  merge60.g.vcf
merge21.g.vcf  merge27.g.vcf  merge36.g.vcf  merge40.g.vcf   merge45.g.vcf  merge51.g.vcf  merge56.g.vcf  merge61.g.vcf

# merge genotype
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs -R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V merge12.g.vcf -V merge22.g.vcf -V merge29.g.vcf -V merge37.g.vcf  -V merge41.g.vcf -V merge46.g.vcf -V merge52.g.vcf -V merge57.g.vcf -V merge62.g.vcf -V merge18.g.vcf -V merge23.g.vcf -V merge30.g.vcf -V merge38.g.vcf  -V merge42.g.vcf -V merge47.g.vcf -V merge53.g.vcf -V merge58.g.vcf -V merge6.g.vcf -V merge19.g.vcf -V merge24.g.vcf -V merge32.g.vcf -V merge39.g.vcf  -V merge43.g.vcf -V merge49.g.vcf -V merge54.g.vcf -V merge20.g.vcf -V merge26.g.vcf -V merge35.g.vcf -V merge3D7.g.vcf -V merge44.g.vcf -V merge55.g.vcf -V merge60.g.vcf -V merge21.g.vcf -V merge27.g.vcf -V merge36.g.vcf -V merge40.g.vcf -V merge45.g.vcf -V merge51.g.vcf -V merge56.g.vcf -V merge61.g.vcf \
--useNewAFCalculator \
--sample_ploidy 1 \
-o merge.GATK.vcf

# individual genotype
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs -R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V MA12.g.vcf -V MA21.g.vcf -V MA26.g.vcf -V MA32.g.vcf -V MA38.g.vcf -V MA41.g.vcf -V MA45.g.vcf -V MA54.g.vcf -V MA58.g.vcf -V MA62.g.vcf -V MA18.g.vcf -V MA22.g.vcf -V MA27.g.vcf -V MA35.g.vcf -V MA39.g.vcf -V MA42.g.vcf -V MA46.g.vcf -V MA51.g.vcf -V MA55.g.vcf -V MA6.g.vcf -V MA19.g.vcf -V MA23.g.vcf -V MA29.g.vcf -V MA36.g.vcf -V MA3D7.g.vcf -V MA43.g.vcf -V MA47.g.vcf -V MA52.g.vcf -V MA56.g.vcf -V MA60.g.vcf -V MA20.g.vcf -V MA24.g.vcf -V MA30.g.vcf -V MA37.g.vcf -V MA40.g.vcf -V MA44.g.vcf -V MA49.g.vcf -V MA53.g.vcf -V MA57.g.vcf -V MA61.g.vcf -V Sample_12.g.vcf -V Sample_23.g.vcf -V Sample_32.g.vcf -V Sample_3D7.g.vcf -V Sample_45.g.vcf -V Sample_52.g.vcf -V Sample_58.g.vcf -V Sample_18.g.vcf -V Sample_24.g.vcf -V Sample_35.g.vcf -V Sample_40.g.vcf -V Sample_46.g.vcf -V Sample_53.g.vcf -V Sample_19.g.vcf -V Sample_26.g.vcf -V Sample_36.g.vcf -V Sample_41.g.vcf -V Sample_47.g.vcf -V Sample_54.g.vcf -V Sample_60.g.vcf -V Sample_20.g.vcf -V Sample_27.g.vcf -V Sample_37.g.vcf -V Sample_42.g.vcf -V Sample_49.g.vcf -V Sample_55.g.vcf -V Sample_61.g.vcf -V Sample_21.g.vcf -V Sample_29.g.vcf -V Sample_38.g.vcf -V Sample_43.g.vcf -V Sample_56.g.vcf -V Sample_62.g.vcf -V Sample_22.g.vcf -V Sample_30.g.vcf -V Sample_39.g.vcf -V Sample_44.g.vcf -V Sample_51.g.vcf -V Sample_57.g.vcf -V Sample_6.g.vcf -V 3D7.g.vcf \
--useNewAFCalculator \
--sample_ploidy 1 \
-o individual.GATK.vcf


# Seperate SNP and Indel
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V individual.GATK.vcf \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
-selectType SNP -o individual.SNP.vcf
	
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V individual.GATK.vcf \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
-selectType INDEL -o individual.INDEL.vcf

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V merge.GATK.vcf \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
-selectType SNP -o merge.SNP.vcf
	
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V merge.GATK.vcf \
-L /master/xli/Index/Known_sites/Core_Genome.intervals \
-selectType INDEL -o merge.INDEL.vcf


# basic filter
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant individual.SNP.vcf \
-o individual.SNP_basic_filter.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
--filterName "BasicFilter" 

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant merge.SNP.vcf \
-o merge.SNP_basic_filter.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
--filterName "BasicFilter" 

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant individual.INDEL.vcf \
-o individual.INDEL_basic_filter.vcf \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
--filterName "BasicFilter"

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
--variant merge.INDEL.vcf \
-o merge.INDEL_basic_filter.vcf \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
--filterName "BasicFilter"


individual.SNP_basic_filter.vcf
individual.INDEL_basic_filter.vcf

merge.SNP_basic_filter.vcf
merge.INDEL_basic_filter.vcf 


java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
-T VariantsToTable \
-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
-V merge.INDEL_basic_filter.vcf \
-F CHROM -F REF -F ALT -F POS -F QUAL -F FILTER -GF GT -GF AD -GF DP \
-o merge.INDEL.table
VariantAnnotator -list

AC=1;AF=0.025;AN=40;BaseQRankSum=0.516;ClippingRankSum=0.00;DP=1934;FS=4.590;MLEAC=1;MLEAF=0.025;MQ=73.76;MQRankSum=0.605;QD=1.63;ReadPosRankSum=0.385;SOR=0.985
AC=1;AF=0.027;AN=37;BaseQRankSum=-9.700e-01;ClippingRankSum=0.00;DP=643;FS=4.012;MLEAC=1;MLEAF=0.027;MQ=81.59;MQRankSum=0.00;QD=9.24;ReadPosRankSum=-3.190e-01;SOR=0.392


AC=1;AF=0.048;AN=21;BaseQRankSum=-3.145e+00;ClippingRankSum=0.00;DP=203;FS=9.551;MLEAC=2;MLEAF=0.095;MQ=39.02;MQRankSum=1.32;QD=7.88;ReadPosRankSum=1.32;SOR=0.970
AC=1;AF=0.040;AN=25;DP=171;FS=0.000;MLEAC=1;MLEAF=0.040;MQ=48.39;QD=17.83;SOR=0.693
AC=7;AF=0.304;AN=23;BaseQRankSum=0.00;ClippingRankSum=0.00;DP=315;FS=5.151;MLEAC=11;MLEAF=0.478;MQ=52.87;MQRankSum=0.431;QD=6.09;ReadPosRankSum=1.02;SOR=0.330

QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0
AC=2;AF=0.069;AN=29;BaseQRankSum=-1.068e+00;ClippingRankSum=0.00;DP=538;FS=160.205;MLEAC=3;MLEAF=0.103;MQ=122.81;MQRankSum=0.00;QD=9.19;ReadPosRankSum=0.122;SOR=1.489

#### VQSR #change parameters according to real dataset

### WARNING: Training with very few variant sites!
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
	-T VariantRecalibrator \
	-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
	-input MA.GATK.basicfilter.SNP.vcf \
	--resource:dbsnp,known=true,training=true,truth=true,prior=15.0 /master/xli/Index/Known_sites/3d7_hb3.combined.final.karo.sort.vcf \
	--resource:dbsnp,known=true,training=true,truth=true,prior=15.0 /master/xli/Index/Known_sites/7g8_gb4.combined.final.karo.sort.vcf \
	--resource:dbsnp,known=true,training=true,truth=true,prior=15.0 /master/xli/Index/Known_sites/hb3_dd2.combined.final.karo.sort.vcf \
	-an QD -an FS -an SOR -an ReadPosRankSum -an MQ \
	-mode SNP \
	-tranche 100.0 -tranche 90 -tranche 80 -tranche 70 \
	-recalFile recalibrate_SNP.recal \
	-tranchesFile recalibrate_SNP.tranches \
	-rscriptFile recalibrate_SNP_plots.R
	

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
    -input round2.vcf \
    -mode SNP \
    --ts_filter_level 99 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -o SNP_VQSR.vcf
	

java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
	-T VariantRecalibrator \
	-R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
    -input MA.GATK.basicfilter.vcf \
	--resource:dbsnp,known=true,training=true,truth=true,prior=15.0 /master/xli/Index/Known_sites/3d7_hb3.combined.final.karo.sort.vcf \
	--resource:dbsnp,known=true,training=true,truth=true,prior=15.0 /master/xli/Index/Known_sites/7g8_gb4.combined.final.karo.sort.vcf \
	--resource:dbsnp,known=true,training=true,truth=true,prior=15.0 /master/xli/Index/Known_sites/hb3_dd2.combined.final.karo.sort.vcf \
	-an QD -an FS -an SOR -an ReadPosRankSum -an MQ \
	-mode INDEL \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    --maxGaussians 4 \
    -recalFile recalibrate_INDEL.recal \
    -tranchesFile recalibrate_INDEL.tranches \
    -rscriptFile recalibrate_INDEL_plots.R
	
java -jar /master/xli/software/GATK/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /master/xli/Index/Pfal32_GATK_index/PlasmoDB-32_Pfalciparum3D7_Genome.fasta \
    -input MA.GATK.basicfilter.vcf \
    -mode INDEL \
    --ts_filter_level 99 \
    -recalFile recalibrate_INDEL.recal \
    -tranchesFile recalibrate_INDEL.tranches \
    -o MA.GATK.basicfilter.Indel_VQSR.vcf 
	
	