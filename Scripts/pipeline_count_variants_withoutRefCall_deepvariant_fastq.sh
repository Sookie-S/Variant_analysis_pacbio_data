#!/bin/sh

###################################################################
#Script Name	: pipeline.sh                                                                                                                                                                              
#Author       	: Soukaina Timouma                                             
#Email         	: soukaina.timouma@manchester.ac.uk                                           
###################################################################


echo "\n\n---------- Reads alignment to reference analysis"


for ref_genome in 'NRRL' 'CBS8639' 'CBS8638'
do
  for reads in 'NRRLreads' 'CBS8639reads' 'CBS8638reads'
  do 
	echo "\n######### $ref_genome vs $reads #########\n"

	echo "-- Reference genome index"	
	pbmm2 index quickstart-testdata/$ref_genome.fasta quickstart-testdata/$ref_genome.mmi
	samtools faidx quickstart-testdata/$ref_genome.fasta

	echo "-- Reads alignment to reference genome"	
	pbmm2 align quickstart-testdata/$ref_genome.mmi quickstart-testdata/"$reads"_CCS.fastq quickstart-testdata/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.bam --preset CCS --sort --rg '@RG\tID:myid\tSM:mysample'

	echo "-- Bam file index"	
	samtools index quickstart-testdata/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.bam
  done
done

echo "\n\n---------- Variant Calling analysis"

for ref_genome in 'NRRL' 'CBS8639' 'CBS8638'
do
  for reads in 'NRRLreads' 'CBS8639reads' 'CBS8638reads'
  do 
	echo "\n######### $ref_genome vs $reads #########\n"


	echo "-- Converting VCF to compressed BCF"
	grep -v "RefCall" quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq.vcf > quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.vcf
	bcftools index quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.vcf
	bcftools view -Ob quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.vcf > quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.bcf

	echo "-- Number of SNPS:"
	bcftools view -v snps quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.vcf | grep -v -c ‘^#’ > quickstart-output/deepvariant_res/All/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall_snps.txt
	echo "-- Number of INDELS:"
	bcftools view -v indels quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.vcf | grep -v -c ‘^#’ > quickstart-output/deepvariant_res/All/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall_indels.txt

	echo "-- Variant calling index"
	bcftools index quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.bcf

	echo "-- Extract variants in CDS"
	bcftools view quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall.bcf -R quickstart-testdata/"$ref_genome"_CDS.bed > quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall_CDS.vcf
	
	echo "-- Number of SNPS:"
	bcftools view -v snps quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall_CDS.vcf | grep -v -c ‘^#’ > quickstart-output/deepvariant_res/CDS/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall_CDS_snps.txt
	echo "-- Number of INDELS:"
	bcftools view -v indels quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall_CDS.vcf | grep -v -c ‘^#’ > quickstart-output/deepvariant_res/CDS/"$ref_genome"_vs_"$reads"_fastq_withoutRefCall_CDS_indels.txt

  done
done
