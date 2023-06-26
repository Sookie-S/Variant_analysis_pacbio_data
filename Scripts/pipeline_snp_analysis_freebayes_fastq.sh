#!/bin/sh

###################################################################
#Script Name	: pipeline.sh                                                                                                                                                                              
#Author       	: Soukaina Timouma                                             
#Email         	: soukaina.timouma@manchester.ac.uk                                           
###################################################################


# conda install -c bioconda pbmm2


for ref_genome in 'NRRL' 'CBS8639' 'CBS8638'
do
  for reads in 'NRRLreads' 'CBS8639reads' 'CBS8638reads'
  do 
	echo "\n######### $ref_genome vs $reads #########\n"

  echo "-- Reference genome index"  
  pbmm2 index quickstart-testdata/$ref_genome.fasta quickstart-testdata/$ref_genome.mmi
  samtools faidx quickstart-testdata/$ref_genome.fasta

  echo "-- Reads alignment to reference genome" 
  pbmm2 align quickstart-testdata/$ref_genome.mmi quickstart-testdata/"$reads"_CCS.fastq quickstart-testdata/"$ref_genome"_vs_"$reads"_fastq.bam --preset CCS --sort --rg '@RG\tID:myid\tSM:mysample'

  echo "-- Bam file index"  
  samtools index quickstart-testdata/"$ref_genome"_vs_"$reads"_fastq.bam

	echo "-- calling variants with freebayes"	
	freebayes -f quickstart-testdata/"$ref_genome".fasta quickstart-testdata/"$ref_genome"_vs_"$reads"_fastq.bam > quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes.vcf

  echo "-- Converting VCF to compressed BCF"
  
  bcftools index quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes.vcf
  bcftools view -Ob quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes.vcf > quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes.bcf
 
  echo "-- Number of SNPS:"
  bcftools view -v snps quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes.vcf | grep -v -c ‘^#’ > quickstart-output/freebayes_res_fastq/All/"$ref_genome"_vs_"$reads"_fastq_freebayes_snps.txt
  echo "-- Number of INDELS:"
  bcftools view -v indels quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes.vcf | grep -v -c ‘^#’ > quickstart-output/freebayes_res_fastq/All/"$ref_genome"_vs_"$reads"_fastq_freebayes_indels.txt

  echo "-- Variant calling index"
  bcftools index quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes.bcf

  echo "-- Extract variants in CDS"
  bcftools view quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes.bcf -R quickstart-testdata/"$ref_genome"_CDS.bed > quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes_CDS.vcf
  
  echo "-- Number of SNPS:"
  bcftools view -v snps quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes_CDS.vcf | grep -v -c ‘^#’ > quickstart-output/freebayes_res_fastq/CDS/"$ref_genome"_vs_"$reads"_fastq_freebayes_CDS_snps.txt
  echo "-- Number of INDELS:"
  bcftools view -v indels quickstart-output/freebayes_res_fastq/"$ref_genome"_vs_"$reads"_fastq_freebayes_CDS.vcf | grep -v -c ‘^#’ > quickstart-output/freebayes_res_fastq/CDS/"$ref_genome"_vs_"$reads"_fastq_freebayes_CDS_indels.txt

done
done