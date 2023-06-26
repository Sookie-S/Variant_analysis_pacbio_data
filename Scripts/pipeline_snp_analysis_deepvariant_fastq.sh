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
	pbmm2 align quickstart-testdata/$ref_genome.mmi quickstart-testdata/"$reads"_CCS.fastq quickstart-testdata/"$ref_genome"_vs_"$reads"_fastq.bam --preset CCS --sort --rg '@RG\tID:myid\tSM:mysample'

	echo "-- Bam file index"	
	samtools index quickstart-testdata/"$ref_genome"_vs_"$reads"_fastq.bam
  done
done

echo "\n\n---------- Variant Calling analysis"

for ref_genome in 'NRRL' 'CBS8639' 'CBS8638'
do
  for reads in 'NRRLreads' 'CBS8639reads' 'CBS8638reads'
  do 
	echo "\n######### $ref_genome vs $reads #########\n"


	echo "-- DeepVariant analysis"

	INPUT_DIR="${PWD}/quickstart-testdata"
	OUTPUT_DIR="${PWD}/quickstart-output/deepvariant_res"

	BIN_VERSION="1.5.0"
	sudo docker run \
	  -v "${INPUT_DIR}":"/input" \
    -v "${OUTPUT_DIR}":"/output" \
	  google/deepvariant:"${BIN_VERSION}" \
	  /opt/deepvariant/bin/run_deepvariant \
	  --model_type=PACBIO \
	  --ref=/input/"$ref_genome".fasta \
	  --reads=/input/"$ref_genome"_vs_"$reads"_fastq.bam \
	  --output_vcf=/output/"$ref_genome"_vs_"$reads"_fastq.vcf \
	  --output_gvcf=/output/"$ref_genome"_vs_"$reads"_fastq.gvcf \
	  --num_shards=$(nproc) \
	  --logging_dir=/output/logs \
	  --dry_run=false
	  
	echo "-- Converting VCF to compressed BCF"
	
	bcftools index quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq.vcf
	bcftools view -Ob quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq.vcf > quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq.bcf

	echo "-- Number of SNPS:"
	bcftools view -v snps quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq.vcf | grep -v -c ‘^#’ > quickstart-output/deepvariant_res/All/"$ref_genome"_vs_"$reads"_fastq_snps.txt
	echo "-- Number of INDELS:"
	bcftools view -v indels quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq.vcf | grep -v -c ‘^#’ > quickstart-output/deepvariant_res/All/"$ref_genome"_vs_"$reads"_fastq_indels.txt



	echo "-- Variant calling index"
	bcftools index quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq.bcf

	echo "-- Extract variants in CDS"
	bcftools view quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq.bcf -R quickstart-testdata/"$ref_genome"_CDS.bed > quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_CDS.vcf
	
	echo "-- Number of SNPS:"
	bcftools view -v snps quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_CDS.vcf | grep -v -c ‘^#’ > quickstart-output/deepvariant_res/CDS/"$ref_genome"_vs_"$reads"_fastq_CDS_snps.txt
	echo "-- Number of INDELS:"
	bcftools view -v indels quickstart-output/deepvariant_res/"$ref_genome"_vs_"$reads"_fastq_CDS.vcf | grep -v -c ‘^#’ > quickstart-output/deepvariant_res/CDS/"$ref_genome"_vs_"$reads"_fastq_CDS_indels.txt

  done
done
