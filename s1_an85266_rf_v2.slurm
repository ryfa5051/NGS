#!/bin/bash
#SBATCH --job-name=s1_an85266_rf
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ryfa5051@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=23:59:00
#SBATCH -p short
#SBATCH --mem=500gb
#SBATCH --output=/Users/ryfa5051/NGS/eo/%j.out
#SBATCH --error=/Users/ryfa5051/NGS/eo/%j.e
#
#directory designation 1.0
sample=an85266_mapped
home=/Users/ryfa5051/NGS
input1=/Users/ryfa5051/NGS/fastq/batch2/180625_M00517_0482_fastq/Sawyer-537-1_S1_L001_R1_001.fastq.gz
input2=/Users/ryfa5051/NGS/fastq/batch2/180625_M00517_0482_fastq/Sawyer-537-1_S1_L001_R2_001.fastq.gz
adapter=$home/trim/adapters.fa
#
#Step 1: Trimming adaptor and low quality reads (note lowered window qual from 4:15 -> 4:10, and MINLEN from 36->30)
module load trimmomatic
java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar PE $input1 $input2 $home/fastq/trimmed/R1.an85266.p.fastq $home/fastq/trimmed/R1.an85266.up.fastq $home/fastq/trimmed/R2.an85266.p.fastq $home/fastq/trimmed/R2.dsamp.an85266.up.fastq ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
#directory designation 2.0
fullfrag_ref=/Users/ryfa5051/NGS/bwa_sam_trial/OM_AN_FRAG_REF_FIN.fa
CDSref=/Users/ryfa5051/NGS/bwa_sam_trial/OM_An_5_genes.fa
pair1_fastq1=$home/fastq/trimmed/R1.an85266.p.fastq
pair2_fastq2=$home/fastq/trimmed/R2.an85266.p.fastq
#
#Step 2: bwa indexing and mapping ((NOTE-take a 2nd look at -R @RG read groups readout to make sure lane designation is ok))
module load bwa
bwa index $fullfrag_ref
#for some reason the output (>) below will not accept any designated strings ($sample for ex.) - full path is required otherwise it outputs as '.sam'
bwa mem -M -R '@RG\tID:batch1\tSM:an85266\tPL:illumina\tLB:lib1\tPU:nano1' $fullfrag_ref $pair1_fastq1 $pair2_fastq2 > /Users/ryfa5051/NGS/bwa_sam_trial/an85266_mapped.sam
#directory designation 3.0 
aligned_an85266=/Users/ryfa5051/NGS/bwa_sam_trial/an85266_mapped.sam
#
#***** START OF samtoolstrialcode.slurm *****
#
#indexing full frag ref seq
module load samtools
#samtools faidx $fullfrag_ref
#sam2bam
samtools view $aligned_an85266 -T $fullfrag_ref -b -o $home/bwa_sam_trial/an85266_mapped.bam
#sort bam file and then index
samtools sort $home/bwa_sam_trial/an85266_mapped.bam -o $home/bwa_sam_trial/an85266_mapped_sorted.bam
samtools index $home/bwa_sam_trial/an85266_mapped_sorted.bam
#pileup to call variants
module load bcftools
bcftools mpileup -Ou -f $fullfrag_ref $home/bwa_sam_trial/an85266_mapped_sorted.bam -o $home/bwa_sam_trial/an85266_mapped_sorted.bcf
#calling and filtering
bcftools call -Ou -m -v $home/bwa_sam_trial/an85266_mapped_sorted.bcf -o $home/bwa_sam_trial/an85266_mapped_sorted_called.bcf
bcftools filter -Oz -e '%QUAL<15' $home/bwa_sam_trial/an85266_mapped_sorted_called.bcf -o $home/bwa_sam_trial/an85266_mapped_sorted_called_flt.vcf.gz
#REST IS VCF 2 FASTA 
module load gatk
module load picard
#index filtered .vcf.gz
bcftools tabix $home/bwa_sam_trial/an85266_mapped_sorted_called_flt.vcf.gz
#
#picard to assign sample name RG
java -jar /opt/picard/2.6.0/picard-2.6.0.jar RenameSampleInVcf \
	INPUT=$home/bwa_sam_trial/an85266_mapped_sorted_called_flt.vcf.gz \
	OUTPUT=$home/bwa_sam_trial/fasta/picard_renamed_an85266.vcf.gz \
	NEW_SAMPLE_NAME=an85266
bcftools tabix $home/bwa_sam_trial/fasta/picard_renamed_an85266.vcf.gz
#making seq directory for full frag file
java -jar /opt/picard/2.6.0/picard-2.6.0.jar CreateSequenceDictionary \
	R=$home/bwa_sam_trial/OM_AN_FRAG_REF_FIN.fa \
	O=$home/bwa_sam_trial/OM_AN_FRAG_REF_FIN.dict
samtools faidx $home/bwa_sam_trial/OM_AN_FRAG_REF_FIN.fa
#gatk FastaAlternateReferenceMaker w/ ambiguity codes
java -jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar \
 	--analysis_type FastaAlternateReferenceMaker \
 	--reference_sequence /Users/ryfa5051/NGS/bwa_sam_trial/OM_AN_FRAG_REF_FIN.fa \
 	--out $home/bwa_sam_trial/fasta/fin/an85266_final_fin.fasta \
 	--variant $home/bwa_sam_trial/fasta/picard_renamed_an85266.vcf.gz \
 	-IUPAC an85266
