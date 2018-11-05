#!/bin/bash
#SBATCH --job-name=3_batch_Miseq
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ryfa5051@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=23:59:00
#SBATCH -p short
#SBATCH --mem=500gb
#SBATCH --output=/Users/ryfa5051/NGS/eo/%x,%j.out
#SBATCH --error=/Users/ryfa5051/NGS/eo/%x,%j.e
#MODULE LOAD
start=$(date +%s)
echo "NGS Processing Start"
echo $(date)
module load trimmomatic bwa samtools bcftools vcftools picard gatk fastx-toolkit/0.0.13 bbmap fastqc
###############################REPLACE_BatchNum_BELOW###############################
batchNum=3_batch
batchpath=/Users/ryfa5051/NGS/$batchNum
mkdir /Users/ryfa5051/NGS/fastq/$batchNum
mkdir $batchpath
#ref seq cp & var:
cp /Users/ryfa5051/NGS/reference/ref_frags_3batch.fa $batchpath
ref_seq=/Users/ryfa5051/NGS/$batchNum/ref_frags_3batch.fa
## potential utilizing dup ref for A3G filtering fasta variable is: ref_seqYESdup6genes=/Users/ryfa5051/NGS/3_batch/duptest_reffrags.fa
#BATCH CORE_FASTQ DIR:
indir=/Users/ryfa5051/NGS/fastq/$batchNum/core_fastq
mkdir $indir
#potentially add rsync step here so manually making the core_fastq directory for each new run isnt nessesary?
#I.E. ftp://http/core_serv7/biof/SEQ/Sawyer/miseq/*REPLACE_WITH_CORE_BATCH_ID* $indir
#DIRECTORY DESIGNATIONS:
outdirTRIM=/Users/ryfa5051/NGS/fastq/$batchNum/trimmed
mkdir $outdirTRIM
indir_sams=$batchpath
out_bams=$batchpath/bams
mkdir $out_bams
metrics=$batchpath/bams/dup_metrics
mkdir $metrics
out_vcfs=$batchpath/vcfs
mkdir $out_vcfs
out_fasta=$batchpath/fastas
mkdir $out_fasta
final_fasta=$out_fasta/final_fasta
mkdir $final_fasta
fastQCdir=/Users/ryfa5051/NGS/fastqc/$batchNum
mkdir $fastQCdir
fastQC_before=$fastQCdir/before
fastQC_after=$fastQCdir/after
mkdir $fastQC_before
mkdir $fastQC_after
#TRIMMING WITH TRIMMOMATIC
echo "*** TRIMMING START" $(date) "***"
for pathandfilename in $(ls $indir/*R1_001.fastq.gz); do
	filename1=$(basename $pathandfilename)
	rootfilename=$(basename $pathandfilename R1_001.fastq.gz)
	filename2=${rootfilename}R2_001.fastq.gz
echo $filename1
echo $filename2
echo $pathandfilename
java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar PE -quiet $indir/$filename1 $indir/$filename2 $outdirTRIM/${filename1//R1_001.fastq.gz/R1.PE.fastq} $outdirTRIM/${filename1//R1_001.fastq.gz/R1.UP.fastq} $outdirTRIM/${filename2//R2_001.fastq.gz/R2.PE.fastq} $outdirTRIM/${filename2//R2_001.fastq.gz/R2.UP.fastq} ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/NexteraPE-PE.fa:2:30:10:1:"false" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
done
#trimmomatic future note: -basename option eliminates the need for 4 file designation
echo "*** TRIMMING COMPLETE *** MAPPING START ***" $(date) "***"
#FASTQC AFTER TRIMMING & REMOVING .ZIPS - .HTMLS ONLY
fastqc -q -t 1728 $indir/*R1_001.fastq.gz --outdir=$fastQC_before/
fastqc -q -t 1728 $indir/*R2_001.fastq.gz  --outdir=$fastQC_before/
fastqc -q -t 1728 $outdirTRIM/*.PE.fastq -o $fastQC_after/
rm $fastQC_before/*.zip
rm $fastQC_after/*.zip
#MAPPING WITH BWA & INDEXING REFERENCE
cd $batchpath
bwa index $ref_seq
samtools faidx $ref_seq
samtools dict -s Nancymaae $ref_seq 
cd ~
#MAPPING WITH BWA MEM
for pathandtrimfilename in $(ls $outdirTRIM/*R1.PE.fastq); do
	mapname1=$(basename $pathandtrimfilename)
	rootmapname=$(basename $pathandtrimfilename R1.PE.fastq)
	mapname2=${rootmapname}R2.PE.fastq
echo $mapname1
echo $mapname2
bwa mem -M -R "@RG\tID:batch2\tSM:${mapname1//_L001_R1.PE.fastq/_}\tPL:illumina\tLB:lib1\tPU:nano1" $ref_seq $outdirTRIM/$mapname1 $outdirTRIM/$mapname2 > $batchpath/${mapname1//_L001_R1.PE.fastq/_mapped.sam}
done
echo $(date), "MAPPING COMPLETE"
#SAM2BAM AND SORTING
for pathandsamfilename in $(ls $indir_sams/*.sam); do
	samname=$(basename $pathandsamfilename)
echo $pathandsamfilename
#PICARD SORTSAM
	java -jar /opt/picard/2.6.0/picard-2.6.0.jar SortSam \
		INPUT=$pathandsamfilename \
		OUTPUT=$out_bams/${samname//.sam/.bam} \
		SORT_ORDER=coordinate 
done
#GATK MARKDUPLICATES
for pathandbamfilename in $(ls $out_bams/*_mapped.bam); do
	basebamname=$(basename $pathandbamfilename)
	java -jar /opt/picard/2.6.0/picard-2.6.0.jar MarkDuplicates \
		I=$pathandbamfilename \
		O=${pathandbamfilename//_mapped.bam/_mapped_mdup.bam} \
		M=$metrics/${basebamname//_mapped.bam/.txt} \
		CREATE_INDEX=true \
		REMOVE_DUPLICATES=true
done
#REMOVE NON-MARKED BAMS
rm --force $out_bams/*_mapped.bam
#holstein phasing algo goes here for future if using non vcf-reliant/trio version
for mapsortbams in $(ls $out_bams/*_mdup.bam); do
	bambasename=$(basename $mapsortbams)
echo $mapsortbams
echo $bambasename
	bcftools mpileup -Ou -d 8000 --adjust-MQ 50 -a AD,ADF,ADR,DP,SP,DP4 -f $ref_seq $mapsortbams | \
	bcftools call -Ou -m -v -P 1.1e-2 - - | \
	bcftools filter -Ov -e '%QUAL<20' - -o $out_vcfs/${bambasename//.bam/.vcf}
done
#adding additional INFO/FORMAT fields for VCF files via bcftools plugins (AF,MAF,HWE etc)
export BCFTOOLS_PLUGINS=/opt/bcftools/1.8/libexec/bcftools
for called_vcfs in $(ls $out_vcfs/*_mdup.vcf); do
	bcftools +fill-tags $called_vcfs -o ${called_vcfs//_mdup.vcf/_mdup.vcf}
done
#REF SEQ DICTIONARY FILE - PICARD
java -jar /opt/picard/2.6.0/picard-2.6.0.jar CreateSequenceDictionary \
	R=$ref_seq \
	O=${ref_seq//.fa/.dict}
#Making Fastas with ambiguity codes
for vcfs in $(ls $out_vcfs/*_mdup.vcf); do
	idSM="$(bcftools query -l $vcfs)"
	vcfbasename=$(basename $vcfs)
java -jar /opt/gatk/3.7.0/GenomeAnalysisTK.jar \
 	--analysis_type FastaAlternateReferenceMaker \
 	--reference_sequence $ref_seq \
 	--out $out_fasta/${vcfbasename//.vcf/.fasta} \
 	--variant $vcfs \
 	-IUPAC "${idSM}"
 done
#FASTA MANIPULATION - RENAMING AND MERGING LIKE FILES BY GENE LABELED BY INDIVIDUAL
for fastas in $(ls $out_fasta/*mdup.fasta); do
	fastabasename=$(basename $fastas)
	rootfilename=$(basename $fastabasename .fasta)
	vcfbasename=${rootfilename}.vcf
	idSM="$(bcftools query -l $out_vcfs/$vcfbasename)"
	sed -i -e "s/1 /${idSM}/;s/2 /${idSM}/;s/3 /${idSM}/;s/4 /${idSM}/;s/5 /${idSM}/;s/:1//g;s/-//" $fastas
	fasta_formatter -w 0 -i $fastas -o $out_fasta/${fastabasename//.fasta/_fin.fasta}
	awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}' $out_fasta/${fastabasename//.fasta/_fin.fasta}
done
#CDS trim each gene from full pcr fragment
cat ~/*TRIMCYP_Nancymaae.fasta | bbduk.sh in=stdin.fasta out=$final_fasta/"${batchNum}"_TRIMCYP_Nancymaae.fasta ftl=62 ftr=1486
cat ~/*TETHERIN_Nancymaae.fasta | bbduk.sh in=stdin.fasta out=$final_fasta/"${batchNum}"_TETHERIN_Nancymaae.fasta ftl=61 ftr=636
cat ~/*CD4_Nancymaae.fasta | bbduk.sh in=stdin.fasta out=$final_fasta/"${batchNum}"_CD4_Nancymaae.fasta ftl=99 ftr=1472
cat ~/*CCR5_Nancymaae.fasta | bbduk.sh in=stdin.fasta out=$final_fasta/"${batchNum}"_CCR5_Nancymaae.fasta ftl=83 ftr=1141
cat ~/*APOBEC3G_Nancymaae.fasta | bbduk.sh in=stdin.fasta out=$final_fasta/"${batchNum}"_APOBEC3G_Nancymaae.fasta ftl=53 ftr=1201
#below is for throwing in retrogene dup if using retro gene as additional gene in reference file for a3g filtering
#cat ~/*A3G2dup_Nancymaae.fasta | bbduk.sh in=stdin.fasta out=$final_fasta/"${batchNum}"_APOBEC3G_Nancymaae.fasta ftl=53 ftr=1201
#linerize FASTAs
for fastas in $(ls $final_fasta/*.fasta); do
	fasta_formatter -w 0 -i $fastas -o ${fastas//.fasta/_final.fasta}
done
#housekeeping
rm ~/*.fasta
rm $out_fasta/*mdup.fasta
rm $final_fasta/*Nancymaae.fasta
#awk CDS trim if we dont want to use bbmap/duk:
#awk '{if (NR%2==0) print substr($0,63); else print $0}' < cat ~/*TRIMCYP_Nancymaae.fasta
# |\
#awk '{if (NR%2==0) print substr($0,); else print $0}' >  $final_fasta/3_batch_TRIMCYP_Nancymaaeawktest.fasta
echo "Variant Calling Pipeline Complete. Check Outputs."
echo $(date)
end=$(date +%s)
startend=$((end-start))
minutes=60
runtime=$((startend/minutes))
echo "Total Runtime (Minutes): " 
printf "%0.1f\n" $runtime
#woo
