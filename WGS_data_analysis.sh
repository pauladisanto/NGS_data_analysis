#Script from Bioinformagician with minor modifications 
#https://www.youtube.com/watch?v=iHkiQvxyr5c&t=16s
##Picard tool: /home/paula/Softwares_externos/picard.jar

mkdir aligned_reads reads scripts results data supporting_files

cd supporting_files  
mkdir hg38
#wget -P \home\... ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
#-P to save it in a sepecific folder

wget -P /home/paula/Curso_NGS/Assignment_5/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P /home/paula/Curso_NGS/Assignment_5/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

#download reference files 

wget -P /home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip /home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/hg38.fa.gz

# index ref - .fai file before running haplotype caller
samtools faidx /home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/hg38.fa

# ref dict - .dict file before running haplotype caller (creating dictionary file)
gatk CreateSequenceDictionary R=/home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/hg38.fa O=/home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/hg38.dict

# download known sites files for BQSR from GATK resource bundle
wget -P /home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# VARIANT CALLING STEPS 

# directories
ref="/home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/hg38.fa"
known_sites="/home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/home/paula/Curso_NGS/Assignment_5/aligned_reads"
reads="/home/paula/Curso_NGS/Assignment_5/reads"
results="/home/paula/Curso_NGS/Assignment_5/results"
data="/home/paula/Curso_NGS/Assignment_5/data"
trimgalore="/home/paula/Curso_NGS/Assignment_5/Trimgalore"

# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

multiqc . #run it in the folder reads

# No trimming required, quality looks okay, but just to play I will do the trimming


# -------------------
# STEP 1: QC - Trimming 
# -------------------

#Run TrimGalore
#unzip the fastq files
#gunzip -c yourfile.gz > yourfile

gunzip ${reads}/*.fastq.gz #be careful and not to remove the gz (files)
mkdir Trimgalore
trim_galore -q 30 ${reads}/*.fastq -o ${trimgalore}/


# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}


# BWA alignment
#You can use the trimmed fastq files intead
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

#-t 4: Specifies the number of threads (CPU cores) to use
#-R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634": Adds a read group (RG) header line to the output SAM/BAM file. This information is useful for downstream analysis and identification of the origin of reads



# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



SRR062634.paired.sam
# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 


#gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file /home/paula/Curso_NGS/Assignment_5/data/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 




# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf



# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf


#Second part 

#ref="/home/paula/Curso_NGS/Assignment_5/supporting_files/hg38/hg38.fa"
#results="/home/paula/Curso_NGS/Assignment_5/results"

#to filter the vcf files follow the specification in gatk https://gatk.broadinstitute.org/hc/en-us

# Filter SNPs
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_snps.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"



# Filter INDELS
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_indels.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"

#Select Variants that PASS filters.
gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_snps.vcf \
	-O ${results}/analysis-ready-snps.vcf



gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_indels.vcf \
	-O ${results}/analysis-ready-indels.vcf


# to exclude variants that failed genotype filters
cat analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > analysis-ready-snps-filteredGT.vcf
cat analysis-ready-indels.vcf| grep -v -E "DP_filter|GQ_filter" > analysis-ready-indels-filteredGT.vcf

#to download the 
#gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download
#gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download


# -------------------
# Annotate Variants - GATK4 Funcotator
# -------------------

# Annotate using Funcotator by default uses gencode
gatk Funcotator \
	--variant ${results}/analysis-ready-snps-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /home/paula/Softwares_externos/GERMLINE_funcotator/funcotator_dataSources.v1.7.20200521g \
	--output ${results}/analysis-ready-snps-filteredGT-functotated.vcf \
	--output-file-format VCF



gatk Funcotator \
	--variant ${results}/analysis-ready-indels-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /home/paula/Softwares_externos/GERMLINE_funcotator/funcotator_dataSources.v1.7.20200521g \
	--output ${results}/analysis-ready-indels-filteredGT-functotated.vcf \
	--output-file-format VCF



# Extract fields from a VCF file to a tab-delimited table

gatk VariantsToTable \
	-V ${results}/analysis-ready-snps-filteredGT-functotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O ${results}/output_snps.table


# sed 's/|/\t/g'  replace the | separation changing the new separation into a tab separation. The line below only selects the heathers
cat analysis-ready-snps-filteredGT-functotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > output_curated_variants.txt

#
cat output_snps.table | cut -f 5 | grep "NBPF1" | sed 's/|/\t/g' >> output_curated_variants.txt


