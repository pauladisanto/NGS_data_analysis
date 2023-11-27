ssh -Y user7@antares.sa.gu.lcl


#to load the packages
module load _program_name/version
module load annovar/20200608 bcftools/1.17 bwa-mem2/2.2.1 delly/1.1.6 fastqc/0.12.1 gatk/4.4.0.0 igv/2.16.1 multiqc/1.14 R/4.3.1 samtools/1.17 picard/3.0.0 trimgalore/0.6.10
module list

#Create directoy
mkdir Exome
#Create a soft link 
ln -s /home/ftp_courses/NGS/Exome/Fastq/* .
ln -s /home/ftp_courses/NGS/Exome/db/* .
#=====================================================================================================
#QC, filtering and mapping

#for i in *fastq; do fastqc $i; done
fastqc *fastq
#to see all the fastqc reports at the same time
multiqc .
#Run TrimGalore
mkdir Trimgalore
trim_galore -q 30 *.fastq -o /home/user7/Exome/Trimgalore
cd Trimgalore
fastqc *_trimmed.fq
#to see all the fastqc reports at the same time
multiqc .


#download multiqc report
#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/Trimgalore/multiqc_report.html .
#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/multiqc_report.html .

#=====================================================================================================
#Map the filtered reads towards chr11.fa using bwa-mem2
mkdir Alignment 

#bwa-mem2 mem /home/user7/Exome/chr11.fa /home/user7/Exome/Trimgalore/Sample.fq  > Sample.bwa.sam

reference="/home/user7/Exome/chr11.fa"
trim_dir="/home/user7/Exome/Trimgalore"

samples=("AS_2" "AS_ctrl")
read_types=("R1" "R2")

for sample in "${samples[@]}"; do
    for read_type in "${read_types[@]}"; do
        input_fq="${trim_dir}/${sample}.chr11.${read_type}_trimmed.fq"
        output_sam="${sample}.chr11.${read_type}_trimmed.bwa.sam"

        bwa-mem2 mem "${reference}" "${input_fq}" > "${output_sam}"
    done
done



#Convert it into BAM, sort it and index it, using samtools

#samtools view -bS  Sample.sam -o  Sample.bam

samples=("AS_2" "AS_ctrl")
read_types=("R1" "R2")

for sample in "${samples[@]}"; do
    for read_type in "${read_types[@]}"; do
        input_sam="${sample}.chr11.${read_type}_trimmed.bwa.sam"
        output_bam="${sample}.chr11.${read_type}_trimmed.bwa.bam"

        samtools view -bS "${input_sam}" -o "${output_bam}"
    done
done


#Run AddOrReplaceReadGroups from picard using the following parameters
#Picard make the indexes the outputs are .rg.sort.bwa.bai and .rg.sort.bwa.bam files

picard_jar="/apps/bio/software/picard/3.0.0/picard.jar"
samples=("AS_2" "AS_ctrl")
read_types=("R1" "R2")

for sample in "${samples[@]}"; do
    for read_type in "${read_types[@]}"; do
        input_bam="${sample}.chr11.${read_type}_trimmed.bwa.bam"
        output_bam="${sample}.chr11.${read_type}_trimmed.rg.sort.bwa.bam"

        java -jar "${picard_jar}" AddOrReplaceReadGroups \
            I="${input_bam}" \
            O="${output_bam}" \
            RGID=1 \
            RGLB=TruSeq \
            RGPL=Illumina \
            RGPU=CGG \
            RGSM="${sample}.chr11.${read_type}" \
            SO=coordinate \
            CREATE_INDEX=true
    done
done


#FLAGSTAT
#samtools flagstat your_aligned.bam > output.flagstat


#IDXSTATS
#samtools idxstats your_aligned.bam > output.idxstats

samples=("AS_2" "AS_ctrl")
read_types=("R1" "R2")

for sample in "${samples[@]}"; do
    for read_type in "${read_types[@]}"; do
        input_bam="${sample}.chr11.${read_type}_trimmed.rg.sort.bwa.bam"
        flagstat_output="${sample}.chr11.${read_type}_trimmed.rg.sort.flagstat"
        idxstats_output="${sample}.chr11.${read_type}_trimmed.rg.sort.idxstats"

        # Run samtools flagstat
        samtools flagstat "${input_bam}" > "${flagstat_output}"

        # Run samtools idxstats
        samtools idxstats "${input_bam}" > "${idxstats_output}"
    done
done


#=======================================================================================================

#On-target coverage

#Since the data was from a targeted approach, let's inspect the on-target coverage plot. GATK's DepthOfCoverage can help us in this task. 

#gatk DepthOfCoverage \
 #  -R db/chr11.fa \
 #  -I YOUR_SAMPLE.metrics  \
 #  -O YOUR_SAMPLE_NAME \
 #  -L chr11.exon.bed \
 #  --output-format TABLE \
 #  --summary-coverage-threshold 10 \
 #  --summary-coverage-threshold 30 \
 #  --summary-coverage-threshold 50    

   #-R db/chr11.fa: Specifies the reference genome file (chr11.fa) against which the sequencing data will be aligned.
   #-I YOUR_SAMPLE.metrics: Specifies the input BAM file (YOUR_SAMPLE.metrics). This BAM file should contain the alignment information of 
   # your sequencing data.
   #-O YOUR_SAMPLE_NAME: Specifies the output prefix or directory name (YOUR_SAMPLE_NAME) for the output coverage files. 
   # The tool will generate multiple output files with this prefix.
   #-L chr11.exon.bed: Specifies the intervals to include in the coverage analysis. 
   # In this case, it's a BED file (chr11.exon.bed) containing the coordinates of exonic regions on chromosome 11.
   #--output-format TABLE: Specifies the output format for the coverage results. In this case, it's set to TABLE format.
   #--summary-coverage-threshold 10, 30, 50: Specifies the coverage thresholds for generating summary statistics. 
   #The tool will calculate coverage at these thresholds (10x, 30x, and 50x) and provide summary statistics for each threshold.


# List of sample names
samples=("AS_2" "AS_ctrl")

# Loop through each sample
for sample in "${samples[@]}"; do
    # Run gatk DepthOfCoverage for R1
    gatk DepthOfCoverage \
        -R /home/ftp_courses/NGS/Exome/db/chr11.fa \
        -I "${sample}.chr11.R1_trimmed.rg.sort.bwa.bam" \
        -O "${sample}.chr11.R1_trimmed.rg.sort_gatk" \
        -L /home/ftp_courses/NGS/Exome/db/chr11.exon.bed \
        --output-format TABLE \
        --summary-coverage-threshold 10 \
        --summary-coverage-threshold 30 \
        --summary-coverage-threshold 50

    # Run gatk DepthOfCoverage for R2
    gatk DepthOfCoverage \
        -R /home/ftp_courses/NGS/Exome/db/chr11.fa \
        -I "${sample}.chr11.R2_trimmed.rg.sort.bwa.bam" \
        -O "${sample}.chr11.R2_trimmed.rg.sort_gatk" \
        -L /home/ftp_courses/NGS/Exome/db/chr11.exon.bed \
        --output-format TABLE \
        --summary-coverage-threshold 10 \
        --summary-coverage-threshold 30 \
        --summary-coverage-threshold 50
done


#I move all the GATK results in to a new folder
mkdir GATK
#mv /path/to/source/folder/*.bam /path/to/destination/folder/
mv /home/user7/Exome/Alignment/*gatk* /home/user7/Exome/Alignment/GATK


#To better understand this data, let's make a plot with this information.
#First we transpose our file (changing columns to rows) and then, we convert it to percentage

#tail -1 SAMPLE_NAME.sample_cumulative_coverage_proportions | \
#sed 's/\t/\n/g' | \
#tail -n +2 | 
#awk 'OFS="\t" {print $1*100,num++}' > SAMPLE_NAME.proportions


# List of sample names
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed")

# Loop through each sample
for sample in "${samples[@]}"; do
    # Run the desired commands for each sample
    tail -1 "${sample}.rg.sort_gatk.sample_cumulative_coverage_proportions" | \
    sed 's/\t/\n/g' | \
    tail -n +2 | \
    awk 'OFS="\t" {print $1*100,num++}' > "${sample}.rg.sort_gatk.proportions"
done

#==============================================================================================
#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/Alignment/GATK/*trimmed.rg.sort_gatk.proportions .
#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/Alignment/GATK/*sample_interval_summary .
#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/Alignment/GATK/*.sample_cumulative_coverage_proportions .

#Visualization in R

# List of sample names
samples <- c("AS_2.chr11.R1_trimmed", "AS_2.chr11.R2_trimmed", "AS_ctrl.chr11.R1_trimmed", "AS_ctrl.chr11.R2_trimmed")

# Loop through each sample
for (sample in samples) {
  # Reading the data
  data <- read.table(paste0(sample, ".rg.sort_gatk.proportions"))

  # Setting name of figure
  png(paste0(sample, ".rg.sort_gatk.coverage.png"))

  # Scatter plot
  plot(data$V2, data$V1, type="l", log="x", ylim=c(0, 100),
       ylab="% On-target", xlab="Depth", main=paste("On-target coverage", sample, ".rg.sort_gatk"))

  # Adding lines at different coverages
  abline(v=c(10, 20, 30, 50), col="gray", lty=2)

  # Adding lines at different % of on-target
  abline(h=c(50, 70, 80), col=c(51, "gray", 51), lty=2)

  # Adding X and % within the graph
  text(c(8, 17, 26, 43), c(100, 100, 100, 100), c("10x", "20x", "30x", "50x"), col="gray")
  text(c(1, 1, 1), c(47, 67, 77), c("50%", "70%", "80%"), col=c(51, "gray", 51))

  # Closing graphical device
  dev.off()
}

#==============================================================================================
#PCR duplicates
#When dealing with variant calling, an important step is to identify possible PCR duplicates, 
#this to remove any bias when we are doing the actual calling.


#-I             -> YOUR_SAMPLE.bwa.sort.fix.bam 
#-O             -> YOUR_SAMPLE.bwa.sort.fix.mkdup.bam
#-M             -> YOUR_SAMPLE.metrics 
#--CREATE_INDEX -> TRUE

mkdir Duplicates
pdw # to get the path of the current working directory

# Array of sample names
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed")

# Picard tool path
picard_path="/apps/bio/software/picard/3.0.0/picard.jar"

# Output directory
output_dir="/home/user7/Exome/Alignment/Duplicates/"

# Loop through each sample
for sample in "${samples[@]}"; do
  input_bam="${sample}.rg.sort.bwa.bam"
  output_bam="${output_dir}${sample}.rg.sort.bwa.mkdup.bam"
  metrics_file="${output_dir}${sample}.rg.sort.bwa.metrics"

  # Run MarkDuplicates
  java -jar "$picard_path" MarkDuplicates \
    I="$input_bam" \
    O="$output_bam" \
    M="$metrics_file" \
    CREATE_INDEX=true
done


#Flagstat and idxstats metrics for this new alignment file
#multiqc to include these statistics

# Array of sample names
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed")

# Output directory
output_dir="/home/user7/Exome/Alignment/Duplicates/"

# Loop through each sample
for sample in "${samples[@]}"; do
  input_bam="${output_dir}${sample}.rg.sort.bwa.mkdup.bam"

  # Run samtools flagstat
  samtools flagstat "$input_bam" > "${output_dir}${sample}.rg.sort.bwa.mkdup.flagstat"

  # Run samtools idxstats
  samtools idxstats "$input_bam" > "${output_dir}${sample}.rg.sort.bwa.mkdup.idxstats"
done

# Run MultiQC on the output directory
multiqc "$output_dir"

#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/Alignment/Duplicates/*multiqc* .

#Base Quality Score Recalibration (BQSR)
#This pre-processing step detects systematic errors made by the sequencing machine when it estimates the accuracy of each base call. 
#Most of the short variant calling tools rely on the base quality score, so it is important to correct these systematic technical errors.

#First, we will produce a recalibration file from our input data using the BaseRecalibrator tool. 
#To speed up the processing time, let's generate a file containing only SNPs from chr11. This is because I am just analysing the Chr11

mkdir BQSR
# keeping the header of the VCF file
grep ^#  /home/ftp_courses/NGS/Exome/db/dbsnp_138.b37.vcf | \
sed 's/ID=11/ID=chr11/' | \
grep -P  "contig=<ID=[1-9,A-Z]" -v | \
grep -P  "<ID=GL" -v | \
sed 's/135006516/57311868/' > title  

# extracting SNPs in chr11
grep ^11 -w  /home/ftp_courses/NGS/Exome/db/dbsnp_138.b37.vcf | \
awk '$2 < 57311868 {print $0}' | \
sed 's/^11/chr11/'  > 11

# merging title and data, formatting to add "chr", removing unwanted rows
cat title 11 > dbsnp_138.b37.chr11.vcf

#Index the file using gatk IndexFeatureFile and:
#-I  -> dbsnp_138.b37.chr11.vcf

gatk IndexFeatureFile -I dbsnp_138.b37.chr11.vcf

#Run the gatk BaseRecalibrator for both samples, with :
#-R            -> chr11_reference_file 
#-I            -> YOUR_SAMPLE.bwa.sort.fix.mkdup.bam 
#-O            -> YOUR_SAMPLE.table 
#--known-sites -> THE_VCF_FILE_YOU_JUST_CREATED

# Define the paths and files
reference="/home/user7/Exome/chr11.fa"
known_sites="dbsnp_138.b37.chr11.vcf"
input_dir="/home/user7/Exome/Alignment/Duplicates"
output_dir="/home/user7/Exome/BQSR"

# List of samples
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed")

# Run BaseRecalibrator for each sample
for sample in "${samples[@]}"; do
    input_bam="$input_dir/${sample}.rg.sort.bwa.mkdup.bam"
    output_table="$output_dir/${sample}.table"

    gatk BaseRecalibrator \
        -R "$reference" \
        -I "$input_bam" \
        -O "$output_table" \
        --known-sites "$known_sites"
done

#=================================================================================================

#Now let's correct the scores of the sample, running gatk ApplyBQSR for hte alignments

#    -R                -> chr11_reference_file 
#    -I                -> YOUR_SAMPLE.bwa.sort.fix.mkdup.bam 
#    -O                -> YOUR_SAMPLE.sort.fix.mkdup.recal.bam  
#    --bqsr-recal-file -> THE_TABLE_FILE_YOU_JUST_CREATED 


reference="/home/user7/Exome/chr11.fa"
input_dir="/home/user7/Exome/Alignment/Duplicates"
output_dir="/home/user7/Exome/BQSR"

# List of samples
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed") 

# Run ApplyBQSR for each sample
for sample in "${samples[@]}"; do
    input_bam="$input_dir/${sample}.rg.sort.bwa.mkdup.bam"
    output_bam="$output_dir/${sample}.rg.sort.bwa.mkdup.recal.bam"
    bqsr_table="$output_dir/${sample}.table"

    gatk ApplyBQSR \
        -R "$reference" \
        -I "$input_bam" \
        -O "$output_bam" \
        --bqsr-recal-file "$bqsr_table"

    echo "Applied BQSR for sample: $sample"
done

#=============================================================================================

#Variant discovery (GATK)
#The Genome Analysis Toolkit or GATK is a software package developed to analyze next-generation resequencing data, 
#focusing on variant discovery and genotyping. For the variant calling we will use the HaplotypeCaller, 
#which is an SNP/indel caller that uses a Bayesian genotype likelihood model to estimate simultaneously 
#the most likely genotypes and allele frequency in a population of N samples.
#Run gatk HaplotypeCaller to do the variant calling, with these parameters:
#-R -> chr11_reference_file 
#-I -> YOUR_SAMPLE.bwa.sort.fix.mkdup.recal.bam 
#-O -> YOUR_SAMPLE.gatk.raw.vcf  

mkdir VCFs

reference="/home/user7/Exome/chr11.fa"
input_dir="/home/user7/Exome/BQSR"
output_dir="/home/user7/Exome/VCFs"

# List of samples
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed")

# Run HaplotypeCaller for each sample
for sample in "${samples[@]}"; do
    input_bam="$input_dir/${sample}.rg.sort.bwa.mkdup.recal.bam"
    output_vcf="$output_dir/${sample}.gatk.raw.vcf"

    gatk HaplotypeCaller \
        -R "$reference" \
        -I "$input_bam" \
        -O "$output_vcf"

    echo "Called variants for sample: $sample"
done

#Let's convert the VCF files into a tab delimited files.
#Use gatk VariantsToTable with:

output_dir="/home/user7/Exome/VCFs"

# List of samples
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed")

# Run VariantsToTable for each sample
for sample in "${samples[@]}"; do
    input_vcf="/home/user7/Exome/VCFs/${sample}.gatk.raw.vcf"
    output_table="$output_dir/${sample}.gatk.raw.table"

    gatk VariantsToTable \
        -O "$output_table" \
        -V "$input_vcf" \
        -F CHROM \
        -F POS \
        -F ID \
        -F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -GF AC \
        -GF GT \
        -GF AD \
        -GF DP

    echo "Converted variants to table for sample: $sample"
done

#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/VCFs/*gatk.raw.table .

#SNP-callers emit a number of false-positive mutations, as any prediction program. 
#Some of these false-positives can be detected and rejected by various statistical tests. 
#Filtering or flagging the variant calls in the VCF file is a useful step to detect such false-positives.
#Run gatk VariantFiltration with:

output_dir="/home/user7/Exome/VCFs"

# List of samples
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed")

# Run VariantFiltration for each sample
for sample in "${samples[@]}"; do
    input_vcf="/home/user7/Exome/VCFs/${sample}.gatk.raw.vcf"
    output_filtered_vcf="$output_dir/${sample}.gatk.raw.flag.vcf"

    gatk VariantFiltration \
        -O "$output_filtered_vcf" \
        -V "$input_vcf" \
        --filter-name "LowDP" \
        --filter-expression "DP < 10.0"

    echo "Filtered variants for sample: $sample"
done


#These are hard filters that are applied when a variant callset is too small or when training sets are not available. 

#======================================================================================================
#Variant annotation (Annovar)

#Now we need to annotate our variants in terms of genomic context, aminoacid change and effect, 
#previous mutation knowledge, etc. so we can evaluate and filter mutations that are irrelevant.

#ANNOVAR is a tool that annotates genetic variants. The overall procedure is to create an input file 
#from the VCF file and then annotate the variant positions with information from several databases.

#Reformat both alignment files with convert2annovar.pl using:

#-format  -> vcf4 
#-outfile -> YOUR_FILE.gatk.raw.flag.annovar
#YOUR_FILE.gatk.raw.flag.vcf 

mkdir AnnovarFiles
# List of samples
samples=("AS_2.chr11.R1_trimmed" "AS_2.chr11.R2_trimmed" "AS_ctrl.chr11.R1_trimmed" "AS_ctrl.chr11.R2_trimmed")

# Run convert2annovar.pl for each sample
for sample in "${samples[@]}"; do
    input_vcf="/home/user7/Exome/VCFs/${sample}.gatk.raw.flag.vcf"
    output_annovar="/home/user7/Exome/AnnovarFiles/${sample}.gatk.raw.flag.annovar"

    convert2annovar.pl -format vcf4 -outfile "$output_annovar" "$input_vcf"

    echo "Converted VCF to ANNOVAR format for sample: $sample"
done


#Add the annotation to the alignments using table_annovar.pl:

#  -buildver  -> hg19 
#  -protocol  -> dbsnp_137_chr11_gatk,refGene,ljb2_all,1000g2012apr_all,clinvar_20200316 
#  -operation -> f,g,f,f,f 
#  -outfile   -> YOUR_FILE.gatk.raw.flag.annovar.table 
#  --nastring -> - 
#  --remove 
#  YOUR_FILE.gatk.raw.flag.annovar 
#  $ANNOVAR_HOME/humandb 

samples=("AS_2" "AS_ctrl")
read_types=("R1" "R2")

for sample in "${samples[@]}"; do
    for read_type in "${read_types[@]}"; do
        input_annovar="${sample}.chr11.${read_type}_trimmed.gatk.raw.flag.annovar"
        output_table="${sample}.chr11.${read_type}_trimmed.gatk.raw.flag.annovar.table"

        table_annovar.pl \
            "${input_annovar}" \
            "$ANNOVAR_HOME/home/ftp_courses/NGS/Exome/db/humandb" \
            -buildver hg19 \
            -protocol dbsnp_137_chr11_gatk,refGene,ljb2_all,1000g2012apr_all,clinvar_20200316 \
            -operation f,g,f,f,f \
            -outfile "${output_table}" \
            --nastring - \
            --remove
    done
done

#This generates the files sample_gatk.raw.flag.annovar.table.hg19_multianno.txt

#==========================================================================================================
#Merge the *hg19_multianno.txt file with the *gatk.raw.flag.table , and save it in a file with extension .tsv. 

# Specify the file names
table_file="AS_2.chr11.R1_trimmed.gatk.raw.flag.annovar"
anno_file="AS_2.chr11.R1_trimmed.gatk.raw.flag.annovar.table.hg19_multianno.txt"

# Specify the output file name
output_file="merged_output.tsv"


# Use paste to merge the two files
paste "$anno_file" "$table_file" > "$output_file"

# Use paste to merge the two files
paste "$table_file" "$anno_file" > "$output_file"

#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/VCFs/*merged_output.tsv .


#==========================================================================================================
#When you have multiple files and positions to investigate in IGV you can save time by creating a bat script. 
#The following script take snapshots of the coordinate chr11:1017135 across AS_2 and Control using both their bams and mutation files.
#new 
#genome hg19 
#load PATH_TO_FILE/AS_2.bwa.sort.fix.mkdup.recal.bam 
#load PATH_TO_FILE/AS_2.gatk.raw.flag.vcf
#load PATH_TO_FILE/AS_ctrl.bwa.sort.fix.mkdup.recal.bam 
#load PATH_TO_FILE/AS_ctrl.gatk.raw.flag.vcf
#snapshotDirectory . 
#squish 
#goto chr11:1017135 
#snapshot 1017135_Exomedata.png 
#exit  
#igv.sh -b bashIGV.bat


new 
genome hg19 
load /home/user7/Exome/BQSR/AS_2.chr11.R1_trimmed.rg.sort.bwa.mkdup.recal.bam
load /home/user7/Exome/VCFs/AS_2.chr11.R1_trimmed.gatk.raw.flag.vcf
load /home/user7/Exome/BQSR/AS_2.chr11.R2_trimmed.rg.sort.bwa.mkdup.recal.bam
load /home/user7/Exome/VCFs/AS_2.chr11.R2_trimmed.gatk.raw.flag.vcf
load /home/user7/Exome/BQSR/AS_ctrl.chr11.R1_trimmed.rg.sort.bwa.mkdup.recal.bam
load /home/user7/Exome/VCFs/AS_ctrl.chr11.R1_trimmed.gatk.raw.flag.vcf
load /home/user7/Exome/BQSR/AS_ctrl.chr11.R2_trimmed.rg.sort.bwa.mkdup.recal.bam
load /home/user7/Exome/VCFs/AS_ctrl.chr11.R2_trimmed.gatk.raw.flag.vcf
snapshotDirectory . 
squish 
goto chr11:6411935 
snapshot SMPD1_6411935_Exomedata.png 
exit  
#=============================================================================
goto chr11:17408630 
snapshot KCNJ11_17408630_Exomedata.png 
exit  

#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/IGV/*png .
#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/BQSR/*_trimmed.rg.sort.bwa.mkdup.recal.bam .
#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/VCFs/*_trimmed.gatk.raw.flag.vcf .


#======================================================================================================================================

#Now let's analyse WGS data to see if we can identify any structural variants. When analyzing tumor samples you normally have two datasets; 
#one normal (germline) sample and one tumor sample (somatic). I will only work with two tumor samples from chromosome 17, 
#these are breast cancer cell lines: HCC1143 and HCC1954.
#Create a soft link to the entire folder /home/ftp_courses/NGS/Exome/SV
mkdir Structural_variants
ln -s /home/ftp_courses/NGS/Exome/SV/* .

#Before running any variant calling, we need to mark duplicates.
#Run MarkDuplicates for each alignment file, with the following parameters:
#-I -> YOUR_SAMPLE.tumor.bam 
#-O -> YOUR_SAMPLE.tumor.mkdup.bam 
#-M -> YOUR_SAMPLE.tumor.metrics 
#--CREATE_INDEX -> true

java -Xmx4g -jar /apps/bio/software/picard/3.0.0/picard.jar MarkDuplicates \
    I=HCC1143.tumor.bam\
    O=HCC1143.tumor.mkdup.bam \
    M=HCC1143.tumor.metrics \
    CREATE_INDEX=true



java -Xmx4g -jar /apps/bio/software/picard/3.0.0/picard.jar MarkDuplicates \
    I=HCC1954.tumor.bam\
    O=HCC1954.tumor.mkdup.bam \
    M=HCC1954.tumor.metrics \
    CREATE_INDEX=true

#Run the variant calling with delly call to find any structural variants:
#delly call \
#-g Homo_sapiens_assembly19.fasta_WITHIN_/db
#-o ALL.tumor.bcf
#list_of_alignments_separated_by_space 
#list_of_alignments_separated_by_space: Replace this with the actual list of your alignment files (BAM files) separated by spaces. 
#These are the files you want to use for variant calling.

delly call \
-g /home/ftp_courses/NGS/Exome/db/Homo_sapiens_assembly19.fasta  \
-o ALL.tumor.bcf \
HCC1143.tumor.mkdup.bam  HCC1954.tumor.mkdup.bam

#The default output is in bcf format.
#Convert the file to vcf format using bcftools:

bcftools view ALL.tumor.bcf > ALL.tumor.vcf

#Take a look at the vcf file.

bcftools view -i 'FILTER=="PASS"' ALL.tumor.vcf > ALL.tumor_filtered.vcf

#Remove all variants annotates as IMPRECISE

bcftools view -i 'SVTYPE=="PRECISE;SVTYPE="' ALL.tumor_filtered.vcf > ALL.tumor_filtered_PRECISE.vcf



#scp -p user7@antares.sa.gu.lcl:/home/user7/Exome/Structural_variants/*png .
