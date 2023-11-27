ssh -Y user7@antares.sa.gu.lcl


#to load the packages
module load _program_name/version
module load blast+/2.14 fastqc/0.12.1 igv/2.16.1 multiqc/1.14 rseqc/5.0.1 R/4.3.1 samtools/1.17 subread/2.0.6

#Alignment QC

#As you have already performed the quality filtering and alignment steps in previous exercises, we are taking on the analysis from the alignment files. So, just to refresh our minds:
#Q. What are the typical steps you perform to obtain the alignment file?

#Note: Due to time limit, we will be working with a subset of the data, that corresponds to chr21. Besides, we will be practicing using the for loop to automate tasks since we are dealing with 12 samples and it would be a lot of typing if you do one sample at the time!
#Create an RNAseq folder in your home directory
#You can directly create a soft link to the actual folders, but let's practice the for loop
#Within the directory, create the folders db, Alignment and Counts.

mkdir RNAseq
mkdir db
mkdir Alignment 
mkdir Counts


#Create links to:

#Alignment files from /home/ftp_courses/NGS/RNAseq/Alignment_chr21/ to your Alignment directory. You can use a for loop for this:

for i in /home/ftp_courses/NGS/RNAseq/Alignment_chr21/*; do ln -s $i Alignment/; done

#Database files from /home/ftp_courses/NGS/RNAseq/db to your db directory. Use a for loop as in the example above.
# unlink symbolic_link_name

for i in /home/ftp_courses/NGS/RNAseq/db/*; do ln -s $i db/; done

#Alignment QC with samtools

#Samtools is a great tool to gather the alignment statistics.
#Go into your Alignment directory and index your BAM files:
#NOTE: Remember that we can use a loop to run a command on multiple samples. See example on how to do it below

for i in *bam; do samtools index $i; done

#Create statistics files using flagstat and idxstats for each sample:

for i in  *.bam; do samtools flagstat $i > ${i%.bam}.flagstat; done

#NOTE: ${i%.bam}.flagstat removes the suffix .bam from the output name and replaces it with the suffix ".flagstat"

for i in  *.bam; do samtools idxstats $i > ${i%.bam}.idxstats; done

#Run fastqc on each sample


#Merge the results with multiqc

#Have a look at the multiqc HTML report.

multiqc_report.html

#scp -p user7@antares.sa.gu.lcl:/home/user7/RNAseq/multiqc_report.html .

#The header of the alignment file contains information of used chromosomes, as well as command lines used to generate the BAM file among others. 
#This can be handy if you don't know how the alignment was created.

#Use samtools view with the -H flag to display the header

#samtools view -H YOUR_BAM_FILE

samtools view -H SRS308866Aligned.sortedByCoord.out_chr21.bam

#Alignment QC using RSeQC

#RSeQC is a program consisting of multiple python scripts that can give you information about the quality of your data.

#The script infer_experiment.py from RSeQC provides information about what kind of library protocol has been used 
#(paired end, single end, stranded, non stranded).
#Before running RSeQC, let's prepare the reference gene model file (-r):
#hg38_RefSeq.bed has the annotation for all genes, but since we are just focusing on chr21, it will be quicker to have a subset of the bed file.

#Create a new file with only genes belonging to chr21 (use grep or awk)
#Save it as hg38_RefSeq_chr21.bed under your db directory

awk '$1 == "chr21" || $1 == "21"' hg38_RefSeq.bed > hg38_RefSeq_chr21.bed

#Run infer_experiment.py. Try to use a for loop to automate your analysis:

#infer_experiment.py -r LOCATION_OF_YOUR_FILE/hg38_RefSeq_chr21.bed -i YOUR_BAM_FILE > YOUR_BAM_FILE.infer_experiment
#SRS308866Aligned.sortedByCoord.out_chr21.bam  SRS308874Aligned.sortedByCoord.out_chr21.bam  SRS308882Aligned.sortedByCoord.out_chr21.bam
#SRS308868Aligned.sortedByCoord.out_chr21.bam  SRS308876Aligned.sortedByCoord.out_chr21.bam  SRS308883Aligned.sortedByCoord.out_chr21.bam
#SRS308870Aligned.sortedByCoord.out_chr21.bam  SRS308878Aligned.sortedByCoord.out_chr21.bam  SRS308885Aligned.sortedByCoord.out_chr21.bam
#SRS308872Aligned.sortedByCoord.out_chr21.bam  SRS308880Aligned.sortedByCoord.out_chr21.bam  SRS308887Aligned.sortedByCoord.out_chr21.bam
#!/bin/bash

# Assuming all BAM files are in the current directory and follow the pattern SRS*.bam
for bam_file in SRS*.bam; do
    output_file="${bam_file}.infer_experiment"

    infer_experiment.py -r /home/user7/RNAseq/db/hg38_RefSeq_chr21.bed -i "$bam_file" > "$output_file"

    echo "Infer experiment results for $bam_file written to $output_file"
done

#Sometimes you can have problems with the RNA being fragmented (e.g degraded). You can detect this by having more coverage over the one 
#end of the gene compared to the other. You can visualize this by using geneBody_coverage.py:

#geneBody_coverage.py -i BAM_FILE1,BAM_FILE2,BAM_FILE3,... -r LOCATION_OF_YOUR_FILE/hg38_RefSeq_chr21.bed -o geneBodyCoverage 

geneBody_coverage.py -i  SRS*.bam -r /home/user7/RNAseq/db/hg38_RefSeq_chr21.bed -o geneBodyCoverage 

#scp -p user7@antares.sa.gu.lcl:/home/user7/RNAseq/Alignment/geneBodyCoverage.geneBodyCoverage.curves.pdf .

#Gene counts
#The next step is to count the amount of reads aligned towards the genes in the reference genome so we can assess gene expression. 
#Let's extract the gene counts for your alignment files of chr21.

#There are multiple tools to extract the gene counts, here we will be using featureCounts, 
#one program from the Subread package which comprises a suite of software programs for processing next-gen sequencing read data.

#Run featureCounts on your bam files.
#Use Homo_sapiens.GRCh38.109.chr21.gtf as the reference gtf file, which contains the feature (exon/genes) coordinates in the reference genome.
#Don't forget to set the flags for paired-end data and for unstranded data
#Save the results in your Counts directory


#featureCounts -p --countReadPairs -s FLAG_FOR_UNSTRANDED_DATA -t exon -a LOCATION_OF_YOUR_FILE/Homo_sapiens.GRCh38.109.chr21.gtf -o COUNTS_DIRECTORY/21.counts YOUR_BAM_FILES_SEPARATED_BY_SPACE

featureCounts -p --countReadPairs -s 0 -t exon -a /home/user7/RNAseq/db/Homo_sapiens.GRCh38.109.chr21.gtf -o /home/user7/RNAseq/Counts/21.counts SRS308866Aligned.sortedByCoord.out_chr21.bam SRS308874Aligned.sortedByCoord.out_chr21.bam SRS308882Aligned.sortedByCoord.out_chr21.bam SRS308868Aligned.sortedByCoord.out_chr21.bam SRS308876Aligned.sortedByCoord.out_chr21.bam SRS308883Aligned.sortedByCoord.out_chr21.bam SRS308870Aligned.sortedByCoord.out_chr21.bam SRS308878Aligned.sortedByCoord.out_chr21.bam SRS308885Aligned.sortedByCoord.out_chr21.bam SRS308872Aligned.sortedByCoord.out_chr21.bam SRS308880Aligned.sortedByCoord.out_chr21.bam SRS308887Aligned.sortedByCoord.out_chr21.bam
#-s FLAG_FOR_UNSTRANDED_DATA: Specify the strand specificity. 
#Replace FLAG_FOR_UNSTRANDED_DATA with one of the following options based on your data:

 #   0 or unstranded: unstranded data
 #   1 or stranded: stranded data
 #   2 or reversely_stranded: reverse stranded data

#For paired-end data, you include the -p option. For example:
mkdir DESeq2
#DESeq2, one of the statistical packages to perform DE analysis, needs a count matrix of the raw counts as input. 
#A count matrix is a table where the columns represents the samples and the rows represents the genes. 
#Create such a matrix from the count file from the featureCounts step.

#copy the R things to my computer
scp -p user7@antares.sa.gu.lcl:/home/user7/RNAseq/Counts/* .

scp -p user7@antares.sa.gu.lcl:/home/ftp_courses/NGS/RNAseq/forDE/* .

#Now you will use DESeq2, an R package that estimates variance-mean dependence in count data from high-throughput sequencing 
#assays and tests for differential expression based on a model using the negative binomial distribution. 
#In this exercise you will try some different design options.
#You have created a count matrix that contains only chr21. For this part of the analysis you will use the entire dataset.
#Copy the files all_counts.txt and sample_description.txt from /home/ftp_courses/NGS/RNAseq/forDE/ to your local computer. Remember you can use scp <origin> <target>.
#Start R in your local computer.
#Select your working directory (where you have the all_counts.txt file) by clicking:
#Session -> Set Working Directory -> Choose directory
#Load the DESeq2 package
#Read the count data (all_counts.txt) and information about the samples (sample_description.txt)
# loading library

library(DESeq2)
# Not installed? Run:
# install.packages("ggplot2")
library(ggplot2)
# Not installed? Run: 
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2")

# reading data
counts <- read.delim("all_counts.txt", row.names=1)
# look at the data
head(counts)

# read metadata
samples <-read.delim("sample_description.txt")

# check what information you have about the samples
samples

# make sure patient number is not numeric
samples$patient <-as.factor(samples$patient)

### Single factor design  

# Run ``DESeqDataSetFromMatrix()`` using ``design = ~treatment``, so we only take the treatment into account for comparing the samples 

ds <- DESeqDataSetFromMatrix(countData = counts, 
                             colData   = samples, 
                             design    = ~treatment)
# check data
data.frame(colData(ds))                                  

# run the DE analysis
dst <- DESeq(ds)                                         
dst

#==========================================================

#Transform the data using vst() (variance stabilizing transformation). This will make the data more suitable for clustering and visualization
#Plot the results of the Principal Component Analysis,using plotPCA()

# blind means the dispersion will not be recalculated
vsd <- vst(dst, blind=FALSE)

# Plotting the PCA
plotPCA(vsd, intgroup=c("patient", "treatment"))

#============================================================
#You can customize the PCA plot, to better visualize the different sample

# obtaining the PCA values
pcaData <- plotPCA(vsd,
                   intgroup   = c("treatment", "patient"),
                   returnData = TRUE)

# Calculating the variance
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plotting the PCA
ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()


#Another way to visualize the clustering of the samples is by applying the dist() function to the transpose of the transformed count matrix 
#to get sample-to-sample distances, and then provide a hierarchical clustering (hc) to the heatmap function based on the sample distances:


# Calculating per sample not per gene
sampleDists <- dist(t(assay(vsd)))

suppressMessages(library("RColorBrewer"))
# Not installed? Run:
#install.packages("RColorBrewer")

# Calculating distances
sampleDistMatrix <- as.matrix(sampleDists)

# Prettifying sample names
rownames(sampleDistMatrix) <- paste(dst$treatment, dst$patient, sep="-")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)

# SELECT ONE PALETTE: Reds, Blues, Greens, Purples, Oranges or Greys
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

suppressMessages(library(pheatmap))
# Not installed? Run:
# install.packages("pheatmap")

# printing heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#========================================================
#Have a look at the results from the DE analysis:

# results generated during the DE analysis
resultsNames(dst)

# extracting DPN_vs_Control
res_DPN <-results(dst, name="treatment_DPN_vs_Control")

# look at 10 most associated for DPN treatment
head(res_DPN[order(res_DPN$padj),], 10)


#Another way to visualize the amount of significant genes is by plotting a histogram of the adjusted pvalues:

res_DPN$padj


#Do the same visualization, but now extracting OHT_vs_Control. Hint: You can just replace each DPN in the code above with OHT

res_OHT <-results(dst, name="treatment_OHT_vs_Control")
head(res_OHT[order(res_OHT$padj),], 10)

res_OHT$padj


# Specify a fixed bin size, e.g., 0.01
bin_size <- 0.01
xmin=0
xmax=1

# Generate breakpoints using seq()
breakpoints <- seq(xmin, xmax, by = bin_size)

# Create the histogram with fixed bin size
hist(
  res_DPN$padj,          # Data to be plotted (vector of values)
  breaks = breakpoints,  # Breakpoints for fixed bin size
  col = "skyblue",       # Fill color of the bars
  border = "slateblue",  # Border color of the bars
  main = "res_DPN",      # Main title of the plot
  xlab = "padj",         # Label for the x-axis
  xlim = c(xmin, xmax)         # x-axis limits
)

# Create the histogram with fixed bin size
hist(
  res_OHT$padj,          # Data to be plotted (vector of values)
  breaks = breakpoints,  # Breakpoints for fixed bin size
  col = "skyblue",       # Fill color of the bars
  border = "slateblue",  # Border color of the bars
  main = "res_OHT",      # Main title of the plot
  xlab = "padj",         # Label for the x-axis
  xlim = c(xmin, xmax)         # x-axis limits
)

#Multifactor design

#Now we will try to do a paired analysis instead, since we have information for the same patient for different treatments.

# setting the design
dsp <- DESeqDataSetFromMatrix(countData = counts, 
                              colData   = samples, 
                              design    = ~patient + treatment)
# running the analysis
dsp <- DESeq(dsp)

# Comparisons automatically generated
resultsNames(dsp)


#Let's extract DPN_vs_Control:
# extracting results
res_DPN_paired <-results(dsp, name="treatment_DPN_vs_Control")  

# selecting with 0.01 as threshold
res_DPN_paired_sig <- res_DPN_paired[ which(res_DPN_paired$padj < 0.01), ]

# checking data
head(res_DPN_paired[order(res_DPN_paired$padj),], 10)

dim(res_DPN_paired[which(res_DPN_paired$padj < 0.05), ])

#Generate a histogram of the adjusted values for this paired design and compare it to the one from the previous design. 


# Specify a fixed bin size, e.g., 0.01
bin_size <- 0.01
xmin=0
xmax=1

# Generate breakpoints using seq()
breakpoints <- seq(xmin, xmax, by = bin_size)


# Create the histogram with fixed bin size
hist(
  res_DPN_paired$padj,          # Data to be plotted (vector of values)
  breaks = breakpoints,  # Breakpoints for fixed bin size
  col = "skyblue",       # Fill color of the bars
  border = "slateblue",  # Border color of the bars
  main = "res_DPN_paired",      # Main title of the plot
  xlab = "padj",         # Label for the x-axis
  xlim = c(xmin, xmax)         # x-axis limits
)

#Check the statistically significant genes:
# selecting with 0.01 as threshold
res_DPN_paired_sig <- res_DPN_paired[ which(res_DPN_paired$padj < 0.01), ]

# checking data
res_DPN_paired_sig

#Follow the same procedure, and extract DPN_vs_Control (you can just replace each DPN in the code above with OHT)
res_OHT_paired <-results(dsp, name="treatment_OHT_vs_Control")  

# selecting with 0.01 as threshold
res_OHT_paired_sig <- res_OHT_paired[ which(res_OHT_paired$padj < 0.01), ]

# checking data
res_OHT_paired_sig


#Investigate if there are common genes:

intersect(rownames(res_DPN_paired_sig),rownames(res_OHT_paired_sig))


#And make a graphical representation:

# Not installed? 
# install.packages("ggvenn") # install via CRAN

library("ggvenn")

# gene names
ggvenn(list(DPN=rownames(res_DPN_paired_sig),
            OHT=rownames(res_OHT_paired_sig)),
       fill_color = c("red", "blue"))

#Not all pairwise comparisons are rendered in the default analysis. Therefore, if you are interested in one 
#comparison that hasn't been calculated, you just need to use the contrast argument. 
#For instance, to check differences between the two treatments DPN and OHT :
# calculating extra comparison
res_DPN_OHT_paired<-results(dsp,contrast=c("treatment","DPN","OHT"))

# checking data
head(res_DPN_OHT_paired[order(res_DPN_OHT_paired$padj),], 10)

#Checking the statistically significant genes:

# filtering based on 0.01 
res_DPN_OHT_paired_sig <- res_DPN_OHT_paired[ which(res_DPN_OHT_paired$padj < 0.01), ]

# checking data
res_DPN_OHT_paired_sig

#Generate an MA plot to visualize the log2 fold changes vs the mean normalized counts of all genes (dots in the plot), 
#where a red dot shows statistically significant genes:


plotMA(res_DPN_OHT_paired, 
       main="res_DPN_OHT_paired",
       cex=0.8, 
       ylim=c(-3,3))

#Modify the ylim parameter so you include all the values. 
#Hint: one way is to find the min/max values of the foldchanges: summary(res_DPN_OHT_paired$log2FoldChange)

#Make a list ordered by padj and then select the top 50 genes:

top=50

# ordering based on padj
res_DPN_OHT_paired_order <- res_DPN_OHT_paired[order(res_DPN_OHT_paired$padj),]

# selecting the top 50
res_DPN_OHT_paired_top <- res_DPN_OHT_paired_order[1:top,]

#checking data
res_DPN_OHT_paired_top

#Generate a heatmap of the transformed normalized counts:

# transforming data 
vsd <- vst(dsp, blind = FALSE)

# selecting data
res_DPN_OHT_paired_top2heatmap <- assay(vsd)[rownames(assay(vsd))%in%rownames(res_DPN_OHT_paired_top),]

# normalizing to the mean value
res_DPN_OHT_paired_top2heatmap_Mean <- res_DPN_OHT_paired_top2heatmap - rowMeans(res_DPN_OHT_paired_top2heatmap)

library(pheatmap)

# printing heatmap
pheatmap(res_DPN_OHT_paired_top2heatmap_Mean, 
         main         = "res_DPN_OHT_paired top 50 genes (padj)",
         cluster_cols = FALSE)

#Add some annotation for an easier interpretation
# Defining the annotation
annotation_col <- data.frame(
  Treatment = rep(c("Ctrl", "DPN", "OHT"),4),
  Patient = rep(c("A1", "A2", "A3","A4"),3)
)

# Defining the colors by group for plotting 
annotation_colors <- list(
  Treatment = c(Ctrl= "red", 
                DPN = "blue", 
                OHT = "green"),
  Patient = c(A1 = "yellow", 
              A2 = "orange", 
              A3 = "black", 
              A4 = "grey")
)

rownames(annotation_col) = samples$sample

pheatmap(res_DPN_OHT_paired_top2heatmap_Mean,
         main              = "res_DPN_OHT_paired top 50 genes (padj)",
         annotation_col    = annotation_col,
         annotation_colors = annotation_colors, 
         cluster_cols      = TRUE)



#Repeat the visualization using the top 20 genes instead.

topp=20

# ordering based on padj


# selecting the top 50
res_DPN_OHT_paired_top_20 <- res_DPN_OHT_paired_order[1:topp,]

#checking data
res_DPN_OHT_paired_top_20

# selecting data
res_DPN_OHT_paired_top20heatmap <- assay(vsd)[rownames(assay(vsd))%in%rownames(res_DPN_OHT_paired_top_20),]

# normalizing to the mean value
res_DPN_OHT_paired_top20heatmap_Mean <- res_DPN_OHT_paired_top20heatmap - rowMeans(res_DPN_OHT_paired_top20heatmap)


# printing heatmap
pheatmap(res_DPN_OHT_paired_top20heatmap_Mean, 
         main         = "res_DPN_OHT_paired top 20 genes (padj)",
         cluster_cols = FALSE)

#Add some annotation for an easier interpretation
# Defining the annotation
annotation_col <- data.frame(
  Treatment = rep(c("Ctrl", "DPN", "OHT"),4),
  Patient = rep(c("A1", "A2", "A3","A4"),3)
)

# Defining the colors by group for plotting 
annotation_colors <- list(
  Treatment = c(Ctrl= "red", 
                DPN = "blue", 
                OHT = "green"),
  Patient = c(A1 = "yellow", 
              A2 = "orange", 
              A3 = "black", 
              A4 = "grey")
)

rownames(annotation_col) = samples$sample

pheatmap(res_DPN_OHT_paired_top20heatmap_Mean,
         main              = "res_DPN_OHT_paired top 20 genes (padj)",
         annotation_col    = annotation_col,
         annotation_colors = annotation_colors, 
         cluster_cols      = TRUE)

#Visualize the counts of a specific gene grouped by treatment. In this case we are visualizing the gene with smallest padj:

# retrieving data with smallest padj
d<-plotCounts (dsp, 
               gene       = which.min(res_DPN_OHT_paired$padj),
               intgroup   = c("treatment","patient"), 
               returnData = TRUE)

# plotting counts from the different treatments
ggplot(d, aes(x=treatment, y=count, color=patient, shape=treatment)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))

  #Remember that you can save these results for processing outside R:

  write.table(res_DPN_OHT_paired, 
            file      = "res_DPN_OHT_paired.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            sep       = "\t", 
            quote     = FALSE,
            na        = "NA")

  #Functional analysis

#A common downstream analysis is to find pathways linked to our significant gene list. 
#Here we will perform an over-representation analysis using the package ReactomePA on the gene list where we compared DPN to OHT treatment.
#We start by importing the result table from DESeq2 (in case you do the analysis at another time):

# reading the data from a file
res <- read.csv("res_DPN_OHT_paired.txt", 
                header = T, 
                sep    = "\t")
head(res)


library(ReactomePA)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)

#We select the significant genes (padj < 0.05) and save them to new dataframe called res_sig:
# Selecting genes based on padjusted
res_sig<-res[ which(res$padj < 0.05), ]
nrow(res_sig) #117

#Reactome does not work with ENSEMBL ids, thus we need to convert them to ENTREZ before performing the pathway analysis:
# Mapping identifiers
res_sig$Entrez <- mapIds(org.Hs.eg.db, 
                         keys    = row.names(res_sig), 
                         column  = "ENTREZID", 
                         keytype = "ENSEMBL")
head(res_sig)


#We extract the fold changes of all the significant genes and we set the ENTREZ ids as name, 
#saving them to a variable called ressigFC. There will be genes that lack and ENTREZ id, these are annotatied as NA. 
#For the analysis to run properly, it is necessary to remove these genes.

# Selecting fold changes
ressigFC<-res_sig$log2FoldChange

# Adding the gene name to the list of fold changes
names(ressigFC)<-res_sig$Entrez

# Removing genes without an ENTREZ id
ressigFC<-ressigFC[!is.na(names(ressigFC))]
length(ressigFC)

#Then we run an over-representation analysis using enrichPathway(). This function uses a hypergeometric 
#model to test if the amount of genes from the input list associated to the pathway is larger then expected. Read more here.

# Running the test
x <- enrichPathway(gene         = names(ressigFC),
                   pvalueCutoff = 0.05, 
                   readable     = T)

# Looking at the results
dfx<-data.frame(x)
knitr::kable(dfx)



#Output columns are described as follows

#    ID - Reactome identifier
#    Description - Pathway description
#    GeneRatio - How many genes in our list hit this pathway of the entire input list
#    BgRatio - How many genes does the pathway contain of the total background
#    Pvalue - significant level
#    p.adj - Adjusted Pvalue, default Benjamini hochberg
#    qvalue - Adjusted Pvalue (FDR)
#    geneID - the genes from the input list that hits the pathway


# Saving the results
write.table(x         = dfx,
            file      = "DPN_OHT_ReactomePA.txt", 
            sep       = "\t",  
            quote     = FALSE, 
            col.names = NA)

#Let's plot the most significant pathways:
edox <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')

p1 <- cnetplot(x, color.params=list(foldChange=ressigFC))
p2 <- cnetplot(x, color.params=list(foldChange=ressigFC), circular = TRUE)

pdf("DPN_OHT_ReactomePA_Plots.pdf", height = 15, width= 15)
p1
p2
dev.off()