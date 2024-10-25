# The procedure described below will be applied in the same way to the data from the MDA and U2OS experiments.


# Using Ubuntu, we concatenate the forward reads from both lanes (L001 and L002) for each sample. We do the same 
# with the reverse reads. To do this, we use the following code (example for a single sample):

# cat EV_1_S13_L001_R1_001.fastq EV_1_S13_L002_R1_001.fastq> EV_1_S13_R1_combined.fastq.gz   (forward reads)
# cat EV_1_S13_L001_R2_001.fastq EV_1_S13_L002_R2_001.fastq> EV_1_S13_R2_combined.fastq.gz   (reverse reads)

# In RStudio, once the forward and reverse reads of both lanes have been concatenated, we create two vectors 
# to save the files with the forward and reverse reads (.fastq) of all the samples:

setwd("C:/Users/juand/OneDrive/Escritorio/TFM")

fastq_files_R1 <- list.files(path = "fastq/", pattern = "R1")
fastq_files_R2 <- list.files(path = "fastq/", pattern = "R2")

# We download the genomic annotation (Homo_sapiens.GRCh38.112.gtf file) and the reference genome 
# (Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa file) from the Ensembl database.

# From the downloaded reference genome we build an index of the reference genome, which is necessary to carry out the mapping 
# (or alignment) of the reads:

library(Rbowtie)

gunzip("genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz", remove=FALSE)
buildindex(basename="index",reference="genome/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa")

# To perform this mapping, we first create a vector for the .bam files that will be generated later. To do this, we replace 
# the ‘_R1.fastq.gz’ endings of the previously created vector with ‘.bam’, saving the result of this in the variable bam_files. 
# As what is preserved in the new vector is what precedes the mentioned ending and this is identical in the R1 and R2 files, it 
# is not necessary to carry out the substitution also with the names of the R2 files.

bam_files <- gsub(pattern = "_R1.fastq.gz",replacement = ".bam", fastq_files_R1)

# We incorporate the route into the generated vectors:

fastq_files_R1 <- paste0("fastq/",fastq_files_R1)
fastq_files_R2 <- paste0("fastq/",fastq_files_R2)
bam_files <- paste0("bam/", bam_files)




### ALIGNMENT

# To perform the alignment or mapping, we use the align function, arguing as input files the constructed index and the combined 
# forward and reverse reads (R1 and R2) of all samples from both experiments. This will generate output .bam files containing the 
# reads aligned to the reference genome and listed in bam_files:

library(Rsubread)

setwd("C:/Users/juand/OneDrive/Escritorio/TFM")
align(index="index",
      readfile1=fastq_files_R1, 
      readfile2 = fastq_files_R2, 
      TH1=2, type=0, 
      input_format="gzFASTQ", 
      output_format="BAM", 
      output_file=bam_files, 
      unique = TRUE, 
      nthreads = 20, 
      PE_orientation = "fr")


# We sort the reads from the .bam files by chromosome and coordinate:

for (bam in bam_files) {
  out <- suppressWarnings(sortBam(bam, "temporal"))
  file.rename(out, bam)
}

# We obtain the statistics of the alignment:

diagnostics <- list()
for (bam in bam_files) {
  total <- countBam(bam)$records
  mapped <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE)))$records
  marked <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE, isDuplicate=TRUE)))$records
  diagnostics[[bam]] <- c(Total=total, Mapped=mapped, Marked=marked)
}
diag.stats <- data.frame(do.call(rbind, diagnostics)) 
diag.stats$Prop.mapped <- diag.stats$Mapped/diag.stats$Total*100
diag.stats$Prop.marked <- diag.stats$Marked/diag.stats$Mapped*100
diag.stats

write.csv(diag.stats, "mapped_statistics.csv")

# We index the .bam files, generating .bambai files:

library(Rsamtools)

indexBam(bam_files)




### GENE EXPRESSION ANALYSIS

# We quantify the aligned reads for each gene in each sample by generating a single count table for both experiments. 
# For this, we use the aligned reads for each sample (.bam files) and the genomic annotation previously obtained:

filenames <- bam_files

library(Rsubread)

fc <- featureCounts(files = filenames,
                    annot.ext = "annotation/Homo_sapiens.GRCh38.112.gtf",
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="exon",
                    GTF.attrType="gene_id",
                    isPairedEnd = TRUE,
                    strandSpecific = 2,
                    nthreads = 20)
write.table(fc$counts, file = "results/raw_counts_RNAseq.tsv", sep = "\t")

seqdata <- read.delim("raw_counts_RNAseq.tsv")
seqdata
head(seqdata)
dim(seqdata)

# We proceed to divide the table of counts in which the gene expression of both experiments (MDA and U2OS) 
# is included into two tables, one for each experiment, thus continuing the analysis of each one independently:

seqdata_MDA.n <- (seqdata[,c(1,2,5,6,7,8)])
seqdata_U2OS <- (seqdata[,c(3,4,9,10)])

# We now filter the count tables to eliminate flat genes (genes not expressed or with low copy number) using as a measure the 
# parameter CPM (counts-per-million), i.e. the number of copies of a gene (loci) per million reads. We will select loci with 
# a number of reads greater than 0.5 CPM (=10-15 copies) in at least two of the ten samples in the study.

# To do this, we first calculate the number of CPMs for each gene using the cpm function:

library(limma)

library(edgeR)

myCPM_MDA <- cpm(seqdata_MDA)
myCPM_U2OS <- cpm(seqdata_U2OS)

head(myCPM_MDA.n)
head(myCPM_U2OS)

# We find out which values in the matrix of CPMs are greater than 0.5:

thresh_MDA <- myCPM_MDA > 0.5 
head(thresh_MDA)

thresh_U20S <- myCPM_U2OS > 0.5 
head(thresh_U20S)

# We make a table with the sum of the rows (genes) that meet the threshold (geneIDs with CPM greater than 0.5):

table(rowSums(thresh_MDA))

table(rowSums(thresh_U20S)) 

# We implement the filter, preserving genes with a CPM greater than 0.5 in at least two samples:

keep_MDA <- rowSums(thresh_MDA) >= 2
countdata.no.flat.MDA <- seqdata_MDA[keep_MDA,]

keep_U2OS <- rowSums(thresh_U20S) >= 2
countdata.no.flat.U2OS <- seqdata_U2OS[keep_U2OS,]

# In this way, we see which genes have been deleted (FALSE) and which are preserved (TRUE):
 
summary(keep_MDA)
summary(keep_U2OS)

# We note how the number of rows (genes) in the new matrix has been significantly reduced compared to the initial counting matrix, 
# as only non-flat genes have been selected:

dim(countdata.no.flat.MDA)
dim(countdata.no.flat.U2OS) 

dim(seqdata_MDA)
dim(seqdata_U2OS)




### MULTIDIMENSIONAL SCALING (MDS)

# We use the MDS (Multidimensional Scaling) technique that reduces the dimensionality of the data to graphically represent the 
# differences between the gene expression profiles of the samples in the form of distances between samples:

logcpm.MDA <- cpm(countdata.no.flat.MDA , log=TRUE)
plotMDS(logcpm.MDA)

logcpm.MDA.final <- logcpm.MDA[ , -c(3, 4)]

logcpm.U2OS <- cpm(countdata.no.flat.U2OS , log=TRUE)
plotMDS(logcpm.U2OS)

# We assign colours, shapes and a legend to the graphs using a ‘samples_info’ table with the experiment information broken down, 
# converting each column to factors and associating, finally, colour to condition and tip type to each replicate.

samples_info_MDA <- read.delim("samples_info_MDA.txt")
samples_info_MDA

samples_info_MDA$Condition <- as.factor(samples_info_MDA$Condition)
samples_info_MDA$Replicate <- as.factor(samples_info_MDA$Replicate)

sample.color.MDA <- c("grey", "red", "blue", "cyan")[samples_info_MDA$Condition]
sample.symbol.MDA <- c(16,17,18)[samples_info_MDA$Replicate]

samples_info_U2OS <- read.delim("samples_info_U2OS.txt")
samples_info_U2OS

samples_info_U2OS$Condition <- as.factor(samples_info_U2OS$Condition)
samples_info_U2OS$Replicate <- as.factor(samples_info_U2OS$Replicate)

sample.color.U2OS <- c("grey", "red", "blue", "cyan")[samples_info_U2OS$Condition]
sample.symbol.U2OS <- c(16,17,18)[samples_info_U2OS$Replicate]

# We generate .pdf files with the MDS representations:

pdf("results/MDS_MDA.pdf", width = 15, height = 15)
plotMDS(logcpm.MDA.final, col=sample.color, pch=sample.symbol, cex=2, cex.axis=2, cex.lab=2)
par(mar = c(15,15,15,15)) 
legend(1,1, xpd = T, fill=c("grey", "red", "blue", "cyan"),legend=levels(samples_info_MDA$Condition),cex=1,)
legend(1, 0.75, xpd = T, pch=c(16,17,18),legend=levels(samples_info_MDA$Replicate),cex=1)
dev.off()

pdf("results/MDS_U2OS.pdf", width = 15, height = 15)
plotMDS(logcpm.U2OS, col=sample.color, pch=sample.symbol, cex=2, cex.axis=2, cex.lab=2)
par(mar = c(15,15,15,15))
legend(1,1, xpd = T, fill=c("grey", "red", "blue", "cyan"),legend=levels(samples_info_U2OS$Condition),cex=1,)
legend(1, 0.75, xpd = T, pch=c(16,17,18),legend=levels(samples_info_U2OS$Replicate),cex=1)
dev.off()




### DIFFERENTIAL GENE EXPRESSION ANALYSIS

library(DESeq2)


## MDA

# We create a dataframe to store the replicate groups and their distribution (experimental design):

samples_MDA <- data.frame(conditions=c("Ctrl", "Ctrl", "MVO", "MVO", "shAEP", "shAEP"))
samples_MDA$conditions <- as.factor(samples_MDA$conditions)

# We create the dds object that combines the filtered quantification with the distribution of the replicas:

dds_MDA <- DESeqDataSetFromMatrix(countdata.no.flat.MDA,colData = samples_MDA, design=~conditions)
dds_MDA

# We use the DESeq function to normalise the data taking into account the different sequencing depth of each sample and the 
# variability between samples, so that these biases are eliminated and the samples are comparable to each other. This will 
# allow us to subsequently identify differentially expressed genes in the samples:

dds_MDA <- DESeq(dds_MDA)
dds_MDA

# To visualise the effect of the applied normalisation, we plotted the overall gene distribution before and after data 
# processing using boxplots:

# Before normalising ->

colnames(countdata.no.flat.MDA) <- c("Control_1","Control_2","MVO_1","MVO_2","shAEP_1","shAEP_2") 

png(filename = "results/boxplot_MDA_pre.png", width = 800, height = 800)

countdata.no.flat.MDA.final <- countdata.no.flat.MDA[ , -c(3, 4)]

boxplot(countdata.no.flat.MDA.final, outline=F,col=rainbow(4),ylab="Raw gene expression (MDA)",
        cex.lab=1.5)

dev.off()

# After normalising ->
  
normalized_counts_MDA <- counts(dds_MDA, normalized=TRUE)

colnames(normalized_counts_MDA) <- c("Control_1","Control_2","MVO_1","MVO_2","shAEP_1","shAEP_2") 

png(filename = "results/boxplot_MDA_post.png", width = 800, height = 800)

normalized_counts_MDA_final <- normalized_counts_MDA[ , -c(3, 4)]

boxplot(normalized_counts_MDA_final, outline=F,col=rainbow(4),ylab="Processed gene expression (MDA)",
        cex.lab=1.5)

dev.off()

# We use the results function to extract the results of the differential expression between the conditions we are interested in, 
# i.e. MVO vs Control and shAEP vs Control, thus obtaining a dataframe that includes parameters such as log2FoldChange, p-value 
# and adjusted p-value:

results_MVO.vs.Ctrl <- results(dds_MDA,contrast = c("conditions","MVO", "Ctrl"))

results_shAEP.vs.Ctrl <- results(dds_MDA,contrast = c("conditions","shAEP", "Ctrl"))

# We make some modifications to the generated dataframes and save them in .tsv format:

final_MVO.vs.Ctrl <- data.frame(
  GENEID = row.names(results_MVO.vs.Ctrl),
  log2BaseMean = log2(results_MVO.vs.Ctrl$baseMean),
  log2Ratio = results_MVO.vs.Ctrl$log2FoldChange,
  STDERR_log2Ratio = results_MVO.vs.Ctrl$lfcSE,
  pvalue = results_MVO.vs.Ctrl$pvalue,
  padjust = results_MVO.vs.Ctrl$padj
)

head(final_MVO.vs.Ctrl)

final_shAEP.vs.Ctrl <- data.frame(
  GENEID = row.names(results_shAEP.vs.Ctrl),
  log2BaseMean = log2(results_shAEP.vs.Ctrl$baseMean),
  log2Ratio = results_shAEP.vs.Ctrl$log2FoldChange,
  STDERR_log2Ratio = results_shAEP.vs.Ctrl$lfcSE,
  pvalue = results_shAEP.vs.Ctrl$pvalue,
  padjust = results_shAEP.vs.Ctrl$padj
)

head(final_shAEP.vs.Ctrl)


write.table (final_MVO.vs.Ctrl, file="results/final_MVO.vs.Ctrl.tsv", sep = "\t", row.names = FALSE)    

write.table (final_shAEP.vs.Ctrl, file="results/final_shAEP.vs.Ctrl.tsv", sep = "\t", row.names = FALSE)    



## U2OS

# We repeat the steps followed in the MDA experiment:

samples_U2OS <- data.frame(conditions=c("Ctrl", "Ctrl", "shCtsL", "shCtsL"))
samples_U2OS$conditions <- as.factor(samples_U2OS$conditions) 

dds_U2OS <- DESeqDataSetFromMatrix(countdata.no.flat.U2OS, colData = samples_U2OS, design=~conditions)
dds_U2OS


dds_U2OS <- DESeq(dds_U2OS)
dds_U2OS



colnames(countdata.no.flat.U2OS) <- c("Control_1","Control_2","shCtsL_1","shCtsL_2") 

png(filename = "results/boxplot_U2OS_pre.png", width = 800, height = 800)

boxplot(countdata.no.flat.U2OS, outline=F,col=rainbow(4),ylab="Raw gene expression (U2OS)",
        cex.lab=1.5)

dev.off()


normalized_counts_U2OS <- counts(dds_U2OS, normalized=TRUE)

colnames(normalized_counts_U2OS) <- c("Control_1","Control_2","shCtsL_1","shCtsL_2") 

png(filename = "results/boxplot_U2OS_post.png", width = 800, height = 800)

boxplot(normalized_counts_U2OS, outline=F,col=rainbow(4),ylab="Processed gene expression (U2OS)",
        cex.lab=1.5)

dev.off()




results_shCtsL.vs.Ctrl <- results(dds_U2OS, contrast = c("conditions", "shCtsL", "Ctrl"))

final_shCtsL.vs.Ctrl <- data.frame(
  GENEID = row.names(results_shCtsL.vs.Ctrl),
  log2BaseMean = log2(results_shCtsL.vs.Ctrl$baseMean),
  log2Ratio = results_shCtsL.vs.Ctrl$log2FoldChange,
  STDERR_log2Ratio = results_shCtsL.vs.Ctrl$lfcSE,
  pvalue = results_shCtsL.vs.Ctrl$pvalue,
  padjust = results_shCtsL.vs.Ctrl$padj
)

head(final_shCtsL.vs.Ctrl)

write.table (final_shCtsL.vs.Ctrl, file="results/final_shCtsL.vs.Ctrl.tsv", sep = "\t", row.names = FALSE)


# We now filter the sets of genes obtained for both experiments in excel, preserving only those with adjusted p-values less 
# than 0.01, thus ensuring that we select only significantly differentially expressed genes.

# We then split the resulting set of genes for each experiment into one whose genes have a log2ratio (log2foldchange) greater 
# than 1 and one whose genes have a log2ratio less than -1. This will allow us to obtain, on the one hand, the genes that are 
# activated in the shCtsL and shAEP samples (where the expression of the protease has been silenced) with respect to the control 
# and, on the other hand, those that are repressed in the same samples.

# We then generate Venn diagrams to verify, by graphical representation, the non-overlap between the figure representing 
# activated and repressed genes in each experiment, thus ensuring that genes are correctly categorised as activated or repressed:

library(VennDiagram)


# MDA

activated_genes_MDA <- read.table("results/final_shAEP.vs.Ctrl_filtered_activated.txt")
activated_genes_MDA <- activated_genes_MDA$V1

repressed_genes_MDA <- read.table("results/final_shAEP.vs.Ctrl_filtered_repressed.txt")
repressed_genes_MDA <- repressed_genes_MDA$V1

png(filename = "results/venn_diagram_MDA.png", width = 800, height = 800)

venn.plot.MDA <- venn.diagram(
  x = list(Set1 = activated_genes_MDA, Set2 = repressed_genes_MDA),
  category.names = c("Set 1", "Set 2"),
  fill = c("blue", "red"),
  alpha = 0.5, 
  cex = 2,    
  cat.cex = 2,  
  filename = NULL 
)

grid.draw(venn.plot.MDA)

dev.off()



# U2OS

activated_genes_U2OS <- read.table("results/final_shCtsL.vs.Ctrl_filtered_activated.txt")
activated_genes_U2OS <- activated_genes_U2OS$V1

repressed_genes_U2OS <- read.table("results/final_shCtsL.vs.Ctrl_filtered_repressed.txt")
repressed_genes_U2OS <- repressed_genes_U2OS$V1

png(filename = "results/venn_diagram_U2OS.png", width = 800, height = 800)

venn.plot.U2OS <- venn.diagram(
  x = list(Set1 = genes_activados_U2OS, Set2 = genes_reprimidos_U2OS),
  category.names = c("Set 1", "Set 2"),
  fill = c("blue", "red"),
  alpha = 0.5,  
  cex = 2,      
  cat.cex = 2,  
  filename = NULL  
)

grid.draw(venn.plot.U2OS)

dev.off()



# We also represent the differentially expressed genes in each experiment by volcanoplot. Here we can distinguish between 
# activated (green) and repressed (red) genes, as well as those not significantly differentially expressed. They are also 
# distributed according to their adjusted p-value.

library(ggplot2)


# MDA

volcano_plot_MDA <- ggplot(final_shAEP.vs.Ctrl, aes(x = log2Ratio, y = -log10(padjust))) +
  geom_point(aes(color = ifelse(padjust < 0.01 & log2Ratio > 1, "Activado", 
                                ifelse(padjust < 0.01 & log2Ratio < -1, "Reprimido", "No significativo")))) +
  scale_color_manual(values = c("Activado" = "green", "Reprimido" = "red", "No significativo" = "gray")) +
  labs(x = "log2 Fold Change", y = "-log10 Adjusted p-value", color = "Estado")

ggsave("results/volcano_plot_MDA.png", plot = volcano_plot_MDA)



#U2OS

volcano_plot_U2OS <- ggplot(final_shCtsL.vs.Ctrl, aes(x = log2Ratio, y = -log10(padjust))) +
  geom_point(aes(color = ifelse(padjust < 0.01 & log2Ratio > 1, "Activado", 
                                ifelse(padjust < 0.01 & log2Ratio < -1, "Reprimido", "No significativo")))) +
  scale_color_manual(values = c("Activado" = "green", "Reprimido" = "red", "No significativo" = "gray")) +
  labs(x = "log2 Fold Change", y = "-log10 Adjusted p-value", color = "Estado")

ggsave("results/volcano_plot_U2OS.png", plot = volcano_plot_U2OS)




### FUNCTIONAL ENRICHMENT

# In this phase of the analysis we elucidate the biological processes enriched in the sets of differentially expressed genes 
# for each experiment, i.e. we identify those processes in which genes found in these sets are involved. These biological 
# processes are identified as GO terms, following a structured vocabulary that describes such processes (gene ontology).

# To this end, on the one hand, the identification, as GO (Gene Ontology) terms, of biological processes overrepresented in 
# our gene sets was carried out using the online bioinformatics tool DAVID (Database for Annotation, Visualization, and 
# Integrated Discovery). On the other hand, the interactive online tool ShinyGO was used to, based on the list of GO terms 
# obtained, perform a deeper visual and functional analysis of the enriched GO terms in the gene set and identify enriched 
# metabolic and signalling pathways using databases such as KEGG.
