
args=commandArgs(trailingOnly=TRUE)


#---------------------------------------------------
#Single-cell CN normalization (Roadmap DHS resize 2kb) and chromVAR analysis (R_chromVAR.Rmd)
#---------------------------------------------------

##Load annotation of Roadmap cell-type promoters and enhancers
DHS_annot_resize <- read.table(args[1], header=FALSE, sep ='\t')
DHS_annot_resize_new <- DHS_annot_resize[,5:ncol(DHS_annot_resize)]
DHS_annot_resize_Prom <- rowSums(DHS_annot_resize_new[,128:254])
DHS_annot_resize_Enh <- rowSums(DHS_annot_resize_new[,1:127])

DHS_matrix_resize <- read.table(args[2], header=FALSE, sep ='\t')

DHS_matrix_resize$chr <- "chr1"
DHS_matrix_resize[DHS_matrix_resize[,1]==1, 5] <- "chr1"
DHS_matrix_resize[DHS_matrix_resize[,1]==2, 5] <- "chr2"
DHS_matrix_resize[DHS_matrix_resize[,1]==3, 5] <- "chr3"
DHS_matrix_resize[DHS_matrix_resize[,1]==4, 5] <- "chr4"
DHS_matrix_resize[DHS_matrix_resize[,1]==5, 5] <- "chr5"
DHS_matrix_resize[DHS_matrix_resize[,1]==6, 5] <- "chr6"
DHS_matrix_resize[DHS_matrix_resize[,1]==7, 5] <- "chr7"
DHS_matrix_resize[DHS_matrix_resize[,1]==8, 5] <- "chr8"
DHS_matrix_resize[DHS_matrix_resize[,1]==9, 5] <- "chr9"
DHS_matrix_resize[DHS_matrix_resize[,1]==10, 5] <- "chr10"
DHS_matrix_resize[DHS_matrix_resize[,1]==11, 5] <- "chr11"
DHS_matrix_resize[DHS_matrix_resize[,1]==12, 5] <- "chr12"
DHS_matrix_resize[DHS_matrix_resize[,1]==13, 5] <- "chr13"
DHS_matrix_resize[DHS_matrix_resize[,1]==14, 5] <- "chr14"
DHS_matrix_resize[DHS_matrix_resize[,1]==15, 5] <- "chr15"
DHS_matrix_resize[DHS_matrix_resize[,1]==16, 5] <- "chr16"
DHS_matrix_resize[DHS_matrix_resize[,1]==17, 5] <- "chr17"
DHS_matrix_resize[DHS_matrix_resize[,1]==18, 5] <- "chr18"
DHS_matrix_resize[DHS_matrix_resize[,1]==19, 5] <- "chr19"
DHS_matrix_resize[DHS_matrix_resize[,1]==20, 5] <- "chr20"
DHS_matrix_resize[DHS_matrix_resize[,1]==21, 5] <- "chr21"
DHS_matrix_resize[DHS_matrix_resize[,1]==22, 5] <- "chr22"
DHS_matrix_resize[DHS_matrix_resize[,1]==23, 5] <- "chrX"
DHS_matrix_resize <- DHS_matrix_resize[,c(5,2:4)]
DHS_matrix_resize_all <- DHS_matrix_resize


##Load the CN tables for the DHS matrix

class_label<-read.table(args[5], sep = '\t', header=T, comment.char = "")
class_label<-class_label[,2]
class_label <- as.matrix(class_label)
class_label2 <- as.matrix(class_label)

#SC_CN_DHS <- read.table(args[3], header=FALSE, sep ='\t', comment.char = "")
##SC_CN_DHS==0 makes NA after normalization --> fix it with 2
#SC_CN_DHS2 <- SC_CN_DHS
#SC_CN_DHS2[SC_CN_DHS==0] <- 2

SC_CN_DHS <- matrix(0, nrow(DHS_matrix_resize_all), nrow(class_label))
for (i in 1:nrow(class_label)){
	CN_filename <- paste0("input_user/CN_matrix_for_DHS_", class_label[i], ".txt")
	tmp <- read.table(CN_filename, header=TRUE, sep ='\t', comment.char = "")
	SC_CN_DHS[,i] <- as.matrix(tmp)
}
SC_CN_DHS2 <- SC_CN_DHS
SC_CN_DHS2[SC_CN_DHS==0] <- 2


##Load the count tables for the DHS matrix

data1_resize<-read.table(args[4], sep = '\t', header=T, comment.char = "")
data1_new_resize<-data1_resize[,4:ncol(data1_resize)];

##Sort the order of single-cells
subclone_list <- read.table("input_user/input_subclonality.txt", header=T, sep ='\t', comment.char = "")

Deeptool_result_name <- as.data.frame(as.matrix(colnames(data1_new_resize)))
Deeptool_result_name$index <- 0
for (j in 1:nrow(Deeptool_result_name)){
	Deeptool_result_name[j,1] <- strsplit(Deeptool_result_name[j,1], ".sort.mdup.sc_pre_mono_sort_for_mark_uniq.bam")[[1]][1]
	Deeptool_result_name[j,2] <- which(subclone_list[,1]==Deeptool_result_name[j,1])
}
data1_new_resize <- data1_new_resize[,order(Deeptool_result_name[,2])]


data_bulk_resize<-cbind(data1_new_resize*2/SC_CN_DHS2)



DHS_count_resize <- data_bulk_resize[DHS_matrix_resize_all[,1]!="chrY" & DHS_matrix_resize_all[,1]!="chrM" & DHS_annot_resize_Enh>0,]
DHS_matrix_resize <- DHS_matrix_resize_all[DHS_matrix_resize_all[,1]!="chrY" & DHS_matrix_resize_all[,1]!="chrM" & DHS_annot_resize_Enh>0,]

#> sum(DHS_matrix_resize_all[,1]!="chrY" & DHS_matrix_resize_all[,1]!="chrM" & DHS_annot_resize_Enh>0)
#[1] 619368
#> sum(DHS_matrix_resize_all[,1]!="chrY" & DHS_matrix_resize_all[,1]!="chrM" & DHS_annot_resize_Prom>0)
#[1] 45963




#class_label<-read.table(args[5], sep = '\t', header=T, comment.char = "")
#class_label<-class_label[,2]
#class_label <- as.matrix(class_label)
#class_label2 <- as.matrix(class_label)



CM<-DHS_count_resize #Input: count table
CM_row<-DHS_matrix_resize[,4]
rownames(CM)<-as.matrix(CM_row)



##---------------------------------------------------------------------
##STEP1: count---------------------------------------------------------
##Loading the package--------------------------------------------------

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)



write.table(DHS_matrix_resize, "result/peak_roadmap_resize_2kb_Enh.bed", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
peakfile <- "result/peak_roadmap_resize_2kb_Enh.bed"
peaks <- getPeaks(peakfile, sort_peaks = TRUE)

#Warning message:
#In getPeaks(peakfile, sort_peaks = TRUE) :
#  Peaks are overlapping!After getting peak counts, peaks can be reduced to non-overlapping set using filterPeaks function

##Extract sorted peaks
df <- data.frame(seqnames=seqnames(peaks),
  starts=start(peaks)-1,
  ends=end(peaks),
  names=c(1:length(peaks))
  )

##Getting sorted peak count

CM_sort <- CM

colnames(CM_sort) <- colnames(CM)
my_counts_matrix <- as.matrix(CM_sort)

fragment_counts <- SummarizedExperiment(assays = list(counts = my_counts_matrix),rowRanges = peaks)
fragment_counts@colData$depth <- colSums(CM)
#fragment_counts@colData$Cell_Type <- conds

##Adding GC content
library(BSgenome.Hsapiens.UCSC.hg38)
fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
head(rowData(fragment_counts))


##---------------------------------------------------------------------
##STEP2: Annotations---------------------------------------------------
##Loading the package--------------------------------------------------
register(SerialParam())
jaspar_motifs <- getJasparMotifs() # default species is human


##Finding motif matches--------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)

# First get filtered counts
counts_filtered <- filterSamples(fragment_counts, min_depth = 1500, min_in_peaks = 0.15, shiny = FALSE)
counts_filtered <- filterPeaks(counts_filtered)



# get motif matches
motif_ix <- matchMotifs(jaspar_motifs, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)



# computing deviations
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)
dev_zscore <- deviationScores(dev)
write.table(dev_zscore, args[6], sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)
save(dev, file=args[7])

## computing variability of each motif
variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE)
write.table(variability, args[8], sep="\t", quote=FALSE,row.names=TRUE, col.names=TRUE)

