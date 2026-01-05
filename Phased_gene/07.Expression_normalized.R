####Expression_normalized.R
options(stringsAsFactors = FALSE)
library(edgeR)
#library(peer)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
#library(SNPRelate)
#----------------------------------------------------------------------------
### functions
"%&%" = function(a, b) { paste0(a, b) }
# Transform rows to a standard normal distribution
inverse_normal_transform = function(x) {
    qnorm(rank(x) / (length(x)+1))
}
#----------------------------------------------------------------------------
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
### data input <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ARGS <- commandArgs(trailingOnly = TRUE)
file_counts = ARGS[1] # Counts file. Row is gene, column is sample; rowname is gene id, colname is sample id
#file_tpm = ARGS[2] # TPM file. Row is gene, column is sample; rowname is gene id, colname is sample id
name=ARGS[2]
#tss_annot_file = ARGS[3] # TSS annotation file
#vcf.fn = ARGS[4] # Input data for genotype PCA. genotype data from imputation (VCF format)
#tis = ARGS[3] # Prefix of output file, like tissue name

if (!file.exists(file_counts)) { stop("Can not find the file_counts") }
#if (!file.exists(file_tpm)) { stop("Can not find the file_tpm") }
#if (!file.exists(tss_annot_file)) { stop("Can not find the tss_annot_file") }
#if (!file.exists(vcf.fn)) { stop("Can not find the vcf.fn") }

setwd("./")
#----------------------------------------------------------------------------
### main program
#----------------------------------------------------------------------------
# Input data for TMM calculation
Counts = read.table(file_counts,header=T,row.names=1,check.names=F)
## 1. prepare TMM
samids = colnames(Counts) # sample id
expr_counts = Counts
expr = DGEList(counts=expr_counts) # counts
nsamples = length(samids) # sample number
ngenes = nrow(expr_counts) # gene number

# calculate TMM
y = calcNormFactors(expr, method="TMM")
TMM = cpm(y,normalized.lib.sizes=T)

# expression thresholds
count_threshold = 10
#tpm_threshold = 0.1
sample_frac_threshold = 0.2
sample_count_threshold = 10

count_th = rowSums(expr_counts >= count_threshold)
ctrl2 = count_th >= (sample_frac_threshold * nsamples)
mask = ctrl2
TMM_pass = TMM[mask,] ##row is gene; column is sample

TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = inverse_normal_transform)) #apply to each row, each row represents one gene, observed values for all the samples
write.table(TMM_pass, paste0(name,".expr_tmm_aFC.txt"), sep = "\t", quote = FALSE)
write.table(TMM_inv, paste0(name,".expr_tmm_inv.txt"), sep = "\t", quote = FALSE)
