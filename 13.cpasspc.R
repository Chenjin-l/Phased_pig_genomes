#!/usr/bin/env Rscript
# =========================================================
# CPASSOC input preparation and cross-trait association
# =========================================================

library(data.table)
library(Matrix)
library(MASS)
library(future)

plan(sequential)               # force single-threaded execution
options(future.globals.maxSize = Inf)  # allow large global objects
# Load CPASSOC functions
source("/home/Mzhou/02.F2/z02data/14regnie/05cpasspc/FunctionSet.R")  

# ----------------------------
# 1. Define trait files and names
# ----------------------------
trait_files <- c(
  "/regnie/02_step2/diff/AntebrachialL/step2_diff_AntebrachialL.regenie",
  "/regnie/02_step2/diff/BrachialL/step2_diff_BrachialL.regenie",
  "/regnie/02_step2/diff/CarcassDiaL/step2_diff_CarcassDiaL.regenie",
  "/regnie/02_step2/diff/CarcassStrL/step2_diff_CarcassStrL.regenie",
  "/regnie/02_step2/diff/CruralL/step2_diff_CruralL.regenie",
  "/regnie/02_step2/diff/FemoralL/step2_diff_FemoralL.regenie",
  "/regnie/02_step2/diff/ScapularL/step2_diff_ScapularL.regenie"
)

trait_names <- c(
  "AntebrachialL","BrachialL","CarcassDiaL",
  "CarcassStrL","CruralL","FemoralL","ScapularL"
)

# ----------------------------
# 2. Read files and compute Z-scores
# ----------------------------
trait_list <- list()
for(i in seq_along(trait_files)){
  dt <- fread(trait_files[i])
  setnames(dt, trimws(names(dt)))  # remove extra spaces
  Z_col <- paste0("Z_", trait_names[i])
  dt[, (Z_col) := BETA / SE]      # Z = beta / SE
  trait_list[[i]] <- dt[, .(ID, CHROM, GENPOS, ALLELE0, ALLELE1, get(Z_col))]
  setnames(trait_list[[i]], "V6", Z_col)  # rename Z column
}

# ----------------------------
# 3. Merge traits by SNP
# ----------------------------
merged <- Reduce(function(x, y){
  merge(x, y, by = c("ID","CHROM","GENPOS","ALLELE0","ALLELE1"), all = TRUE)
}, trait_list)

# Save CPASSOC input
input_file <- "/home/Mzhou/02.F2/z02data/14regnie/05cpasspc/03data/CPASSOC_input_7traits.txt"
fwrite(merged, input_file, sep = "\t", quote = FALSE)
cat("✅ CPASSOC input saved to:", input_file, "\n")

# ----------------------------
# 4. Load input for CPASSOC
# ----------------------------
df <- fread(input_file)

# Extract Z-score columns
z_cols <- grep("^Z_", names(df), value = TRUE)
Z <- as.matrix(df[, ..z_cols])

# Sample sizes
SampleSize <- c(
  AntebrachialL = 589,
  BrachialL     = 589,
  CarcassDiaL   = 591,
  CarcassStrL   = 591,
  CruralL       = 587,
  FemoralL      = 589,
  ScapularL     = 589
)

# ----------------------------
# 5. Filter invalid SNPs
# ----------------------------
keep <- apply(Z, 1, function(x) all(is.finite(x)) && var(x) > 0)
Z <- Z[keep, ]
df <- df[keep, ]
cat("SNPs retained after filtering:", nrow(Z), "\n")

# ----------------------------
# 6. Correlation matrix
# ----------------------------
CorrMat <- cor(Z, use = "pairwise.complete.obs")
CorrMat <- CorrMat + diag(1e-8, ncol(CorrMat))  # avoid singularity

# ----------------------------
# 7. Gamma parameter estimation
# ----------------------------
cat("Estimating gamma parameters...\n")
para <- EstimateGamma(N = 1e4, SampleSize = SampleSize, CorrMatrix = CorrMat)

# ----------------------------
# 8. Compute SHet statistic
# ----------------------------
cat("Running SHet...\n")
SHet_stat <- SHet(X = Z, SampleSize = SampleSize, CorrMatrix = CorrMat)
p_SHet <- pgamma(q = SHet_stat - para[3], shape = para[1], scale = para[2], lower.tail = FALSE)

# ----------------------------
# 9. Compute SHom statistic
# ----------------------------
cat("Running SHom...\n")
SHom_stat <- SHom(X = Z, SampleSize = SampleSize, CorrMatrix = CorrMat)
p_SHom <- pchisq(SHom_stat, df = 1, lower.tail = FALSE)

# ----------------------------
# 10. Save CPASSOC output
# ----------------------------
df_out <- df
df_out$p_SHet <- p_SHet
df_out$p_SHom <- p_SHom

out_file <- "/home/Mzhou/02.F2/z02data/14regnie/05cpasspc/03data/CPASSOC_output_7traits.txt"
fwrite(df_out, out_file, sep = "\t")

