#!/usr/bin/env Rscript
# Haplotype-specific genetic effect estimation
# This script computes haplotype-specific effects for SNP-gene pairs
# using allele-resolved expression data.

library(dplyr)

# -------- USER SETTINGS --------
absolute_paths_file <- "absolute_paths.list"  # File containing full paths to .snpeffect files
## *.snpeffect files
FORMAT  GT      Geneid  ENSSSCG00000000058
F2-1000Hap1     0|0     F2-1000Hap1     0.760757019581693
F2-1000Hap2     0|0     F2-1000Hap2     -1.23464019689639
F2-1001Hap1     0|0     F2-1001Hap1     1.02872517973089
F2-1001Hap2     1|1     F2-1001Hap2     -1.0877520893075

output_dir <- "./haplotype_results"           # Output directory for results
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# -------------------------------

# Read the list of files to process
absolute_paths <- readLines(absolute_paths_file)
cat(sprintf("Found %d .snpeffect files to process\n", length(absolute_paths)))

# ---------- FUNCTIONS ----------

# Compute haplotype-specific Z-test
calc_z_test <- function(group1, group2) {
  if (length(group1) < 2 || length(group2) < 2) {
    # Not enough samples to compute Z-test
    return(list(z = NA, p = NA, effect = NA, se = NA,
                n1 = length(group1), n2 = length(group2)))
  }
  # Compute effect size (mean difference)
  effect <- mean(group1) - mean(group2)
  # Standard error of effect
  se <- sqrt((sd(group1)^2 / length(group1)) + (sd(group2)^2 / length(group2)))
  # Z-statistic
  z_score <- effect / se
  # Two-sided p-value
  p_value <- 2 * pnorm(-abs(z_score))
  list(z = z_score, p = p_value, effect = effect, se = se,
       n1 = length(group1), n2 = length(group2))
}

# Parse gene, SNP, and allele information from file path
parse_info_from_path <- function(file_path) {
  file_name <- basename(file_path)
  file_name_no_ext <- sub("\\.snpeffect$", "", file_name)
  parts <- strsplit(file_name_no_ext, "_")[[1]]

  if (length(parts) >= 5) {
    # Standard format: chr1_123456_A_T_ENSG00000000001
    gene_id <- parts[length(parts)]
    snp_id <- paste(parts[1:(length(parts)-1)], collapse = "_")
    allele0 <- parts[length(parts)-2]
    allele1 <- parts[length(parts)-1]
  } else {
    # Non-standard or missing gene ID
    gene_id <- "UNKNOWN"
    snp_id <- file_name_no_ext
    allele0 <- ifelse(length(parts)>=2, parts[length(parts)-1], "NA")
    allele1 <- ifelse(length(parts)>=1, parts[length(parts)], "NA")
  }

  list(gene_id = gene_id, snp_id = snp_id, allele0 = allele0,
       allele1 = allele1, file_path = file_path)
}

# Process a single .snpeffect file
process_single_file <- function(file_path) {
  info <- parse_info_from_path(file_path)

  # Try reading the file
  data <- tryCatch({
    read.table(file_path, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, fill = TRUE,
               strip.white = TRUE, blank.lines.skip = TRUE,
               comment.char = "")
  }, error = function(e) return(NULL))
  if (is.null(data) || nrow(data) == 0 || ncol(data) < 4) return(NULL)

  # Prepare data frame with haplotype and genotype group
  df <- data %>%
    select(Sample = 1, Genotype = 2, GeneID = 3, Expression = 4) %>%
    mutate(
      Haplotype = ifelse(grepl("Hap1", Sample), "Hap1", "Hap2"),
      Genotype_Group = ifelse(substr(Genotype,1,1) == "0", "0-type", "1-type")
    ) %>%
    filter(!is.na(Expression))

  # Extract expression values by haplotype and genotype
  hap1_0 <- df %>% filter(Haplotype=="Hap1", Genotype_Group=="0-type") %>% pull(Expression)
  hap1_1 <- df %>% filter(Haplotype=="Hap1", Genotype_Group=="1-type") %>% pull(Expression)
  hap2_0 <- df %>% filter(Haplotype=="Hap2", Genotype_Group=="0-type") %>% pull(Expression)
  hap2_1 <- df %>% filter(Haplotype=="Hap2", Genotype_Group=="1-type") %>% pull(Expression)

  # Compute haplotype-specific effects and Z-tests
  hap1_res <- calc_z_test(hap1_0, hap1_1)
  hap2_res <- calc_z_test(hap2_0, hap2_1)

  # Compute difference between Hap1 and Hap2 effects
  if (!is.na(hap1_res$se) & !is.na(hap2_res$se)) {
    eff_diff <- hap1_res$effect - hap2_res$effect
    se_diff <- sqrt(hap1_res$se^2 + hap2_res$se^2)
    z_diff <- eff_diff / se_diff
    p_diff <- 2 * pnorm(-abs(z_diff))
  } else {
    eff_diff <- se_diff <- z_diff <- p_diff <- NA
  }

  # Return result as a data frame
  data.frame(
    gene = info$gene_id, snp = info$snp_id, ALLELE0 = info$allele0, ALLELE1 = info$allele1,
    Hap1_Z = round(hap1_res$z,4), p_Hap1_Z = format(hap1_res$p, scientific = TRUE, digits = 4),
    Hap2_Z = round(hap2_res$z,4), p_Hap2_Z = format(hap2_res$p, scientific = TRUE, digits = 4),
    Zdiff = round(z_diff,4), p_Zdiff = format(p_diff, scientific = TRUE, digits = 4),
    Hap1_effect = round(hap1_res$effect,4), Hap2_effect = round(hap2_res$effect,4),
    Effect_diff = round(eff_diff,4),
    Hap1_n = paste(hap1_res$n1,hap1_res$n2,sep="/"),
    Hap2_n = paste(hap2_res$n1,hap2_res$n2,sep="/"),
    stringsAsFactors = FALSE
  )
}

# ---------- MAIN SCRIPT ----------

# Process all files and combine results
results <- lapply(absolute_paths, function(f) {
  cat("Processing:", basename(f), "\n")
  process_single_file(f)
})
results <- do.call(rbind, results)

# Save full results
write.csv(results, file.path(output_dir,"haplotype_comparison_full.csv"),
          row.names = FALSE, quote = FALSE)

# Save simplified results
simple_res <- results %>% select(gene, snp, ALLELE0, ALLELE1,
                                 Hap1_Z, p_Hap1_Z, Hap2_Z, p_Hap2_Z,
                                 Zdiff, p_Zdiff, Hap1_n, Hap2_n)
write.csv(simple_res, file.path(output_dir,"haplotype_comparison_simple.csv"),
          row.names = FALSE, quote = FALSE)

cat("Done! Results saved in", output_dir, "\n")
