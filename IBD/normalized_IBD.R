library(dplyr)

# --------------------------------------------------
# 1. Read chromosome length file
# --------------------------------------------------
# Input file format: Chr  Length  #chr_length.txt from run_KING_IBD.sh
chr_lengths <- read.table(
  "chr_length.txt",
  header = TRUE,
  sep = " ",
  col.names = c("Chr", "Length")
)

# Compute the average length for each chromosome
# (in case multiple records exist per chromosome)
chr_avg_length <- chr_lengths %>%
  group_by(Chr) %>%
  summarise(avg_length = mean(Length)) %>%
  arrange(Chr)  # Ensure chromosomes are ordered from 1 to 18

# --------------------------------------------------
# 2. Process IBD files for chromosomes 1–18
# --------------------------------------------------
for (chr in 1:18) {
  
  # Input IBD file for the current chromosome
  input_file <- paste0(chr, "norma.ibd")
  
  # Output file with length-normalized IBD values
  output_file <- paste0(chr, "norma_norm.ibd")
  
  # Read IBD data
  # Columns represent:
  # pair: individual pair ID
  # i: index or marker ID
  # paternal: paternal haplotype IBD length
  # maternal: maternal haplotype IBD length
  # individual: individual-level IBD length
  data <- read.table(input_file, header = FALSE, sep = " ")
  colnames(data) <- c("pair", "paternal", "maternal", "individual")
  
  # Retrieve the average length of the current chromosome
  current_length <- chr_avg_length$avg_length[chr_avg_length$Chr == chr]
  
  # --------------------------------------------------
  # Normalize IBD values by chromosome length
  # --------------------------------------------------
  data$paternal_norm   <- data$paternal   / current_length
  data$maternal_norm   <- data$maternal   / current_length
  data$individual_norm <- data$individual / current_length
  
  # Compute log10 ratio of paternal to maternal normalized IBD
  # A small constant is added to avoid division by zero
  data$log_ratio_norm <- log10(
    (data$paternal_norm + 1e-10) /
    (data$maternal_norm + 1e-10)
  )
  
  # --------------------------------------------------
  # Save normalized results
  # --------------------------------------------------
  # The output file contains both original and normalized columns
  write.table(
    data,
    file = output_file,
    sep = " ",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  cat("Processed:", input_file, "-> Saved as", output_file, "\n")
}
