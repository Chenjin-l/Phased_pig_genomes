library(dplyr)
library(purrr)

# --------------------------------------------------
# Define file path pattern (chromosomes 1–18)
# --------------------------------------------------
file_pattern <- paste0(1:18, "norma_norm.ibd")

# --------------------------------------------------
# Function to read and select required columns
# --------------------------------------------------
read_and_select <- function(file) {
  
  # Read input file (no header assumed)
  data <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  # Check column number and assign column names
  if (ncol(data) >= 4) {
    colnames(data)[1:4] <- c(
      "pair",
      "paternal_norm",
      "maternal_norm",
      "individual_norm"
    )
    
    # Return only the required columns
    return(data[, 1:5])
    
  } else {
    warning("File ", file, " has fewer than 5 columns and will be skipped.")
    return(NULL)
  }
}

# --------------------------------------------------
# Read and merge all chromosome files with error handling
# --------------------------------------------------
combined_data <- map_dfr(
  file_pattern,
  function(file) {
    tryCatch({
      df <- read_and_select(file)
      if (!is.null(df)) {
        df %>%
          mutate(chr = sub("norma_norm.ibd", "", basename(file)))
      }
    }, error = function(e) {
      message("Error processing file ", file, ": ", e$message)
      return(NULL)
    })
  }
)

# --------------------------------------------------
# Sanity check
# --------------------------------------------------
if (is.null(combined_data) || nrow(combined_data) == 0) {
  stop("No data were successfully read. Please check input file formats.")
}

# --------------------------------------------------
# Compute genome-wide mean normalized IBD values
# --------------------------------------------------
combined_data <- map_dfr(file_pattern, read_and_select) %>%
  group_by(pair) %>%
  summarise(
    mean_paternal   = sum(paternal_norm)   / 18,
    mean_maternal   = sum(maternal_norm)   / 18,
    mean_individual = sum(individual_norm) / 18
  )

# --------------------------------------------------
# Inspect results
# --------------------------------------------------
head(combined_data)

# --------------------------------------------------
# Save output
# --------------------------------------------------
write.table(
  combined_data,
  "mean_normalized_values.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
