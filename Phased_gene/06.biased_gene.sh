
## =========================================================
## Generate R script for phased_counts analysis
## =========================================================

OUT=phased_counts.R

cat << 'EOF' > ${OUT}
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

hap1_file <- "merged_Hap1_counts.txt"
hap2_file <- "merged_Hap2_counts.txt"
out_file  <- "Imprinting_test_with_category.txt"


hap1 <- read.table(hap1_file, header=TRUE, check.names=FALSE)
hap2 <- read.table(hap2_file, header=TRUE, check.names=FALSE)

genenames <- as.character(hap1[,1])

hap1 <- hap1 %>% select(-1) %>% mutate_all(as.numeric) %>% as.matrix() + 0.001
hap2 <- hap2 %>% select(-1) %>% mutate_all(as.numeric) %>% as.matrix() + 0.001

total_expr <- hap1 + hap2
low_expr <- total_expr < 5

paternal_ratio <- hap1 / total_expr
paternal_ratio[low_expr] <- NA

keep <- apply(
  paternal_ratio,
  1,
  function(x) sum(is.na(x)) < 0.9 * length(x)
)

paternal_ratio <- paternal_ratio[keep, ]
genenames <- genenames[keep]

gene_stat <- apply(paternal_ratio, 1, function(x) {

  x <- na.omit(x)

  if (length(x) < 58) {
    return(c(Nind=length(x), Mean=NA, T=NA, P=NA,
             Prop_pat=NA, Prop_mat=NA, Major=NA))
  }

  mean_ratio <- mean(x)

  if (sd(x) < 1e-10) {
    tval <- NA
    pval <- NA
  } else {
    tt <- t.test(x, mu=0.5)
    tval <- as.numeric(tt$statistic)
    pval <- as.numeric(tt$p.value)
  }

  prop_pat <- mean(x > 0.5)
  prop_mat <- mean(x < 0.5)
  major    <- max(prop_pat, prop_mat)

  return(c(
    Nind=length(x),
    Mean=mean_ratio,
    T=tval,
    P=pval,
    Prop_pat=prop_pat,
    Prop_mat=prop_mat,
    Major=major
  ))
})

gene_stat <- as.data.frame(t(gene_stat))
gene_stat$Genename <- genenames

gene_stat$Category <- "Consistently_balanced"

sig <- !is.na(gene_stat$P) & gene_stat$P < 1e-8

gene_stat$Category[sig & gene_stat$Major >= 0.8] <- "Stable_biased"
gene_stat$Category[sig & gene_stat$Major <  0.8] <- "Variably_biased"

gene_stat <- gene_stat %>%
  select(
    Genename,
    Nind,
    Mean,
    T,
    P,
    Prop_pat,
    Prop_mat,
    Major,
    Category
  )

write.table(
  gene_stat,
  file=out_file,
  sep="\\t",
  quote=FALSE,
  row.names=FALSE
)

cat("Done. Output:", out_file, "\\n")
EOF

Rscript Imprinting_analysis_phased_counts.R
