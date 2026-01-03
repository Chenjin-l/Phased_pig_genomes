
Inpath=path/to/F2-*featureCounts.txt

sed '1d' ${Inpath}/F2-503.Hap1.featureCounts.txt|cut -f1,6 > F2ind.quantile
pastestr='paste -d " " F2ind.quantile'
echo "Genename ratio"> biased.genes
for i in `ls ${Inpath}/F2-*.featureCounts.txt`;do
    id=`basename $i .featureCounts.txt`
    cut -f7 ${Inpath}/${id}.featureCounts.txt|sed '1,2d' >id.countall
    cut -f7 ${Inpath}/${id}.Hap1.featureCounts.txt|sed '1,2d'>id.countHap1
    cut -f7 ${Inpath}/${id}.Hap2.featureCounts.txt|sed '1,2d'>id.countHap2
    paste id.countall id.countHap1 id.countHap2|awk '{if(($2+$3)>5) {print $1,$2/($2+$3)*$1,$3/($2+$3)*$1} \
    else {print $1,$1/2,$1/2}}'|sed "1i${id} ${id}Hap1 ${id}Hap2" >${id}.count
    #Summary for Parents biased expression 
    cut -f1 F2ind.quantile |paste id.countall id.countHap1 id.countHap2 - |awk '{if(($2+$3)>5) print $4,($2+0.1)/($3+0.1)}' >> biased.genes
    rm id.countall id.countHap1 id.countHap2
    pastestr=$pastestr' '${id}'.count'
done

## =========================================================
## Generate R script for imprinting analysis
## =========================================================

OUT=Imprinting_analysis_phased_counts.R

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
