
INPUT_FILE=path/to/bfile        #paternal or maternal
PHENO_FILE=expression_gene.txt  #paternal or maternal gene_expression.txt
OUTPUT_DIR=path/to/heritability/out

# To obtain kinship for paternal or maternal
gcta64 --bfile $INPUT_FILE --make-grm --out GRM_PREFIX

num_cols=$(awk 'NR==1{print NF}' "$PHENO_FILE")

for ((col=3; col<=num_cols; col++)); do
    echo "正在分析表型列: $col"

    awk -v col="$col" '{print $1,$2,$col}' "$PHENO_FILE" > "$OUTPUT_DIR/col${col}_temp.pheno"

    gcta64 --reml --grm-gz "GRM_PREFIX" --pheno "$OUTPUT_DIR/col${col}_temp.pheno" --reml-pred-rand --out "$OUTPUT_DIR/col${col}"
# 清理临时文件
    rm "$OUTPUT_DIR/col${col}_temp.pheno"
done
