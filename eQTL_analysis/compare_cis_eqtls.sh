#!/bin/bash
# Compare individual vs phased cis-eQTL

# Configurable paths
IND_FILE="/path/to/F2_incis.cis_qtl.txt.gz"   # individual eQTL result/paternal eQTL result
PHASED_FILE="/path/to/F2_pm_cis.cis_qtl.txt.gz" # phased eQTL result /maternal eQTL result
OUT_DIR="/path/to/output_directory"            # output directory

mkdir -p $OUT_DIR

# Filtered intermediate files
INDIVIDUAL_FILTER="$OUT_DIR/individual_fit.txt"
PHASED_FILTER="$OUT_DIR/phased_fit.txt"

# 1. Extract relevant columns and filter qval < 1e-5
# Columns: phenotype_id, variant_id, slope, slope_se, qval
# -----------------------------
# Individual eQTL
zcat $IND_FILE | awk 'NR==1{print "phenotype_id variant_id_individual slope_individual slope_se_individual qval_individual"; next}
{if($(NF-1)<1e-5) print $1,$7,$14,$15,$(NF-1)}' > $INDIVIDUAL_FILTER

# Phased eQTL
zcat $PHASED_FILE | awk 'NR==1{print "phenotype_id variant_id_pm slope_pm slope_se_pm qval_pm"; next}
{if($(NF-1)<1e-5) print $1,$7,$14,$15,$(NF-1)}' > $PHASED_FILTER

# 2. Sort files for join
# -----------------------------
sort -k1,1 $INDIVIDUAL_FILTER > $OUT_DIR/individual_sorted.txt
sort -k1,1 $PHASED_FILTER > $OUT_DIR/phased_sorted.txt

# 3. Merge by gene (phenotype_id) and classify
# -----------------------------
join -a1 -a2 -e "NA" -o 0,1.2,1.5,1.3,1.4,2.2,2.5,2.3,2.4 \
    $OUT_DIR/individual_sorted.txt $OUT_DIR/phased_sorted.txt | \
awk 'BEGIN{
    print "phenotype_id variant_id_individual qval_individual slope_indiv slope_se_indiv variant_id_pm qval_pm slope_pm slope_se_pm source"
}
{
    indiv_q=$3; phased_q=$7;
    source="NA";
    if($2=="NA") source="phased_only";
    else if($6=="NA") source="individual_only";
    else {source="both"; if(indiv_q<phased_q) source="both_indiv"; else source="both_phased"}
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,source
}' OFS='\t' > $OUT_DIR/eqtl_compare_final.txt

echo "Done! Output saved to $OUT_DIR/eqtl_compare_final.txt"
