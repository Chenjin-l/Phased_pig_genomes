
#Compute LD between paternal and maternal eQTLs
# Paths 
OUT_DIR="/path/to/output_directory"
PLINK_BFILE="/path/to/plink_bfile"

# Files
INPUT="${OUT_DIR}/needld.snp"  # genes with independent paternal and maternal eQTLs
OUTPUT="${OUT_DIR}/ld_results.csv"
TMP="/tmp/ld_tmp"

mkdir -p "$TMP"
echo "phenotype_id,variant_id_p,variant_id_m,ld_r2" > "$OUTPUT"

while read pid snp_p snp_m; do
    [[ "$pid" == "phenotype_id" ]] && continue
    chr=${snp_p%%_*}; chr=${chr#chr}
    out="$TMP/${snp_p}_${snp_m}"

    plink --bfile "$PLINK_BFILE" --chr "$chr" \
          --ld "$snp_p" "$snp_m" --out "$out" >/dev/null 2>&1

    r2=$(grep "R-sq" "$out.log" | sed -E 's/.*= ([0-9.]+).*/\1/' || echo NA)
    echo "$pid,$snp_p,$snp_m,$r2" >> "$OUTPUT"
    rm -f "$out".{log,ld,nosex}
done < "$INPUT"

rm -rf "$TMP"
