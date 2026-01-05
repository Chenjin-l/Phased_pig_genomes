
PLINK_PREFIX="/your/path/to/genotypes"      # PLINK genotype prefix (bed/bim/fam)
EXPR_FILE="/your/path/to/expression.bed"    # Gene expression matrix

PREFIX="/your/path/to/output_prefix"        # Output prefix
SINGULARITY_IMAGE="/your/path/to/tensorqtl_GPU.simg"  # Singularity image
SEED=123                                    # Random seed for reproducibility

# Run cis/trans eQTL
singularity exec $SINGULARITY_IMAGE python3 -m tensorqtl \
    $PLINK_PREFIX \
    $EXPR_FILE \
    $PREFIX \
    --mode cis \  # cis_nominal or trans
    --seed $SEED

# Generate Python script to convert Parquet to CSV
cat << 'EOF' > parquet_to_csv.py
#!/usr/bin/env python3
# Convert TensorQTL Parquet output to CSV

import sys
import pandas as pd

if len(sys.argv) != 3:
    print("Usage: python parquet_to_csv.py <input.parquet> <output.csv>")
    sys.exit(1)

inp, out = sys.argv[1], sys.argv[2]

df = pd.read_parquet(inp)
df.to_csv(out, index=False)

print(f"Converted {inp} → {out}")
print(df.head())
EOF

python parquet_to_csv.py F2_cis_qtl.parquet result.txt.gz

