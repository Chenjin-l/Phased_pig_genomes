

#!/usr/bin/env python3
# Usage: python parquet_to_csv.py F2_cis_qtl.parquet result.txt.gz
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


