#!/usr/bin/env python3
Description:
This script classifies parental-origin effects (POE) for SNPs based on REGENIE Step 2 results.
For each SNP, maternal and paternal haplotype effects (Z-scores) are compared to assign POE type.

Input:
  - Diff, maternal, and paternal REGENIE output files per phenotype.

Output:
  - Classified POE CSV for each phenotype.
  - Detailed results including Z-scores and effect sizes.
"""

import os
import re
import pandas as pd

# ----------------------------
# Paths
# ----------------------------
BASE_DIR = "path/to/regnie/02_step2"
DIFF_DIR = os.path.join(BASE_DIR, "diff")
MAT_DIR  = os.path.join(BASE_DIR, "maternal")
PAT_DIR  = os.path.join(BASE_DIR, "paternal")
OUTPUT_BASE_DIR = os.path.join(BASE_DIR, "POE_results")

# ----------------------------
# Phenotype list
# ----------------------------
PHENOTYPES = [
    "AntebrachialL", "BrachialL", "CarcassDiaL",
    "CarcassStrL", "CruralL", "FemoralL", "ScapularL"
]

# ----------------------------
# Thresholds
# ----------------------------
SNP_TOTAL = 21273979
P_DIFF_THRESHOLD = 0.05 / SNP_TOTAL  # ≈ 2.35e-9

# ----------------------------
# Helper functions
# ----------------------------
def clean_id(x):
    """Clean SNP IDs by stripping spaces."""
    return re.sub(r'\s+', '', str(x).strip())

def classify_poe(row):
    """Classify POE type based on maternal/paternal Z-scores."""
    ZMat, ZPat = row['ZMat'], row['ZPat']
    if pd.isna(ZMat) or pd.isna(ZPat):
        return "Unclassified"

    # Bipolar effect: opposite signs, similar magnitude
    if ZMat * ZPat < 0 and abs(ZMat) > 0.5*abs(ZPat) and abs(ZPat) > 0.5*abs(ZMat):
        return "Bipolar"
    # Maternal dominant effect
    if abs(ZMat) >= 2*abs(ZPat):
        return "Maternal effect"
    # Maternal asymmetric
    if ZMat*ZPat > 0 and abs(ZMat) > abs(ZPat) and abs(ZPat) > 0.5*abs(ZMat):
        return "Maternal asymmetric"
    # Paternal dominant effect
    if abs(ZPat) >= 2*abs(ZMat):
        return "Paternal effect"
    # Paternal asymmetric
    if ZMat*ZPat > 0 and abs(ZPat) > abs(ZMat) and abs(ZMat) > 0.5*abs(ZPat):
        return "Paternal asymmetric"

    return "Unclassified"

# ----------------------------
# Main processing loop
# ----------------------------
for pheno in PHENOTYPES:
    print(f"\n=== Processing phenotype: {pheno} ===")

    # Input files
    diff_file = os.path.join(DIFF_DIR, pheno, f"step2_diff_{pheno}.regenie")
    mat_file  = os.path.join(MAT_DIR,  pheno, f"step2_mat_{pheno}.regenie")
    pat_file  = os.path.join(PAT_DIR,  pheno, f"step2_pat_{pheno}.regenie")

    # Output directory
    out_dir = os.path.join(OUTPUT_BASE_DIR, pheno)
    os.makedirs(out_dir, exist_ok=True)
    output_file  = os.path.join(out_dir, f"POE_classified_{pheno}.csv")
    detailed_file = os.path.join(out_dir, f"POE_detailed_{pheno}.csv")

    # Read files
    diff_df = pd.read_csv(diff_file, delim_whitespace=True)
    mat_df  = pd.read_csv(mat_file,  delim_whitespace=True)
    pat_df  = pd.read_csv(pat_file,  delim_whitespace=True)

    # Clean SNP IDs
    for df in [diff_df, mat_df, pat_df]:
        df['ID'] = df['ID'].apply(clean_id)

    # Compute P_diff and filter candidates
    diff_df['P_diff'] = 10 ** (-diff_df['LOG10P'])
    candidates = diff_df[diff_df['P_diff'] < P_DIFF_THRESHOLD].copy()
    if candidates.empty:
        print(f"No significant SNPs for {pheno}, skipping.")
        continue
    print(f"Candidate SNPs passing threshold: {len(candidates)}")

    # Merge maternal/paternal data
    mat_df_renamed = mat_df[['ID','BETA','SE']].rename(columns={'BETA':'BETA_mat','SE':'SE_mat'})
    pat_df_renamed = pat_df[['ID','BETA','SE']].rename(columns={'BETA':'BETA_pat','SE':'SE_pat'})
    candidates = candidates.merge(mat_df_renamed, on='ID', how='left')
    candidates = candidates.merge(pat_df_renamed, on='ID', how='left')

    # Remove missing values
    candidates.dropna(subset=['BETA_mat','SE_mat','BETA_pat','SE_pat'], inplace=True)
    print(f"Remaining SNPs with full maternal/paternal data: {len(candidates)}")

    # Compute Z-scores
    candidates['ZMat'] = candidates['BETA_mat'] / candidates['SE_mat']
    candidates['ZPat'] = candidates['BETA_pat'] / candidates['SE_pat']

    # Classify POE
    candidates['POE_type'] = candidates.apply(classify_poe, axis=1)

    # Display POE summary
    poe_counts = candidates['POE_type'].value_counts()
    print("POE classification summary:")
    for poe_type, count in poe_counts.items():
        print(f"  {poe_type}: {count} ({count/len(candidates)*100:.1f}%)")

    # Save outputs
    cols_to_save = ['CHROM','GENPOS','ID','ALLELE0','ALLELE1','P_diff',
                    'BETA_mat','SE_mat','ZMat','BETA_pat','SE_pat','ZPat','POE_type']
    candidates.to_csv(detailed_file, index=False, float_format='%.6g')
    candidates[cols_to_save].to_csv(output_file, index=False, float_format='%.6g')

    print(f"Output saved for {pheno}: {len(candidates)} SNPs -> {output_file}")
