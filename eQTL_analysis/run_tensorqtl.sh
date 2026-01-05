
PLINK_PREFIX="/your/path/to/genotypes"      # PLINK genotype prefix (bed/bim/fam)
EXPR_FILE="/your/path/to/expression.bed"    # Gene expression matrix  ration(paternal_counts/all_counts);indivual-level;phased-level;POE(paternal,maternal)

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



