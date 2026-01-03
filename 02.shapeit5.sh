#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# SHAPEIT5 Pedigree-aware phasing pipeline
#
# Inputs:
#   - Filtered VCF (bgzipped + indexed)
#   - Pedigree file (F2.pedigree.fam)
#   - Chromosome list (optional, defaults to chr1..chr18, X, Y, M)
#
# Outputs:
#   - phased VCFs per chromosome
###############################################################################

THREADS=50
VCF=$1
PED=$2
REF=F2.pedigree.fam.revise
OUTDIR=$4

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# Define chromosomes if not provided
chroms=($(seq -f "chr%g" 1 18) chrX chrY chrM)

for chr in "${chroms[@]}"; do
    echo "[$(date)] Processing ${chr}"

    # 1. Filter VCF: missing rate < 0.2, MAF > 0.01, biallelic, minGQ 10, minDP 2
    vcftools --gzvcf ${VCF} \
        --chr ${chr} \
        --max-missing 0.8 \
        --maf 0.01 \
        --min-alleles 2 --max-alleles 2 \
        --minGQ 10 --minDP 2 \
        --recode --recode-INFO-all \
        --stdout | bgzip -c > ${chr}.filtered.vcf.gz

    tabix -p vcf ${chr}.filtered.vcf.gz

    # 2. Reheader (optional, if sample names need correction)
    # bcftools reheader -s rename.IDS2 ${chr}.filtered.vcf.gz -o ${chr}.rename.vcf.gz

    # 3. Fill AC/AN tags
    bcftools +fill-tags ${chr}.filtered.vcf.gz -- -t AC,AN | bgzip -c > ${chr}.addAC.vcf.gz
    tabix -p vcf ${chr}.addAC.vcf.gz

    # 4. SHAPEIT5 pedigree-aware phasing
    singularity exec shapeit5.simg phase_common_static \
        --input ${chr}.addAC.vcf.gz \
        --pedigree ${PED} \
        --region ${chr} \
        --thread ${THREADS} \
        --log ${chr}.shapeit.log \
        --output ${chr}.phased.bcf

    # 5. Convert to VCF
    bcftools view ${chr}.phased.bcf | bgzip -c > ${chr}.phased.vcf.gz
    tabix -p vcf ${chr}.phased.vcf.gz
done

echo "[$(date)] SHAPEIT5 phasing completed for all chromosomes."
