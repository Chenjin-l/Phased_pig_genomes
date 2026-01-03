#!/bin/bash

BFILE="path/to/bfile"         #diff/paternal/maternal bfile
PHENO_FILE="path/to/bfile/bone.txt"
COVAR_FILE="path/to/bfile/sex_batch.txt"
COVARS="Sex,Batch"
QT_FLAG="--qt --apply-rint"    # Quantitative trait with rank-based inverse normal transformation
THREADS=12                     # Number of threads
OUT_DIR="path/to/regnie/out"
PHENO="$1"                     # Phenotype column name (pass as argument)

mkdir -p $OUT_DIR/00_step1_prune
mkdir -p $OUT_DIR/01_step1
mkdir -p $OUT_DIR/02_step2/diff/$PHENO


# Step 0: LD pruning to extract independent SNPs
plink \
  --bfile $BFILE \
  --geno 0.05 \
  --indep-pairwise 1000 50 0.2 \
  --maf 0.05 \
  --out $OUT_DIR/00_step1_prune/step1_prune \
  --threads 8

plink \
  --bfile $BFILE \
  --extract $OUT_DIR/00_step1_prune/step1_prune.prune.in \
  --make-bed \
  --out $OUT_DIR/00_step1_prune/fit_pheno_step1 \
  --threads 8

# Step 1: Running REGENIE Step 1
regenie \
  --step 1 \
  --bed $OUT_DIR/00_step1_prune/fit_pheno_step1 \
  --phenoFile $PHENO_FILE \
  --covarFile $COVAR_FILE \
  --covarCol $COVARS \
  --catCovarList $COVARS \
  $QT_FLAG \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix $OUT_DIR/01_step1/step1_tmp \
  --threads $THREADS \
  --out $OUT_DIR/01_step1/step1_bone

# Step 2: Running REGENIE Step 2
regenie \
  --step 2 \
  --bed /home/Mzhou/02.F2/z02data/02vcf/01phenovcf/07impvcf/imp \
  --covarFile $COVAR_FILE \
  --covarCol $COVARS \
  --catCovarList $COVARS \
  --phenoFile $PHENO_FILE \
  --phenoCol $PHENO \
  --pred $OUT_DIR/01_step1/step1_bone_pred.list \
  --bsize 200 \
  --threads $THREADS \
  --out $OUT_DIR/02_step2/diff/$PHENO/step2_diff

