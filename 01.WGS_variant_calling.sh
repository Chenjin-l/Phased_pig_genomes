# Input arguments
# ========================
ref=/path/to/Sscrofa11.1.fa
input=/path/to/input
output=/path/to/output
gatk=/path/to/gatk

ID=$1          # sample ID
CHR=$2         # chromosome or interval (e.g. chr1)

mkdir -p ${output}
cd ${output}

# ========================
# 1. Alignment (BWA-MEM)
# ========================
RG="@RG\\tID:${ID}\\tPL:ILLUMINA\\tSM:${ID}\\tLB:${ID}\\tPU:1"

bwa mem \
    -t 20 \
    -R "${RG}" \
    ${ref} \
    ${input}/${ID}_1.fastq.gz \
    ${input}/${ID}_2.fastq.gz \
| samtools view -@ 10 -hSb \
> ${ID}.bam

# ========================
# 2. Sort and index BAM
# ========================
samtools sort -@ 10 \
    -o ${ID}.sorted.bam \
    ${ID}.bam

samtools index ${ID}.sorted.bam
rm ${ID}.bam

# ========================
# 3. Remove PCR duplicates (Sambamba)
# ========================
sambamba markdup \
    -t 10 \
    --remove-duplicates \
    ${ID}.sorted.bam \
    ${ID}.sorted.markdup.bam

samtools index ${ID}.sorted.markdup.bam

# ========================
# 4. GATK HaplotypeCaller (GVCF)
# ========================
${gatk} HaplotypeCaller \
    --native-pair-hmm-threads 30 \
    -R ${ref} \
    -I ${ID}.sorted.markdup.bam \
    -ERC GVCF \
    -L ${CHR} \
    -O ${ID}.${CHR}.g.vcf.gz

tabix -p vcf ${ID}.${CHR}.g.vcf.gz


###############################################################################
# Steps below are cohort-level analyses
# Run AFTER all samples and chromosomes are completed
###############################################################################

# ========================
# 5. Joint genotyping (GLnexus)
# ========================
# Merge all per-sample GVCFs into a cohort VCF

glnexus_cli \
    --config gatk \
    *.g.vcf.gz \
> cohort.bcf

bcftools view cohort.bcf -Oz -o cohort.vcf.gz
tabix -p vcf cohort.vcf.gz

# ========================
# 6. Variant filtering
# ========================
# Filtering criteria:
#   i)   Read depth >= 3
#   ii)  Genotype quality >= 10
#   iii) Biallelic variants only
#   iv)  MAF > 0.05
#   v)   Missing genotype rate < 20%

bcftools view \
    -m2 -M2 -v snps,indels \
    -i 'INFO/DP>=3 && FORMAT/GQ>=10 && INFO/MAF>0.05 && F_MISSING<0.2' \
    cohort.vcf.gz \
    -Oz -o cohort.filtered.vcf.gz

tabix -p vcf cohort.filtered.vcf.gz
