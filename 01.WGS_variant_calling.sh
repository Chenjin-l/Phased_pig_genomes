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

# 1. BWA-MEM alignment
# --------------------------
RG="@RG\\tID:${SAMPLE}\\tPL:ILLUMINA\\tSM:${SAMPLE}\\tLB:${SAMPLE}\\tPU:1"

bwa mem -t 20 -R "${RG}" \
    ${ref} \
    ${input}/${SAMPLE}_1.fastq.gz \
    ${input}/${SAMPLE}_2.fastq.gz \
| samtools view -@ 10 -hSb > ${SAMPLE}.bam

# --------------------------
# 2. Sort and index BAM
# --------------------------
samtools sort -@ 10 -o ${SAMPLE}.sorted.bam ${SAMPLE}.bam
samtools index ${SAMPLE}.sorted.bam
rm ${SAMPLE}.bam

# --------------------------
# 3. Remove PCR duplicates (Sambamba)
# --------------------------
sambamba markdup -t 10 --remove-duplicates \
    ${SAMPLE}.sorted.bam ${SAMPLE}.sorted.markdup.bam
samtools index ${SAMPLE}.sorted.markdup.bam

# --------------------------
# 4. GATK HaplotypeCaller (GVCF)
# --------------------------
${gatk} HaplotypeCaller \
    --native-pair-hmm-threads 30 \
    -R ${ref} \
    -I ${SAMPLE}.sorted.markdup.bam \
    -ERC GVCF \
    -L ${CHR} \
    -O ${SAMPLE}.${CHR}.g.vcf.gz

tabix -p vcf ${SAMPLE}.${CHR}.g.vcf.gz


# --------------------------
# 5.GLnexus joint genotyping, filtering
# --------------------------
chroms=$(seq -f "chr%g" 1 18) 

for CHR in "${chroms[@]}"; do
    echo "[$(date)] Processing ${CHR}"

    # --------------------------
    # 1. Prepare per-chr GVCF list
    # --------------------------
    ls ${input}/*/*${CHR}.*g.vcf.gz > ${CHR}.list

    # --------------------------
    # 2. Joint genotyping with GLnexus
    # --------------------------
    singularity exec -B ${input}:${input} /home/guilu/software/glnexus.sif \
        glnexus_cli --config gatk --list ${CHR}.list --dir ${CHR}.GLnexus.DB

    bcftools view ${CHR}.GLnexus.DB/cohort.bcf | bgzip -@ ${THREADS} -c > ${CHR}.cohort.vcf.gz
    tabix -p vcf ${CHR}.cohort.vcf.gz

    # --------------------------
    # 3. Per-chr filtering & rename
    # --------------------------
    singularity exec /home/guilu/software/087.vcftools.sif \
        vcftools --gzvcf ${CHR}.cohort.vcf.gz \
                 --max-missing 0.8 --maf 0.01 --min-alleles 2 --max-alleles 2 \
                 --recode --recode-INFO-all --stdout | bgzip -c > ${CHR}.filtered.vcf.gz
    tabix -p vcf ${CHR}.filtered.vcf.gz

    bcftools reheader -s rename.IDS2 ${CHR}.filtered.vcf.gz -o ${CHR}.rename.vcf.gz
    bcftools +fill-tags ${CHR}.rename.vcf.gz -- -t AC,AN | bgzip -c > ${CHR}.addAC.vcf.gz
    tabix -p vcf ${CHR}.addAC.vcf.gz
done
