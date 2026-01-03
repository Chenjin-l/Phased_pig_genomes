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


