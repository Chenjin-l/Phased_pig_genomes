#!/usr/bin/env bash
set -euo pipefail


# --------------------------
# Parameters
# --------------------------
OUTPUT=/path/to/output
INPUT=/path/to/phased_BAMs
GTF=/home/jxlabgdp/02.GenomeIndex/03.Sus11.1/04.Muscle-GTF/Sus_scrofa.Sscrofa11.1.100.revise.gtf
THREADS=40

# Select BAM list file
BAMLIST=$1  # Pass as argument: F2.JR.Phased.bamfiles / Hap1Hap2.bamfiles / Hap1Hap2.bamfiles2
echo "Processing BAM list: ${BAMLIST}"

# --------------------------
# Loop over all samples in the list
# --------------------------
while read -r ID BAM; do
    echo "[$(date)] Processing sample: ${ID}"
    echo "Input BAM: ${BAM}"

    # --------------------------
    # Optional: FPKM / TPM quantification with StringTie
    # --------------------------
    # Uncomment if normalized expression is needed
    # singularity exec -B ${INPUT}:${INPUT} /home/guilu/software/139.stringtie2.sif stringtie \
    #     -p ${THREADS} \
    #     -e -B \
    #     -G ${GTF} \
    #     -o ${OUTPUT}/${ID}.gtf \
    #     -A ${OUTPUT}/${ID}.tsv \
    #     ${INPUT}/${BAM} \
    #     --fr

    # --------------------------
    # Gene-level raw counts with featureCounts
    # --------------------------
    singularity exec -B ${INPUT}:${INPUT} /home/guilu/software/star_2.7.9a_subread.simg featureCounts \
        -T ${THREADS} \
        -p \
        -t gene \
        -g gene_id \
        -a ${GTF} \
        -o ${OUTPUT}/${ID}.featureCounts.txt \
        -s 2 \
        ${INPUT}/${BAM}

    echo "[$(date)] Finished sample: ${ID}"
done < ${BAMLIST}

echo "All samples in ${BAMLIST} processed."
