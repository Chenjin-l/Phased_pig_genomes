
line0=${PBS_ARRAYID}                       # line number in bamfiles list
OUTPUT=$2                      # output directory
ref=susScr11.fa                         # reference genome fasta
INPUT=${CHR}.shapiet.vcf.gz                       # directory with ${CHR}.shapiet.vcf.gz
THREADS=20               # threads

mkdir -p ${OUTPUT}
cd ${OUTPUT}

ID=$(sed -n "${line0}p" ${OUTPUT}/F2.JR.bamfiles | awk '{print $1}')
BAM=$(sed -n "${line0}p" ${OUTPUT}/F2.JR.bamfiles | awk '{print $2}')
ID2=$(echo $ID | sed 's/.$//' | fgrep -wf - ${OUTPUT}/F2.infor | awk '{print "F2-"$2}')

chroms=( $(awk '{print $1}' ${ref}.fai | sed -n 1,20p | tr '\n' ' ') )
echo "[$(date)] Chromosomes: ${chroms[*]}"

for CHR in "${chroms[@]}"; do
{
    echo "[$(date)] Processing ${CHR}"

    # Extract sample from phased VCF
    bcftools view -s ${ID2} ${INPUT}/${CHR}.shapiet.vcf.gz -Oz --threads=${THREADS} \
        -o ${ID2}.${CHR}.vcf.gz
    tabix -p vcf ${ID2}.${CHR}.vcf.gz

    # Add read group
    RG="@RG\\tID:${ID2}\\tPL:RNAseq\\tSM:${ID2}\\tLB:${ID2}\\tPU:1"

    # Extract per-chr BAM
    samtools view -h ${BAM} -@ ${THREADS} ${CHR} \
        | samtools addreplacerg -r ${RG} - -@${THREADS} \
        | samtools view -Sb -@${THREADS} \
        > ${ID}.${CHR}.bam
    samtools index ${ID}.${CHR}.bam -@${THREADS}

    # Haplotag with WhatsHap
    singularity exec ~/01.Biosoft/03.Sif/whatshap.simg \
        /opt/conda/envs/whatshap-env/bin/whatshap haplotag \
        --regions ${CHR} \
        -o ${ID}.${CHR}.phased.bam \
        -r ${ref} \
        ${ID2}.${CHR}.vcf.gz \
        ${ID}.${CHR}.bam
    samtools index ${ID}.${CHR}.phased.bam -@${THREADS}

} &
done && wait

ls ${ID}.*.phased.bam > ${ID}.bamfiles2
samtools cat -b ${ID}.bamfiles2 -o ${ID2}.RNAphased.bam
samtools index -@${THREADS} ${ID2}.RNAphased.bam
