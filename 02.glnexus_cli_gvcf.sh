
# --------------------------
# GLnexus joint genotyping, filtering，rename
# --------------------------
chroms=$(seq -f "chr%g" 1 18) 

for CHR in "${chroms[@]}"; do
    echo "[$(date)] Processing ${CHR}"

    ls ${input}/*/*${CHR}.*g.vcf.gz > ${CHR}.list

    singularity exec -B ${input}:${input} /home/guilu/software/glnexus.sif \
        glnexus_cli --config gatk --list ${CHR}.list --dir ${CHR}.GLnexus.DB

    bcftools view ${CHR}.GLnexus.DB/cohort.bcf | bgzip -@ ${THREADS} -c > ${CHR}.cohort.vcf.gz
    tabix -p vcf ${CHR}.cohort.vcf.gz


    bcftools view \
        -m2 -M2 -v snps,indels \
        -i 'INFO/DP>=3 && FORMAT/GQ>=10 && INFO/MAF>0.05 && F_MISSING<0.2' \
        ${CHR}.cohort.vcf.gz \
        -Oz -o ${CHR}.filt.vcf.gz
    tabix -p vcf ${CHR}.filt.vcf.gz

    singularity exec /home/guilu/software/087.vcftools.sif \
        vcftools --gzvcf ${CHR}.filt.vcf.gz \
                 --max-missing 0.8 --maf 0.01 --min-alleles 2 --max-alleles 2 \
                 --recode --recode-INFO-all --stdout | bgzip -c > ${CHR}.filtered.vcf.gz
    tabix -p vcf ${CHR}.filtered.vcf.gz
    

    bcftools reheader -s rename.IDS2 ${CHR}.filtered.vcf.gz -o ${CHR}.rename.vcf.gz
    bcftools +fill-tags ${CHR}.rename.vcf.gz -- -t AC,AN | bgzip -c > ${CHR}.addAC.vcf.gz
    tabix -p vcf ${CHR}.addAC.vcf.gz
done
