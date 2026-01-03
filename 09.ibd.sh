out=path/to/output
pvcf="/home/Mzhou/02.F2/z02data/02vcf/02genevcf"   # phased genotype VCF files
pmap=path/to/map       # genetic map files

export _JAVA_OPTIONS="-Djava.io.tmpdir=/home/Mzhou/02.F2/log/tmp"

# 遍历每条染色体
for chr in {1..18}; do

    # 运行 Hap-IBD
    java -Xmx50g -jar /home/Mzhou/02.F2/software/hap-ibd.jar \
        gt=${pvcf}/chr${chr}.vcf.gz \
        map=${pmap}/1chr${chr}.chr.map \
        min-output=0.01 \
        min-seed=0.01 \
        max-gap=1000 \
        min-markers=200 \
        min-extend=0.01 \
        out=${chr}.out
