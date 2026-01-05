out=path/to/output
pvcf=path/to/${CHR}.shapiet.vcf.gz  # phased genotype VCF files
pmap=path/to/map                    # genetic map files

export _JAVA_OPTIONS="-Djava.io.tmpdir=/home/Mzhou/02.F2/log/tmp"

# 遍历每条染色体
for chr in {1..18}; do

    # 运行 Hap-IBD
    java -Xmx50g -jar /home/Mzhou/02.F2/software/hap-ibd.jar \
        gt=${pvcf}/${CHR}.shapiet.vcf.gz \
        map=${pmap}/chr${chr}.chr.map \
        min-output=0.01 \
        min-seed=0.01 \
        max-gap=1000 \
        min-markers=200 \
        min-extend=0.01 \
        out=${chr}.out

zcat ${chr}.out.ibd.gz | awk '$2==1 && $4==1 {
        pair = ($1 < $3) ? $1 "_" $3 : $3 "_" $1;
        sum[pair] += $NF;
    } END { for(p in sum) print p, sum[p]; }' > ${chr}.paternal_ibd.txt
zcat ${chr}.out.ibd.gz | awk '$2==2 && $4==2 {
        pair = ($1 < $3) ? $1 "_" $3 : $3 "_" $1;
        sum[pair] += $NF;
    } END { for(p in sum) print p, sum[p]; }' > ${chr}.maternal_ibd.txt

join <(sort ${chr}.paternal_ibd.txt) <(sort ${chr}.maternal_ibd.txt) | \
        awk '{mean=($2+$3)/2; printf "%s %.6f %.6f %.6f %.6f\n",$1,$2,$3,$2+$3,mean}' \
        > ${chr}.normalized_ibd.txt    
done
