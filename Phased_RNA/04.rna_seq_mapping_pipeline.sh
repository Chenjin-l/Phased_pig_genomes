
line0=${PBS_ARRAYID}

INPUT=`pwd`
OUTPUT=`pwd`
ref=~/02.GenomeIndex/03.Sus11.1/susScr11.hisat2

uid="MZ"
ls_date=`date +m%d%H%M%S`
temp=/tmpdisk/${uid}_${ls_date}_${PBS_ARRAYID}
mkdir ${temp}
cd ${temp}

ID=`sed -n "${line0}p" ${INPUT}/Muscle-RNAseq.fqfiles|awk '{print $1}'`
fq1=`sed -n "${line0}p" ${INPUT}/Muscle-RNAseq.fqfiles|awk '{print $2}'`
fq2=`sed -n "${line0}p" ${INPUT}/Muscle-RNAseq.fqfiles|awk '{print $3}'`
##1.filter
fastp --thread=40  -l 75 -c -i $fq1  -o ${ID}_RNAClean_1.fq.gz -I  $fq2  -O ${ID}_RNAClean_2.fq.gz -j ${ID}.json -h ${ID}.html

##2.Mapping
/home/jxlabgdp/01.Biosoft/hisat2-2.2.0/hisat2 -x $ref -p 100 -I 0 --qc-filter -X 500 --dta -1 --rna-strandness RF \
${ID}_RNAClean_1.fq.gz -2  ${ID}_RNAClean_2.fq.gz -S ${ID}.sam
samtools sort -@30 ${ID}.sam -o ${ID}.sorted.bam
samtools index -@30 ${ID}.sorted.bam
sambamba markdup  -t 10  ${ID}.sorted.bam ${ID}.sorted.markdup.bam --overflow-list-size 600000
samtools index ${ID}.sorted.markdup.bam -@30


cp ${ID}.sorted.bam* ${ID}.sorted.markdup.bam*  $OUTPUT
rm -r ${temp}
