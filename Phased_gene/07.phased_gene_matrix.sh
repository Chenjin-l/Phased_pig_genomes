
Inpath=path/to/F2-*featureCounts.txt

sed '1d' ${Inpath}/F2-503.Hap1.featureCounts.txt|cut -f1,6 > F2ind.quantile
pastestr='paste -d " " F2ind.quantile'
echo "Genename ratio"> biased.genes
for i in `ls ${Inpath}/F2-*.featureCounts.txt`;do
    id=`basename $i .featureCounts.txt`
    cut -f7 ${Inpath}/${id}.featureCounts.txt|sed '1,2d' >id.countall
    cut -f7 ${Inpath}/${id}.Hap1.featureCounts.txt|sed '1,2d'>id.countHap1
    cut -f7 ${Inpath}/${id}.Hap2.featureCounts.txt|sed '1,2d'>id.countHap2
    paste id.countall id.countHap1 id.countHap2|awk '{if(($2+$3)>5) {print $1,$2/($2+$3)*$1,$3/($2+$3)*$1} \
    else {print $1,$1/2,$1/2}}'|sed "1i${id} ${id}Hap1 ${id}Hap2" >${id}.count
    #Summary for Parents biased expression 
    cut -f1 F2ind.quantile |paste id.countall id.countHap1 id.countHap2 - |awk '{if(($2+$3)>5) print $4,($2+0.1)/($3+0.1)}' >> biased.genes
    rm id.countall id.countHap1 id.countHap2
    pastestr=$pastestr' '${id}'.count'
done
