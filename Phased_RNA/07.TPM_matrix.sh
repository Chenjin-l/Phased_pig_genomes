#TPM matrix
name1=$1
ls ./ |grep  "tsv$" >all.sample.geneCounts.files
for i in `cat all.sample.geneCounts.files `
do
        name=`echo $i |awk -F . '{print $1}'`
        csvtk uniq -f 1 -Tt $i >$i.rmdup
        cat $i.rmdup | grep -Ev "chrX|chrY|chrM|NW" |awk '{print $1,$NF}' |sed 's/gene://g' |sed "s/TPM/$name/g" |tr " " "\t" |csvtk uniq -f 1 -Tt  >${i%%.*}.TMP.tmp
done
csvtk join -Tt `ls *TMP.tmp` >$name1.TPM.txt

