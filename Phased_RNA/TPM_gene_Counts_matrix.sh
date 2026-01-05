#TPM matrix
name1=F2
ls ./ |grep  "tsv$" >all.sample.geneCounts.files
for i in `cat all.sample.geneCounts.files `
do
        name=`echo $i |awk -F . '{print $1}'`
        csvtk uniq -f 1 -Tt $i >$i.rmdup
        cat $i.rmdup | grep -Ev "chrX|chrY|chrM|NW" |awk '{print $1,$NF}' |sed 's/gene://g' |sed "s/TPM/$name/g" |tr " " "\t" |csvtk uniq -f 1 -Tt  >${i%%.*}.
done
csvtk join -Tt `ls *TMP.tmp` >$name1.TPM.txt

#gene Counts matrix
ls ./ |grep  "featureCounts.txt$" >all.sample.featureCounts.files
for i in `cat all.sample.featureCounts.files`
do
         dirname `fgrep "Geneid" $i |awk '{print $NF}'`  |sed 's/\//\\\//g' |awk '{print $0"\\\/"}' >$i.dir.tmp
         dir=`cat $i.dir.tmp`
         bufix_raw=`fgrep "Geneid" $i |awk -F "/" '{print $NF}'`
         bufix=`fgrep "Geneid" $i |awk -F "/" '{print $NF}' |awk -F . '{print $1}'`     
         cat $i | grep -Ev "chrX|chrY|chrM|NW" | awk '{print $1,$NF}' | sed "s/$dir//g" |sed "s/$bufix_raw/$bufix/g" |sed 's/gene://g' |fgrep -v "#" |tr " " "

done
csvtk join -Tt `ls *FC.tmp` >$name1.GeneCounts.txt
#rm *tmp
