
bfile=path/to/bfile   #

king -b $bfile --ibdseg 

awk '{print $1"_"$3,$7}' king.seg > individual.king.out

awk '{print $1,$5}' kingallsegs.txt > chr_length.txt
