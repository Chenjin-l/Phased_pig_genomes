# ---------------------------
# Environment setup
# ---------------------------
export _JAVA_OPTIONS="-Djava.io.tmpdir=/home/Mzhou/02.F2/2025-12-23/log/tmp"
conda activate rnaqc   # Activate conda environment if needed
export RTG_JAR=/home/Mzhou/anaconda3/envs/rnaqc/share/rtg-tools-3.12.1-0/RTG.jar
export PATH=/home/Mzhou/anaconda3/envs/rnaqc/bin:$PATH

# ---------------------------
# File paths
# ---------------------------
VCF="/home/Mzhou/02.F2/2025-7-7/vcf/merge.vcf.gz"
SDF="/home/Mzhou/00.Genome/susScr11_analysis_set.sdf"
TRIO_INFO="/home/Mzhou/02.F2/2025-7-7/vcf/mender/1empt_new.txt"
OUTPUT_DIR="/home/Mzhou/02.F2/2025-7-7/vcf/mender/output"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# Get current line based on PBS array ID
line0=${PBS_ARRAYID}
line=$(sed -n "${line0}p" $TRIO_INFO)

# ---------------------------
# Iterate over each F2 sample, generate PED file, and run RTG Mendelian
# ---------------------------

SON=$(echo $line | awk '{print $1}')
FATHER=$(echo $line | awk '{print $2}')
MOTHER=$(echo $line | awk '{print $3}')
SEX=$(echo $line | awk '{print $4}')

PED_FILE="${SON}.ped"

# Generate PED file
cat <<EOM > $PED_FILE
1 $SON $FATHER $MOTHER $SEX 0
1 $FATHER 0 0 1 0
1 $MOTHER 0 0 2 0
EOM

# Extract trio VCF
bcftools view -s "$SON,$FATHER,$MOTHER" -Oz -o ${SON}.trio.vcf.gz $VCF
bcftools index ${SON}.trio.vcf.gz

# Run RTG Mendelian analysis
rtg mendelian \
    -i ${SON}.trio.vcf.gz \
    -o ${SON}.trio.annotated.vcf.gz \
    --pedigree $PED_FILE \
    -t $SDF \
    | tee ${SON}.rtg.log

# Extract low concordance report (<98%)
awk '
  /Family:/ { family = $0 }
  /Concordance/ {
    match($0, /\( *([0-9.]+)\%\)/, a)
    percent = a[1]
    if (percent < 98) {
        print family
        print $0
    }
  }
' ${SON}.rtg.log > ${SON}.low_concordance.txt

# Clean up intermediate files
rm -f $PED_FILE ${SON}.trio.vcf.gz ${SON}.trio.vcf.gz.csi


