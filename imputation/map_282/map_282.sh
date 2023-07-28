
## do this all on a cbsu machine cbsulm05, with 64 cpu

## clean up anmes to get a thign to call files
sed -e 's/^.*clean_data\///g' ~/transfer/282_cbsufsrv4_paths.txt | sed -e 's/_1.clean.fq.gz//g' | sed -e 's/_2.clean.fq.gz//g' | uniq > 282_fqnames.txt

scp cbsufeschotte2:/local/workdir/mcs368/maize_TE_load_fitness/genomes_and_annotations/tes/NAM.EDTA2.0.0.MTEC02052020.TElib.fa .

CPU=64

while read name
do

R1source=/data2/maizediv/illumina/RawSeqData/WGS/Zea/*/data_release/clean_data/${name}_1.clean.fq.gz
R2source=/data2/maizediv/illumina/RawSeqData/WGS/Zea/*/data_release/clean_data/${name}_2.clean.fq.gz

if [ ! -f ${name}.NAMTElib.bam ]
then

scp cbsufsrv4:$R1source .
scp cbsufsrv4:$R2source .
echo "reads transferred from $name"

R1=${name}_1.clean.fq.gz
R2=${name}_2.clean.fq.gz

minimap2 -ax sr -t $CPU NAM.EDTA2.0.0.MTEC02052020.TElib.fa $R1 $R2 | samtools view -Sb - | samtools sort - -O bam -o ${name}.NAMTElib.bam

scp ${name}.NAMTElib.bam cbsufeschotte2:/workdir/mcs368/maize_TE_load_fitness/282_load/

echo "mapped sequence back to cbsufeschotte2"

fi

done < 282_fqnames.txt
