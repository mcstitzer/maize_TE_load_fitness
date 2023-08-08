

# wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cds.fa.gz

function map_fq(){  minimap2 -ax sr -t 2 Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cds.fa ${1}_1.clean.fq.gz ${1}_2.clean.fq.gz | samtools view -Sb - | samtools sort - -O bam -o ${1}.B73CDS.bam ; }

export -f map_fq
cat 282_fqnames.txt | parallel --progress --jobs 31 map_fq {}


for i in *B73CDS.bam; do samtools index $i & done

for i in *B73CDS.bam; do samtools flagstat $i > ${i%.*}.B73geneflagstat.txt & done
