
## on cbsu
export PYTHONPATH=/programs/liftoff-1.6.3.2/lib64/python3.9/site-packages:/programs/liftoff-1.6.3.2/lib/python3.9/site-packages
export PATH=/programs/liftoff-1.6.3.2/bin:$PATH

B73=/workdir/mcs368/NAM_genome_and_annotation_Jan2020_release/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fasta
PHZ51=/workdir/mcs368/merritt_anchorwave/5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta

liftoff -p 20 -g /workdir/mcs368/NAM_genome_and_annotation_Jan2020_release/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 -o 5_Buckler-PHZ51_mecatErrorCorrected.contigs.liftoff.gff3 $PHZ51 $B73
