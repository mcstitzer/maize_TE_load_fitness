
#singularity pull EDTA.sif docker://oushujun/edta:2.0.0
#
#export LC_ALL=C
#singularity exec --no-home ./EDTA.sif EDTA.pl --genome 5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta --cds genome.cds.fa --curatedlib NAM.EDTA2.0.0.MTEC02052020.TElib.fa --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 10


## or dang, why not just repeatmask with the NAM library - that's got to be closer...
export PATH=/programs/RepeatMasker_4-1-0:$PATH
RepeatMasker -e ncbi -pa 24 -q -no_is -norna -nolow -div 40 -lib NAM.EDTA2.0.0.MTEC02052020.TElib.fa -cutoff 225 5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta
/programs/RepeatMasker_4-1-0/util/rmOutToGFF3.pl 5_Buckler-PHZ51_mecatErrorCorrected.contigs.fasta.out > 5_Buckler-PHZ51_mecatErrorCorrected.contigs.RepeatMasker.gff3
