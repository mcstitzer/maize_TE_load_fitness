taxon=Buckler-PHZ51
readsIn=mecatErrorCorrectedReads_try4/$taxon/1-consensus/cns_final.fasta

canu -d canuAssembly/$taxon -p ${taxon}_mecatErrorCorrected -trim-assemble -saveReads=true -genomeSize=2.5g -ovlMerThreshold=500 -useGrid=false -gnuplot=/usr/bin/gnuplot -corrected -pacbio $readsIn
