
## weird way to run emboss stuff on cbsu
singularity run -B $PWD --pwd $PWD /programs/EMBOSS/emboss.sif infoseq -auto -only -name -length -noheading NAM.EDTA2.0.0.MTEC02052020.TElib.fa > NAM.EDTA2.0.0.MTEC02052020.TElib.contigsizes.txt
