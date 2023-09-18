
#### bedtools coverage, make sure entire UMR is inside of TE (-F 1), then only keep tes that have umr in them (awk)

## b73 names are different
wget https://de.cyverse.org/anon-files//iplant/home/shared/NAM/NAM_genome_and_annotation_Jan2021_release/DNA_METHYLATION_UMRs/NAM_UMRs_on_B73v5_coordinates/Zm-B73-REFERENCE-NAM-5.0_B73v5.ref_meth_UMR.bed
bedtools coverage -a <( zcat ../genomes_and_annotations/tes/new_edta_tes/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.split.gff3.gz | sed -e "s/B73_chr//g" ) -b Zm-B73-REFERENCE-NAM-5.0_B73v5.ref_meth_UMR.bed -F 1 | awk '$10 > 0' > B73-TE-UMR.bed



for nam in B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 Il14H Ki11 Ki3 Ky21 M162W M37W Mo18W Ms71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8 

do
if [ ! -f Zm-${nam}-REFERENCE-NAM-1.0.ref_meth_UMR.bed ]
then
wget https://de.cyverse.org/anon-files///iplant/home/shared/NAM/NAM_genome_and_annotation_Jan2021_release/DNA_METHYLATION_UMRs/NAM_UMRs/Zm-${nam}-REFERENCE-NAM-1.0.ref_meth_UMR.bed
fi
bedtools coverage -a <( zcat ../genomes_and_annotations/tes/new_edta_tes/${nam}.pseudomolecules-v*.fasta.mod.EDTA.TEanno.split.gff3.gz | sed -e "s/${nam}_chr//g" ) -b Zm-${nam}-REFERENCE-NAM-1.0.ref_meth_UMR.bed -F 1 | awk '$10 > 0' > ${nam}-TE-UMR.bed

done

## oh7b is annoying
wget https://de.cyverse.org/anon-files///iplant/home/shared/NAM/NAM_genome_and_annotation_Jan2021_release/DNA_METHYLATION_UMRs/NAM_UMRs/Zm-Oh7B-REFERENCE-NAM-1.0.ref_meth_UMR.bed
bedtools coverage -a <( zcat ../genomes_and_annotations/tes/new_edta_tes/Oh7b.pseudomolecules-v*.fasta.mod.EDTA.TEanno.split.gff3.gz | sed -e "s/Oh7b_chr//g" ) -b Zm-Oh7B-REFERENCE-NAM-1.0.ref_meth_UMR.bed -F 1 | awk '$10 > 0' > Oh7B-TE-UMR.bed
