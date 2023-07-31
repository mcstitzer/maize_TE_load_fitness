
### Filter out reference ranges we can't trust

1. to filter reference ranges by b73 genotyping quality, and by range of haplotype length
  - Use genotyping of 38 b73 samples to remove refranges that are hard to genotype, or might be a different flavor of B73 `filter_refranges_B73variation.R`
  - Use assembly haplotypes of NAM to determine when there's a greater than 1Mb difference between largest and smallest allele of a haplotype `filter_refranges_minMaxSize.R`
    - also combines with b73genotyping results to generate final filtering file
  - Filtering file of which refranges to KEEP: `refranges_B73correctlygenotypedAND1Mbrangeremoved.2022-03-22.txt`
  
  
### Generate matrix of haplotype length, TE length, each repeat class length

2. `count_tes_in_haplotypes.R` outputs a table of each haplotype with columns for reference range ID, haplotype length, TE bp, and bp of each repeat class
 - uses the header from hvcf files to get regions for each haplotype
 - `zcat nam_unmerged_haplotypes.vcf.gz | grep "##" > nam_unmerged_haplotypes.haplines`
 - outputs `allNAM_hapids.TEbp.sup.2022-05-31.txt`, also from an updated TE annotation `allNAM_hapids.TEbpUpdate.sup.2022-07-12.txt`. I'll move forward with this updated annotation
3. `hapids_to_taxa.R` takes those hapids and tells which are present in a given RIL
 - This is all done through the PHG (maize_1_0), these are haploid mappings of NAM RIL GBS data, run against all NAM parents
 - Info on connecting to BRAPI server here [https://bitbucket.org/bucklerlab/rphg/wiki/Home#markdown-header-brapi-access]


### Generate PCs for population structure

By Merritt, documented here [https://github.com/merkhaiburch/misc_scripts/blob/master/tassel/MDS_pcs.sh]

Briefly, in TASSEL, subsampled genotypes to 400k random sites per chromosome, merged all chromosomes together, calculated a distance matrix with all the defaults, and calculated MDS's

