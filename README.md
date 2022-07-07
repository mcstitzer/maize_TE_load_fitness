# maize_TE_load_fitness
Do TEs decrease maize fitness?

1. Download annotations from maizegdb, `genomes_and_annotations/`
2. Impute parental genotypes onto RILs using the PHG `imputation/`
  - Output as haplotype size, TE content for each reference range, and summed for the individual RIL genome-wide
4. Correlate individual TE and haplotype counts to phenotypes `phenotypes/`
5. Correlate TE and haplotype counts to segregation distortion in NAM families `segregation_distortion/`
6. Correlate TE and haplotype counts to other features of the region `phg_distance_matrices/`
