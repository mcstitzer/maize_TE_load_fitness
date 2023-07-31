# maize_TE_load_fitness
Do TEs decrease maize fitness?



1. Download genomes and annotations, in `genomes_and_annotations/`
2. Download phenotypes, in `phenotypes/`
3. Count TEs and impute parental genotypes onto RILs using the PHG `imputation/`
  - Output haplotype size and TE content for each reference range, and summed for each RIL genome-wide
  - Measure and categorize TEs in haplotypes based on genes and methylation
4. Correlate individual TE and haplotype counts to phenotypes, in `models/`
5. Figures in `figures/`
