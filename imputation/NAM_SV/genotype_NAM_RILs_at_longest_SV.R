


## read in haplotypes for each ril
all.haps=read.table('../phg_hapids_namrils_dbout.txt.gz', header=T)
all.haps[all.haps==-1]=NA

## for B73, find refranges that overlap longest sv
## then make a table of which parent they have in this region
## ideally, the majority of refranges will have one parent or the other (not expecting recombination in the sv)

