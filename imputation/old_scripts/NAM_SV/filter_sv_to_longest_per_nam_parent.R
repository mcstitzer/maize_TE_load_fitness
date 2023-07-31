library(rtracklayer)
library(stringr)
library(data.table)

files=list.files(pattern='*.bed')
genomes=str_split_fixed(files, '-', 3)[,2]

allbionano=lapply(genomes, function(x){
  a=fread(files[genomes==x])
  a$length=as.numeric(gsub('length=', '', a$V4) )
  a$genome=x
  a=a[a$V6=='method=Bionano-PAF' & a$length>1e6,]
  return(a)
})

sv=do.call(rbind, allbionano)
## get longest sv for each genome
## for all buut cml277, it's an inversion, for cml277 it's a deletion
longest=sv %>% group_by(genome) %>% filter(length==max(length)) %>% data.frame()


write.table(longest, 'longest_sv_per_nam.txt', quote=F, sep='\t', row.names=F, col.names=T)




