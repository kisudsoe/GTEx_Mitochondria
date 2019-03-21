#!/usr/bin/env Rscript
# This file is for filtering Regulome data

## Command Arg Parameters ##
# CMD command1: Rscript regulome.r data/seedSNP_1817.bed db/RegulomeDB.dbSNP132.Category1.txt.rds db/RegulomeDB.dbSNP132.Category2.txt.rds
# CMD command2: Rscript regulome.r random/seeds/seedSNP_1817.bed db/RegulomeDB.dbSNP132.Category1.txt.rds db/RegulomeDB.dbSNP132.Category2.txt.rds
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript regulome.r [SNP_BED_file_path] [regulome_file_path_1] [regulome_file_path_2] [regulome_file_path_3]
  - Argument [SNP_file_path] is mendatory.
  - At least one and up to three [regulome_file_path] arguments are needed to be supplied.
  - The SNPs with high Regulome score is filtering by as >=2b.'
if(length(args)<2|length(args)>4) stop(hmsg)
f.path = NULL
for(i in 1:length(args)) {
	if(i==1) snp.path = unlist(args[i])
	else if(i>=2) f.path = c(f.path,unlist(args[i]))
}

# System parameter
library(tools)
library(plyr)
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
snp = read.delim(snp.path,header=F)
colnames(snp) = c('chr','start','end','ann')
snpids = gsub('(.*)_.*','\\1',snp[,4])
cat(paste0('Input SNPs number = ',length(snpids),'\n'))

n = length(f.path); cat(paste0('Input regulome files = ',n,'\n'))
reg.li = list()
for(i in 1:n) {
	cat(paste0('\nRead input file ',i,': ',f.path[i],'\n'))
	if(file_ext(f.path[i])=='gz') reg.li[[i]] = try(read.delim(gzfile(f.path[i]),header=F))
    else if(file_ext(f.path[i])=='rds') reg.li[[i]] = try(readRDS(f.path[i]))
	else reg.li[[i]] = try(read.delim(f.path[i],header=F))
	colnames(reg.li[[i]]) = c('chr','id','rsid','description','level')
	cat(paste0('Table row = ',dim(reg.li[[i]])[1],', col = ',dim(reg.li[[i]])[2],'\n'))
	cat(paste0(pdtime(t0,2),'\n'))
}
reg.df = ldply(reg.li) # plyr
reg_1f_only = subset(reg.df,level%in%c("1f"))
reg_2b = subset(reg.df,level%in%c("1a","1b","1c","1d","1e","1f","2a","2b"))
cat(paste0('\nRegulome score >=2b, SNPs = ',nrow(reg_2b),'\n'))
cat(paste0('Functional motifs (1f_only-2b) = ',nrow(reg_2b)-nrow(reg_1f_only),'\n'))

snpids.1f_only = subset(reg_1f_only,rsid%in%snpids)
snpids.2b      = subset(reg_2b,     rsid%in%snpids)
cat(paste0('  - Regulome >=2b SNPs = ',nrow(snpids.2b),'\n'))
cat(paste0('  - SNPs with functional motifs (1f_only-2b) = ',
                nrow(snpids.2b)-nrow(snpids.1f_only),'\n'))
f.name1 = paste0('data/regulome_',nrow(snpids.2b),'.tsv')
write.table(snpids.2b,f.name1,row.names=F,col.names=T,quote=T,sep='\t')
cat(paste0('\n>> File write: ',f.name1,'\n'))

snp.rsid= data.frame(snp,rsid=snpids)
snp.bed = subset(snp.rsid,rsid%in%snpids.2b$rsid)
f.name2 = paste0('data/snp_',nrow(snp.bed),'_regulome2b.bed')
write.table(snp.bed[,1:4],f.name2,row.names=F,col.names=F,quote=T,sep='\t')
cat(paste0('>> File write: ',f.name2,'\n'))
cat(paste0(pdtime(t0,1),'\n'))
##################
## Function end ##
##################