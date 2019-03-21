#!/usr/bin/env Rscript
# This file is for downloading gene location data from Ensembl biomaRt

## Command Arg Parameters ##
# Rscript src/biomart_gene.r
args = commandArgs(trailingOnly=T)
if(length(args)>0) stop('Rscript biomart_gene.r, no argument is needed.')

# system parameters
suppressMessages(library(biomaRt))
source('src/pdtime.r') # pdtime(time,1/2); 1= Job done, 2= Job process
t0 = Sys.time(); cat(paste0('Process initiate at ',t0,'\n\n'))

#########################
## Function start here ##
#########################
# Ensembl biomart (grch37)
hg19_gene = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
					dataset="hsapiens_gene_ensembl", path="/biomart/martservice")
#print(subset(listAttributes(hg19_gene),page=='feature_page')[1:2,])
#print(listFilters(hg19_gene))
gene_attr = c('chromosome_name','start_position','end_position','ensembl_gene_id','external_gene_name')
genes = getBM(attributes = gene_attr,
			  filters    = '',
			  values     = '',
			  mart       = hg19_gene)
genes = as.data.frame(genes)
cat(paste0('\tEnsembl table, rows= ',dim(genes)[1],' cols= ',dim(genes)[2],'\n'))
f.name1 = paste0('db/ensembl_gene_ann.tsv')
write.table(genes,f.name1,row.names=F,col.names=T,quote=F,sep='\t')
cat(paste0('\tFile write: ',f.name1,'\n'))

genes_filt = subset(genes,chromosome_name%in%c(1:22,'X','Y'))
chr_ = paste0('chr',genes_filt[,1])
genes_ = data.frame(chr_,genes_filt[,2:5])
colnames(genes_) = c('chr','start','end','ENSGid','Symbol')
cat(paste0('\n\tFilter result, rows= ',dim(genes_)[1],' cols= ',dim(genes_)[2],'\n'))
f.name2 = paste0('db/ensembl_gene.bed')
write.table(genes_[,1:4],f.name2,row.names=F,col.names=F,quote=F,sep='\t')
cat(paste0('\tFile write: ',f.name2,'\n'))
cat(paste0('\n',pdtime(t0,1),'\n'))
##################
## Function end ##
##################