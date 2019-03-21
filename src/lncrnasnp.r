#!/usr/bin/env Rscript
# This file is for filtering RoadMap data

## Command Arg Parameters ##
# CMD command: Rscript lncrnasnp.r data/seedSNP_1817.bed db/lncRNASNP2_snplist.txt.rds db/lncrnas.txt.rds db/lncrna-diseases_experiment.txt.rds
args = commandArgs(trailingOnly=T)
hmsg = 'Rscript lncrnasnp.r [SNP_BED_file_path] [lncRNAsnp2_SNP_list_file_path] [lncRNAsnp2_lncRNA_list_file_path] [lncRNAsnp2_diseases_list_file_path]
  - [SNP_BED_file_path] is a mendatory argument for your SNP list with BED format.
  - [lncRNAsnp2_SNP_list_file_path] is a mendatory argument for SNP-lncRNA pair list.
  - [lncRNAsnp2_lncRNA_list_file_path] is a mendatory argument for lncRNA annotation.
  - [lncRNAsnp2_diseases_list_file_path] is an optional argument for experimental validate lncRNA-associated diseases.'
if(length(args)<3|length(args)>4) stop(hmsg)
snp_path    = args[1]
snplnc_path = args[2]
ann_path    = args[3]
if(length(args)>3) dis_path = args[4]

# System parameter
source('src/pdtime.r')
t0 = Sys.time()

#########################
## Function start here ##
#########################
# 1. Read files..
cat("\n(1/3) Read files..\n")
snp = read.delim(snp_path,header=F)
colnames(snp) = c('chr','start','end','snp_id')
dbsnp = gsub('(.*)_.*','\\1',snp[,4])
snp.df = cbind(snp,dbsnp)
cat(paste0('  - ',snp_path,'; ',pdtime(t0,2),'\n'))

snplnc = readRDS(snplnc_path)
cat(paste0('  - ',snplnc_path,'; ',pdtime(t0,2),'\n'))

ann = readRDS(ann_path)
colnames(ann)[1] = 'lncRNA'
cat(paste0('  - ',ann_path,'; ',pdtime(t0,2),'\n'))
df = data.frame(path=c(snp_path,   snplnc_path,   ann_path),
				nrow=c(dim(snp)[1],dim(snplnc)[1],dim(ann)[1]),
				ncol=c(dim(snp)[2],dim(snplnc)[2],dim(ann)[2]))

if(length(args)>3) { # Disease option
	dis = readRDS(dis_path)
	cat(paste0('  - ',dis_path,'; ',pdtime(t0,2),'\n'))
	df_dis= data.frame(path=dis_path,nrow=dim(dis)[1],ncol=dim(dis)[2])
	df = rbind(df,df_dis)
}
knitr::kable(df)
cat(paste0('\n  - ',pdtime(t0,2),'\n'))

# 2. Overlapping lncRNA to my SNP list and binding annotation..
cat("\n(2/3) Overlapping lncRNA to my SNP list and binding annotation..\n")
snp_lnc = merge(snp.df,snplnc,by='dbsnp')
df2 = data.frame(lncRNA=c(length(unique(snp_lnc$lncRNA))),
				 SNPs  =c(length(unique(snp_lnc$snp_id))))
knitr::kable(df2)
f_name1 = paste0('data/snp_',length(unique(snp_lnc$snp_id)),'_lncrnasnp.bed')
write.table(unique(snp_lnc[,2:5]),f_name1,row.names=F,col.names=F,quote=T,sep='\t')
cat(paste0('\n>> File write: ',f_name1,'\n'))

# 3. Annotating SNPs in lncRNAs
cat("\n(3/3) Annotating SNPs in lncRNAs\n")
snp_lnc_ann = merge(snp_lnc[,c(1,6)],ann,by='lncRNA')
if(length(args)>3) { # Disease option
	snp_lnc_ann_dis = merge(snp_lnc_ann,dis,by='lncRNA',all.x=T)
	f_name2 = paste0('data/lncrnasnp_',length(unique(snp_lnc_ann_dis$dbsnp)),'.tsv')
	write.table(snp_lnc_ann_dis,f_name2,row.names=F,col.names=T,quote=F,sep='\t')
	cat(paste0('\n>> File write: ',f_name2,'\n'))
} else {
	cat("\n>> Disease associated SNP option does not processed.\n")
    f_name2 = paste0('data/lncrnasnp_',length(unique(snp_lnc_ann$dbsnp)),'.tsv')
	write.table(snp_lnc_ann,f_name2,row.names=F,col.names=T,quote=F,sep='\t')
}
cat(paste0('\n',pdtime(t0,1),'\n'))
##################
## Function end ##
##################