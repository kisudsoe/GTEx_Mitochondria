bedtools sort -i data/gwas_rsid_4864.bed > data/gwas_rsid_4864_sort.bed
bedtools closest -d -a data/gwas_rsid_4864_sort.bed -b db/roadmap_enh_merge.bed > data/roadmap_dist.tsv
#Rscript src/bedtools_closest_roadmap.r data/roadmap_dist.tsv

bedtools closest -d -a data/gwas_rsid_4864_sort.bed -b db/encode_tfbs_merge.bed > data/encode_dist.tsv
#Rscript src/bedtools_closest_roadmap.r data/encode_dist.tsv

#Rscript src/regulome.r data/gwas_rsid_4864_sort.bed db/RegulomeDB.dbSNP132.Category1.txt.rds db/RegulomeDB.dbSNP132.Category2.txt.rds

#Rscript src/lncrnasnp.r data/gwas_rsid_4864_sort.bed db/lncRNASNP2_snplist.txt.rds db/lncrnas.txt.rds db/lncrna-diseases_experiment.txt.rds