suppressMessages(library(biomaRt))
biomart_gene_search = function(gene_input,db_version='hg19',only_chr=TRUE) {
    cat(paste0('Input length = ',length(gene_input)))
    if(db_version=='hg19') {
        cat(' >> biomart search: hg19\n')
        hg_gene = useMart(biomart='ENSEMBL_MART_ENSEMBL',host='grch37.ensembl.org',
                          dataset='hsapiens_gene_ensembl',path='/biomart/martservice')
    } else if(db_version=='hg38') {
        cat(' >> biomart search: hg38\n')
        hg_gene = useMart('ensembl',dataset='hsapiens_gene_ensembl')
    } else stop('db_version argument is either "hg19" or "hg38". Default is "hg19".')
    gene_attr = c(
    	'chromosome_name','start_position','end_position',
        'ensembl_gene_id','hgnc_symbol')
    genes = getBM(attributes = gene_attr,
                  filters    = 'hgnc_symbol',
                  values     = gene_input,
                  mart       = hg_gene)
    genes = as.data.frame(genes)
    cat(paste0('  >> Result: Ensembl table with rows= ',dim(genes)[1],', cols= ',dim(genes)[2],'\n'))
    
    # Filter only chromosomes
    if(only_chr==TRUE) {
        genes = subset(genes,chromosome_name%in%c(1:22,'X','Y'))
        cat(paste0('  >> Filter: Ensembl table with rows= ',dim(genes)[1],', cols= ',dim(genes)[2],'\n'))
    }
    return(genes)
}

biomart_ensg_search = function(ensg_ids,db_version='hg19',only_chr=TRUE) {
    cat(paste0('Input length = ',length(ensg_ids)))
    if(db_version=='hg19') {
        cat(' >> biomart search: hg19\n')
        hg_ensg = useMart(biomart='ENSEMBL_MART_ENSEMBL',host='grch37.ensembl.org',
                          dataset='hsapiens_gene_ensembl',path='/biomart/martservice')
    } else if(db_version=='hg38') {
        cat(' >> biomart search: hg38\n')
        hg_ensg = useMart('ensembl',dataset='hsapiens_gene_ensembl')
    } else stop('db_version argument is either "hg19" or "hg38". Default is "hg19".')
    gene_attr = c('chromosome_name','start_position','end_position',
                  'ensembl_gene_id','hgnc_symbol')
    genes = getBM(attributes = gene_attr,
                  filters    = 'ensembl_gene_id',
                  values     = ensg_ids,
                  mart       = hg_ensg)
    genes = as.data.frame(genes)
    cat(paste0('  >> Result: Ensembl table with rows= ',dim(genes)[1],', cols= ',dim(genes)[2],'\n'))
    
    # Filter only chromosomes
    if(only_chr==TRUE) {
        genes = subset(genes,chromosome_name%in%c(1:22,'X','Y'))
        cat(paste0('  >> Filter: Ensembl table with rows= ',dim(genes)[1],', cols= ',dim(genes)[2],'\n'))
    }
    return(genes)
}

biomart_rsid_search = function(rsids,db_version='hg19',only_chr=TRUE) {
    cat(paste0('Input length = ',length(rsids)))
    if(db_version=='hg19') {
        cat(' >> biomart search: hg19\n')
        hg_snp = useMart(biomart='ENSEMBL_MART_SNP',host='grch37.ensembl.org',
                         dataset='hsapiens_snp',path='/biomart/martservice')
    } else if(db_version=='hg38') {
        cat(' >> biomart search: hg38\n')
        hg_snp = useMart('ensembl',dataset='hsapiens_gene_ensembl')
    } else stop('db_version argument is either "hg19" or "hg38". Default is "hg19".')
    snp_attr = c('refsnp_id','ensembl_gene_stable_id',
                 'chr_name','chrom_start','chrom_end')
    snps = getBM(attributes = snp_attr,
                 filters    = 'snp_filter',
                 values     = rsids,
                 mart       = hg_snp)
    snps = as.data.frame(snps)
    cat(paste0('  >> Result: Ensembl table with rows= ',dim(snps)[1],', cols= ',dim(snps)[2],'\n'))
    
    # Filter only chromosomes
    if(only_chr==TRUE) {
        snps = subset(snps,chr_name%in%c(1:22,'X','Y'))
        cat(paste0('  >> Filter: Ensembl table with rows= ',dim(snps)[1],', cols= ',dim(snps)[2],'\n'))
    }
    return(snps)
}