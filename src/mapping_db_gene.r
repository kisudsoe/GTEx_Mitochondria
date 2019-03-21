suppressMessages(library(plyr))
mapping_go_gene = function(
        go,                                      # Input data; GO result table
        pval, Term_top=0, Term_size=0, category, # Criteria
        dir                                      # File out dir path
) {
    if(Term_size>0) go_sub = subset(go, Category%in%category & PValue<pval & Pop.Hits<Term_size)
    else go_sub = subset(go, Category%in%category & PValue<pval)
    cat(paste0('Input dim, rows= ',dim(go_sub)[1],' cols=',dim(go_sub)[2],'\n'))
    
    # 1. Filter GO results by criteria
    cat(paste0('(1/3) Filter GO results by criteria\n'))
    groups = unique(go$Group); n=length(groups)
    go_li = lapply(c(1:n), function(i) {
        cat(paste0('  - (',i,'/',n,') ',groups[i],' : '))
        terms_df = subset(go_sub,Group==groups[i])
        terms_df = terms_df[order(terms_df$PValue),] # ordered by pvalue
        if(Term_top>0) terms_ = terms_df[c(1:Term_top),]
        else terms_ = terms_df
        cat(paste0(paste(dim(terms_),collapse=', '),'\n'))
        return(terms_)
    })
    names(go_li) = groups
    #return(go_li)
    
    # 2. Generate GO-Gene table
    cat(paste0('\n(2/3) Generate GO-Gene table\n'))
    go_gene_li = lapply(c(1:n), function(i) {
        m = length(go_li[[i]]$Genes)
        go_gene_term = lapply(c(1:m),function(j) {
            genes = unique(unlist(strsplit(as.character(go_li[[i]]$Genes[j]),', ')))
            go_gene = data.frame(
                Genes     = genes,
                Term      = rep(go_li[[i]]$Term[j],length(genes)),
                Enrichment= rep(go_li[[i]]$Fold.Enrichment[j],length(genes)),
                FDR       = rep(go_li[[i]]$FDR[j],length(genes)),
                Term_size = rep(go_li[[i]]$Pop.Hits[j],length(genes)),
                Input_n   = rep(go_li[[i]]$Count[j],length(genes))
                                 )
            return(go_gene)
        })
        go_gene = ldply(go_gene_term,data.frame) # plyr
        go_gene = go_gene[order(go_gene$Genes,go_gene$Term_size),]
        cat(paste0('  - (',i,'/',n,') ',groups[i],' : Gene_n = ',
                   length(unique(go_gene$Genes)),'\n'))
        #print(dim(go_gene))
        return(go_gene)
    })
    names(go_gene_li) = groups
    #return(go_gene_li)
    
    # 3. Write tsv files
    cat(paste0('\n(3/3) Write tsv files\n'))
    lapply(c(1:n), function(i) {
        f_name = paste0(dir,'/go_gene_',category,'_',groups[i],'.tsv')
        write.table(go_gene_li[[i]],f_name,sep='\t',row.names=F,quote=F)
        cat(paste0('  - (',i,'/',n,') ',f_name,'\n'))
    })
    return(go_gene_li)
}

suppressMessages(library(dplyr))
summary_go_gene = function(dir,files) {
    # Input paths are results of 'mapping_go_gene' function.
    paths = paste0(dir,'/',files)
    n = length(paths)
    cat(paste0('Write file:\n'))
    k=lapply(c(1:n), function(i) {
        #cat(paste0('(',i,'/',n,') ',paths[i],'\n'))
        go_gene = read.delim(paths[i])
        df = go_gene %>% group_by(Genes) %>% slice(which.min(Term_size)) # Select only least term size
        df$Gene_n = apply(df,1,function(row) {
            sum(df$Term==row[2],na.rm=T)
        })
        df = df[order(-df$Gene_n),]
        
        names  = unlist(strsplit(files[i],'_'))
        f_name = paste0(names[c(4,7:length(names)-1)],collapse='_')
        f_path = paste0(dir,'_summ/',f_name,'.tsv')
        write.table(df,f_path,sep='\t',row.names=F)
        cat(paste0('  (',i,'/',n,') ',f_path,'\n'))
    })
    return()
}


mapping_reactome = function(
        reactome,                              # Input data; GO result table
        pval, Term_top=0, Term_size, category, # Criteria
        dir                                    # File out dir path
) {
    cat(paste0('Input dim, rows= ',dim(reactome)[1],' cols=',dim(reactome)[2],'\n'))
    re_sub = subset(reactome,Entities.pValue<pval)
    cat(paste0('Filtered dim, rows= ',dim(re_sub)[1],' cols=',dim(re_sub)[2],'\n'))
    
    m = length(re_sub$Submitted.entities.found)
    term_gene_li = lapply(c(1:m),function(i) {
        genes = unique(unlist(strsplit(as.character(re_sub$Submitted.entities.found[i]),';')))
        term_gene = data.frame(
            Genes     = genes,
            Term      = rep(re_sub$Pathway.name[i],length(genes)),
            Term_size = rep(re_sub$X.Entities.total[i],length(genes))
                            )
        return(term_gene)
    })
    term_gene_df = ldply(term_gene_li,data.frame) # plyr
    term_gene_df = term_gene_df[order(term_gene_df$Genes,
                                      term_gene_df$Term_size),]
    cat(paste0('  - Covered gene_n = ',length(unique(term_gene_df$Genes)),'\n'))
    #print(dim(term_gene_df))
    
    # 3. Write tsv files
    f_name = paste0(dir,'/Pathway_gene.tsv')
    write.table(term_gene_df,f_name,sep='\t',row.names=F,quote=F)
    return(term_gene_df)
}