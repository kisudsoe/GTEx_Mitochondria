suppressMessages(library(piano))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
gsea = function(
    organs,organs_df,gwas_df_ann,
    myGSC,  # trait, efo, cat (only for Adipose_Subcutaneous)
    myStat, # fisher, wilcoxon, mean, median, page
    min=1   # limit min gene number
) {
    for(i in 1:length(organs)) {
        cat(organs[i])
        # Filtering data by each organ
        organ_which1 = which(organs_df[,i]=='Y')
        organ_which2 = which(gwas_df_ann$Organs==organs[i])
        organ_which  = intersect(organ_which1,organ_which2)
        gwas_df_ann_sub = gwas_df_ann[organ_which,] %>% unique
        n = dim(gwas_df_ann_sub)
        cat(paste0('; rows = ',n[1],', cols = ',n[2],'\n'))
        
        # Preparing pre-dataset
        cols = c(
            'gene_id','Age.Coef','p.value','DISEASE.TRAIT',
            'Primary_EFO','Category2','old','young')
        col_which = which(colnames(gwas_df_ann_sub)%in%cols)
        gwas_df_df = gwas_df_ann_sub[,col_which] %>% unique
        colnames(gwas_df_df) = c(
            'gene','coef','pval','trait','efo','cat','old','young')
        gwas_df_df$old_t = apply(gwas_df_df[,c(7,4)],1,paste0,collapse='_')
        gwas_df_df$young_t = apply(gwas_df_df[,c(8,4)],1,paste0,collapse='_')
        
        # Preparing input dataset
        if(myStat%in%c('wilcoxon')) {
            gwas_df_coef = data.frame(
                gene = gwas_df_df$gene,
                coef = gwas_df_df$coef
            ) %>% unique
            myGeneData = data.frame(coef=gwas_df_coef$coef)
            rownames(myGeneData) = gwas_df_coef$gene
            myMethod = 'nullDist'# 'geneSampling'
            myadj = 'fdr'
        } else {
            gwas_df_pval = data.frame(
                gene = gwas_df_df$gene,
                pval = gwas_df_df$pval
            ) %>% unique
            myGeneData = data.frame(pval=gwas_df_pval$pval)
            rownames(myGeneData) = gwas_df_pval$gene
            myMethod = 'nullDist'
            myadj = 'bonferroni'
        }
        
        # Preparing geneset groups
        if(myGSC == 'trait') {
            gwas_df_trait_ = data.frame(
                gene  = gwas_df_df$gene,
                trait = gwas_df_df$trait
            ) %>% unique
            myGsc = loadGSC(gwas_df_trait_)
            
            gwas_df_coef = data.frame(
                gene = gwas_df_df$gene,
                coef = gwas_df_df$coef,
                trait = gwas_df_df$trait
            ) %>% unique
        } else if(myGSC == 'efo') {
            gwas_df_efo_ = data.frame(
                gene  = gwas_df_df$gene,
                trait = gwas_df_df$efo
            ) %>% unique
            myGsc = loadGSC(gwas_df_efo_)
            
            gwas_df_coef = data.frame(
                gene = gwas_df_df$gene,
                coef = gwas_df_df$coef,
                trait = gwas_df_df$efo
            ) %>% unique
        } else if(myGSC == 'cat') {
            gwas_df_cat_ = data.frame(
                gene  = gwas_df_df$gene,
                cat   = gwas_df_df$cat
            ) %>% unique
            myGsc = loadGSC(gwas_df_cat_)
            
            gwas_df_coef = data.frame(
                gene = gwas_df_df$gene,
                coef = gwas_df_df$coef,
                trait = gwas_df_df$cat
            ) %>% unique
        } else if(myGSC == 'net_old') {
            gwas_df_old = data.frame(
                gene  = gwas_df_df$gene,
                old   = gwas_df_df$old
            ) %>% unique
            myGsc = loadGSC(gwas_df_old)
            
            gwas_df_coef = data.frame(
                gene = gwas_df_df$gene,
                coef = gwas_df_df$coef,
                trait = gwas_df_df$old
            ) %>% unique
        } else if(myGSC == 'net_old_trait') {
            gwas_df_young = data.frame(
                gene  = gwas_df_df$gene,
                old   = gwas_df_df$old_t
            ) %>% unique
            myGsc = loadGSC(gwas_df_young)
            
            gwas_df_coef = data.frame(
                gene = gwas_df_df$gene,
                coef = gwas_df_df$coef,
                trait = gwas_df_df$old_t
            ) %>% unique
        } else if(myGSC == 'net_young') {
            gwas_df_young = data.frame(
                gene  = gwas_df_df$gene,
                young = gwas_df_df$young
            ) %>% unique
            myGsc = loadGSC(gwas_df_young)
            
            gwas_df_coef = data.frame(
                gene = gwas_df_df$gene,
                coef = gwas_df_df$coef,
                trait = gwas_df_df$young
            ) %>% unique
        } else if(myGSC == 'net_young_trait') {
            gwas_df_young = data.frame(
                gene  = gwas_df_df$gene,
                young = gwas_df_df$young_t
            ) %>% unique
            myGsc = loadGSC(gwas_df_young)
            
            gwas_df_coef = data.frame(
                gene = gwas_df_df$gene,
                coef = gwas_df_df$coef,
                trait = gwas_df_df$young_t
            ) %>% unique
        }
            
        # run GSA
        gsaRes = runGSA(
            geneLevelStats = myGeneData, # gwas_df_coef, gwas_df_pval
            #directions = gwas_df_coef,
            geneSetStat = myStat,      
            signifMethod = myMethod,     # geneSampling, samplePermutation, nullDist
            adjMethod = myadj,           # bonferroni, none
            gsc = myGsc,                 # gwas_df_trait, 
            gsSizeLim=c(min,Inf),        # gene set size limit
            nPerm=1e4,
            gseaParam=1,
            verbose=F
        )
        #gsaRes %>% print
        
        # save as tsv
        f_name1 = paste0('data/gsea/gsea_',organs[i],'_',myGSC,'_',myStat,'.tsv')
        gsa_df_ = GSAsummaryTable(gsaRes)#, save=T, file=f_name1)
        #cat(paste0(' - file write: ',f_name1),'\n')
        
        # Significance criteria
        if(myStat%in%c('gsea','wilcoxon')) {
            gsa_df_$log_p = -log(gsa_df_$`p (non-dir.)`,10)
            signif = 0.05; sig = paste0('p < ',signif)
            gsa_df_$Significant <- ifelse(gsa_df_$`p (non-dir.)`<signif,sig,"Not Sig")
            mycol = c("grey","red")
        } else {
            gsa_df_$log_p = -log(gsa_df_$`p adj (non-dir.)`,10)
            top = 5; sig = paste0('Top ',top)
            gsa_df_ = gsa_df_[order(gsa_df_$`p adj (non-dir.)`),]
            gsa_df_$Rank = c(rep(sig,top),rep('Not top',nrow(gsa_df_)-top))
            #signif = 1e-50; sig = paste0('bonferroni < ',signif)
            #gsa_df_$Significant <- ifelse(order(gsa_df_$`p adj (non-dir.)`)<top,
            #                              sig,"Not top")
            #mycol = c("red","grey")
            mycol = c("grey","red")
        }
        
        # age.coef mean by TRAIT
        coef_mean = aggregate(gwas_df_coef[,2],list(gwas_df_coef$trait),mean)
        colnames(coef_mean) = c('Name','Age_Coef')
        gsa_df = merge(gsa_df_,coef_mean,by='Name',all.x=T)
        colnames(gsa_df)[2] = 'Gene_n'
        
        # save as tsv file
        write.table(gsa_df,f_name1,quote=F,sep='\t',row.names=F)
        
        # ggplot (x: average aging_coeff, y: significance, circle size: gene_n)
        f_name2 = paste0('data/gsea/gsea_',organs[i],'_',myGSC,'_',myStat,'.png')
        g=ggplot(gsa_df,aes(x=Age_Coef,y=log_p,size=Gene_n))+theme_bw()+
            geom_point(color='Gray')+ #aes(color=Rank)
            labs(y='-log10(p)')+
            #scale_color_manual(values = mycol)+
            geom_vline(xintercept=0,linetype='dashed',color='black')+
            geom_text_repel(
                #data = subset(gsa_df,`-log10(p)`>-log(signif,10)),
                data = gsa_df[order(-gsa_df$log_p),] %>% head(10),
                aes(label = Name),
                size = 3,
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.3, "lines")
            )
        if(myStat%in%c('gsea','wilcoxon')){
            g+geom_hline(yintercept=-log(signif,10),linetype="dashed",color="black")
        } else g
        ggsave(f_name2,width=6,height=5,units='in')
        cat(paste0(' - plot draw: ',f_name2),'\n')
    }
}