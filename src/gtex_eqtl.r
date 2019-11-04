suppressMessages(library(dplyr))
gtex_eqtl_genotype_p = function(
    gtex_exp_df      = NULL, # DB
    gtex_genotype_df = NULL, # DB
    eqtl_rsid        = NULL, # variantID-Rsid
    rsid_gene        = NULL, # Rsid-GeneSymbol
    statistics		 = '',   # c('anova','lm')
    f_write          = T
) {
    rsid_gene = rsid_gene%>%unique
    n = nrow(rsid_gene)
    paste0('eQTL-gene pairs: ',n,'\n')%>%cat
    #dat_rows = gtex_exp_df[,c(1:3,4,6,10)] #<- SAMPID,SUBJID,SEX,AGE,AGE_db,SMTSD_SYM
    
    # If Statistics == anova,
    if(statistics=='anova') {
    	paste0('\nCalculating anova p by genotypes\n')%>%cat
    	geno_p_li = lapply(c(1:n),function(i) {
			rsid = rsid_gene$Rsid[i]
			variant_id_which = which(eqtl_rsid$RSID %in% rsid)
			variant_id = eqtl_rsid[variant_id_which,1]%>%as.character
			
			gene_name = rsid_gene$GeneSymbol[i]
			colid = which(colnames(gtex_exp_df) %in% gene_name)
			
			var_which = which(colnames(gtex_genotype_df) %in% variant_id)
			if(var_which%>%length==0) {
				paste0('* ',gene_name,'-',rsid,' -> ',variant_id,
					' -> [ERROR] No found variant. Check Tabix code.\n'
					#'   Please revise the Tabix code from the dbGaP GTEx data.\n\n'
				)%>%cat
				anova_df = data.frame(
					organ       = NA,
					rsid        = rsid,
					gene        = gene_name,
					gt_anova_p  = NA,
					gt_male_p   = NA,
					gt_female_p = NA
				)
			} else if(var_which%>%length==1) {
				gtex_geno_df_1 = gtex_genotype_df[,c(1,var_which)]
				gtex_exp_df_1 = merge(gtex_exp_df,gtex_geno_df_1,by="SUBJID",all.x=T)
				colnames(gtex_exp_df_1)[ncol(gtex_exp_df_1)] = 'eQTL'

				organs = levels(gtex_exp_df_1$SMTSD_SYM)
				anova_df=NULL
				for(j in 1:length(organs)) {
					dat_organ = subset(gtex_exp_df_1, SMTSD_SYM %in% organs[j])
					#dat_organ$eQTL%>%levels%>%print
					a = dat_organ[,colid]%>%unlist%>%as.numeric
					b = dat_organ$eQTL%>%unlist%>%as.character
					tryCatch({ #try body
						anova = aov(a~b)
						print(summary(anova))
						anova_p = summary(anova)[[1]][["Pr(>F)"]][[1]]
					}, error=function(cond) {
						anova_p = NA
					})
					
					# split by SEX - Male
					dat_organ_m = subset(dat_organ,SEX=='Male'); #dim(dat_organ_m)%>%print
					a_m = dat_organ_m[,colid]%>%unlist%>%as.numeric
					b_m = dat_organ_m$eQTL%>%unlist%>%as.character
					tryCatch({ #try body
						anova = aov(a_m~b_m)
						anova_m_p = summary(anova)[[1]][["Pr(>F)"]][[1]]
					}, error=function(cond) {
						anova_m_p = NA
					})
					
					# split by SEX - Female
					dat_organ_f = subset(dat_organ,SEX=='Female'); #dim(dat_organ_f)%>%print
					a_f = dat_organ_f[,colid]%>%unlist%>%as.numeric
					b_f = dat_organ_f$eQTL%>%unlist%>%as.character
					tryCatch({ #try body
						anova = aov(a_f~b_f)
						anova_f_p = summary(anova)[[1]][["Pr(>F)"]][[1]]
					}, error=function(cond) {
						anova_f_p = NA
					})
					
					# summary table
					anova_ = data.frame(
						organ       = organs[j],
						rsid        = rsid,
						gene        = gene_name,
						gt_anova_p  = anova_p,
						gt_male_p   = anova_m_p,
						gt_female_p = anova_f_p
					)
					if(j==1) { anova_df = anova_
					} else anova_df = rbind(anova_df,anova_)
				} # end for
			} else '\n\n-> Matched eQTL ID is more than 2.\n'%>%cat
			return(anova_df)
		})
    }
    
    # If Statistics == lm (linear model),
    if(statistics=='lm') {
    	paste0('\nCalculating linear regression p\n')%>%cat
    	geno_p_li = lapply(c(1:n),function(i) {
			rsid = rsid_gene$Rsid[i]
			variant_id_which = which(eqtl_rsid$RSID %in% rsid)
			variant_id = eqtl_rsid[variant_id_which,1]%>%as.character
			
			gene_name = rsid_gene$GeneSymbol[i]
			colid = which(colnames(gtex_exp_df) %in% gene_name)
			
			var_which = which(colnames(gtex_genotype_df) %in% variant_id)
			if(var_which%>%length==0) {
				paste0('* ',gene_name,'-',rsid,' -> ',variant_id,
					' -> [ERROR] No found variant. Check Tabix code.\n'
					#'   Please revise the Tabix code from the dbGaP GTEx data.\n\n'
				)%>%cat
				anova_df = data.frame(
					organ   = NA,
					rsid    = rsid,
					gene    = gene_name,
					gt_lm_p = NA
				)
			} else if(var_which%>%length==1) {
				gtex_geno_df_1 = gtex_genotype_df[,c(1,var_which)]
				gtex_exp_df_1 = merge(gtex_exp_df,gtex_geno_df_1,by="SUBJID",all.x=T)
				colnames(gtex_exp_df_1)[ncol(gtex_exp_df_1)] = 'eQTL'
				
				organs = levels(gtex_exp_df_1$SMTSD_SYM)
				anova_df=NULL
				for(j in 1:length(organs)) {
					dat_organ = subset(gtex_exp_df_1, SMTSD_SYM %in% organs[j])
					exp  = dat_organ[,colid]%>%unlist%>%as.numeric # expression
					geno = dat_organ$eQTL%>%unlist%>%as.character # genotype
					geno%>%table%>%print
					gender=dat_organ$SEX%>%unlist%>%as.character # gender
					gender%>%table%>%print
					age  = dat_organ$AGE_db%>%unlist%>%as.numeric # age
					age%>%summary%>%print
					tryCatch({ #try body
						fit = lm(exp~geno+gender+age+geno*age+gender*age+geno*gender)
						summary(fit)$coefficients %>% print
						summary(fit)$r.squared %>% print
						lm_p = summary(fit)[[1]][["Pr(>F)"]][[1]]
					}, error=function(cond) {
						lm_p = NA
					})
					
					# summary table
					anova_ = data.frame(
						organ    = organs[j],
						rsid     = rsid,
						gene     = gene_name,
						gt_lm_p  = lm_p
					)
					if(j==1) { anova_df = anova_
					} else anova_df = rbind(anova_df,anova_)
				} # end for
			} else '\n\n-> Matched eQTL ID is more than 2.\n'%>%cat
			return(anova_df)
		})
    }
    
    gt_p_df = plyr::ldply(geno_p_li,data.frame) #plyr
    #gt_p_df = Reduce(function(x,y) merge(x,y,by="SUBJID",all=T), geno_li)
    'Iteration done!\n'%>%cat
    
    if(f_write==T) {
        f_name = paste0('data/gtex_eqtl_anova.tsv')
        write.table(gt_p_df,f_name,sep='\t',row.names=F)
        paste0('file write done: ',f_name,'\n')%>%cat
    }
    return(gt_p_df)
}

gtex_genotype_table = function(
    gtex_genotype = NULL,
    eqtl_rsid = NULL
) {
    n = gtex_genotype%>%nrow; #n%>%print
    geno_li = lapply(c(1:n),function(i) { #4544:4548
        row = gtex_genotype[i,]; #row%>%print
        Ref = paste0('Ref_',row$REF,'/',row$REF); #Ref%>%print
        Alt = paste0('Alt_',row$ALT,'/',row$ALT); #Alt%>%print
        Het = paste0('Het_',row$REF,'/',row$ALT); #Het%>%print

        Ref_subjid = strsplit(row$HomSamplesRef%>%as.character,'\\,')[[1]]
        Ref_n = Ref_subjid%>%length
        Alt_subjid = strsplit(row$HomSamplesAlt%>%as.character,'\\,')[[1]]
        Alt_n = Alt_subjid%>%length
        Het_subjid = strsplit(row$HetSamples   %>%as.character,'\\,')[[1]]
        Het_n = Het_subjid%>%length
        #Ref_df%>%length%>%print; Ref_df%>%head%>%print
        genotype_df = data.frame(SUBJID=Ref_subjid,GENOTYPE=rep(Ref,Ref_n))
        genotype_df = rbind(
            genotype_df,
            data.frame(SUBJID=Alt_subjid,GENOTYPE=rep(Alt,Alt_n)))
        genotype_df = rbind(
            genotype_df,
            data.frame(SUBJID=Het_subjid,GENOTYPE=rep(Het,Het_n)))
        #row$RSID%>%as.character%>%print
        colnames(genotype_df) = c('SUBJID',row$ID%>%as.character)
        #genotype_df%>%dim%>%print
        #genotype_df%>%head%>%print
        return(genotype_df)
    })
    #geno_df = plyr::ldply(geno_li,data.frame) #plyr
    geno_df = Reduce(function(x,y) merge(x,y,by="SUBJID",all=T), geno_li)
    return(geno_df)
}

gtex_eqtl_gene_table = function(
    gtex_exp_df      = NULL,
    gtex_genotype_df = NULL,
    eqtl_rsid        = NULL,
    gene_name        = '',
    eqtl_id          = '',
    f_write          = F
) {
    paste0('Gene: ',gene_name,' and eQTL ID: ',eqtl_id,'\n')%>%cat
    variant_id_which = which(eqtl_rsid$RSID%in%eqtl_id)
    variant_id = eqtl_rsid[variant_id_which,1]%>%as.character
    paste0('Searched GTEx variant ID: ',variant_id,'\n')%>%cat
    
    paste0('\n1. Extract columns')%>%cat
    dat_rows = gtex_exp_df[,c(1:3,4,6,10)] #<- SAMPID,SUBJID,SEX,AGE,AGE_db,SMTSD_SYM
    colid    = which(colnames(gtex_exp_df)%in%gene_name)
    dat_val  = gtex_exp_df[,colid]%>%as.numeric #<-gene_name
    #colnames(dat_rows)%>%print
    '...done!\n'%>%cat
    
    paste0('2. Combine variant genotypes')%>%cat
    var_which = which(colnames(gtex_genotype_df)%in%variant_id)
    gtex_geno_df_1 = gtex_genotype_df[,c(1,var_which)]
    if(length(var_which)==0) {
        '\n\n[STOP] No found variant. Please revise the Tabix code from the dbGaP GTEx data.'%>%cat
        return() }
    dat_1 = merge(dat_rows,gtex_geno_df_1,by="SUBJID",all.x=T)
    colnames(dat_1)[ncol(dat_1)] = 'eQTL'
    #dat_1%>%dim%>%print; dat_1%>%head(3)%>%print
    
    gtex_exp_df_1 = merge(gtex_exp_df,gtex_geno_df_1,by="SUBJID",all.x=T)
    colnames(gtex_exp_df_1)[ncol(gtex_exp_df_1)] = 'eQTL'
    #gtex_exp_df_1%>%dim%>%print; gtex_exp_df_1%>%head(3)%>%print
    #gtex_exp_df_1%>%colnames%>%print
    '...done!\n'%>%cat
    
    paste0('3. Calculate correlations r\n')%>%cat
    cor_test = cor.test(
        gtex_exp_df$AGE_db,
        as.numeric(gtex_exp_df[,colid]))
    r_p = c(cor_test$estimate,cor_test$p.value)
    names(r_p) = c('estimate','p.value')
    
    gtex_exp_df_m = subset(gtex_exp_df,SEX=='Male')
    gtex_exp_df_f = subset(gtex_exp_df,SEX=='Female')
    cor_m = cor.test(
        gtex_exp_df_m$AGE_db,
        as.numeric(gtex_exp_df_m[,colid]))
    rp_m = c(cor_m$estimate,cor_m$p.value)
    cor_f = cor.test(
        gtex_exp_df_f$AGE_db,
        as.numeric(gtex_exp_df_f[,colid]))
    rp_f = c(cor_f$estimate,cor_f$p.value)
    names(rp_m) = c('Male_r','Male_p')
    names(rp_f) = c('Female_r','Female_p')
    
    paste0(' - Calculate correlations r by organs')%>%cat
    organs = levels(gtex_exp_df$SMTSD_SYM)
    dat_s=NULL; cor_df=NULL
    for(i in 1:length(organs)) {
        dat_organ = subset(gtex_exp_df_1,SMTSD_SYM%in%organs[i])
        geno = dat_organ$eQTL%>%levels; #print(geno)
        #dat_ref = subset(dat_organ,eQTL==geno[1])
        #dat_alt = subset(dat_organ,eQTL==geno[2])
        #dat_het = subset(dat_organ,eQTL==geno[3])
        a = dat_organ[,colid]%>%unlist%>%as.numeric
        b = dat_organ$eQTL%>%unlist%>%as.character
        anova = aov(a~b); #anova%>%summary%>%print #<-ERROR! if()
        anova_p = summary(anova)[[1]][["Pr(>F)"]][[1]]; #anova_p%>%print
        
        dat_age_  = dat_organ$AGE_db
        dat_val_  = dat_organ[,colid] #<-gene_name
        dat_s_    = scale(as.numeric(unlist(dat_val_)))
        dat_s     = c(dat_s,dat_s_)
        cor_test  = cor.test(as.numeric(dat_age_),as.numeric(dat_val_))
        
        dat_org_m = subset(dat_organ,SEX=='Male')
        if(nrow(dat_org_m)==0) {
            cor_r_male=NA
            cor_p_male=NA
        } else cor_test_m = cor.test(
            as.numeric(dat_org_m$AGE_db),
            as.numeric(dat_org_m[,colid])) #<-gene_name
        dat_org_f = subset(dat_organ,SEX=='Female')
        if(nrow(dat_org_f)==0) {
            cor_r_female=NA
            cor_p_female=NA
        } else cor_test_f = cor.test(
            as.numeric(dat_org_f$AGE_db),
            as.numeric(dat_org_f[,colid])) #<-gene_name
        
        cor_ = data.frame(
            organ=organs[i],
            cor_r=cor_test$estimate,
            cor_p=cor_test$p.value,
            cor_r_male=cor_test_m$estimate,
            cor_p_male=cor_test_m$p.value,
            cor_r_female=cor_test_f$estimate,
            cor_p_female=cor_test_f$p.value,
            gt_anova_p=anova_p
        )
        if(i==1) { cor_df = cor_
        } else cor_df = rbind(cor_df,cor_)
    }
    '...done!\n'%>%cat
    
    paste0('4. Summary tables as list')%>%cat
    dat = cbind(dat_1,value.s=dat_s,value=dat_val)
    
    #cor_df[,c(1,3,8)]%>%head(5)%>%print
    cor_gene_ = subset(
        cor_df, gt_anova_p<1e-4)#cor_p<1e-3
    cor_gene = cor_gene_[order(cor_gene_$gt_anova_p),] #cor_gene_$cor_p
    #cor_gene%>%dim%>%print
    
    if(f_write==T) {
        f_name = paste0('data/gene_',gene_name,'_cor.tsv')
        write.table(cor_gene,f_name,sep='\t',row.names=F)
        print(paste0('file write done: ',f_name))
    }
    '...done!\n'%>%cat
    return(list(dat, cor_gene, r_p, rp_m, rp_f))
}

suppressMessages(library(ggplot2))
gtex_eqtl_gene_plot = function(
    dat_wo_na = NULL,
    dat_list  = NULL,
    gene_name = '',
    eqtl_id   = '',
    plot_type = '' #c('all','violin','dot')
) {
    cor_rp1 = dat_list[[2]][1,]; #print(cor_rp1)
    organ = cor_rp1$organ
    cor_rp1 = unlist(cor_rp1)
    paste0('  a. Filter for ',organ,'\n')%>%cat
    #dat = subset(dat_list[[1]],SMTSD_SYM%in%organ); #dat%>%dim%>%print
    dat = subset(dat_wo_na,SMTSD_SYM%in%organ)
    r   = cor_rp1[2]; p   = cor_rp1[3]
    r_m = cor_rp1[4]; p_m = cor_rp1[5]
    r_f = cor_rp1[6]; p_f = cor_rp1[7]
    anova_p = cor_rp1[8]
    
    paste0('  b. Draw plots','\n\n')%>%cat
    tb = doBy::summaryBy(value~SEX+eQTL_f,dat,FUN=length)%>%t
    if(plot_type=='all'||'violin') {
        p1 = ggplot(dat,aes(x=AGE,y=value.s))+#theme_bw()+
            geom_violin(trim=F)+
            stat_summary(
                fun.data=mean_sdl,
                fun.args=list(mult=1),
                geom='crossbar',
                width=0)+
            stat_summary(
                fun.y=median,
                fun.ymin=median,
                fun.ymax=median,
                geom='crossbar',
                width=0.4)+
            labs(title=paste0(
                gene_name,'(',organ,');\n  r= ',r%>%round(3),', p = ',p%>%round(4)
            ),y='Gene expression')
        p2 = ggplot(dat,aes(x=AGE,y=value.s))+#theme_bw()+
            geom_violin(aes(colour=SEX),trim=F)+
            stat_summary(
                aes(colour=SEX),
                fun.data=mean_sdl,
                fun.args=list(mult=1),
                geom='crossbar',
                width=0,
                position=position_dodge(0.9))+
            stat_summary(
                aes(colour=SEX),
                fun.y=median,
                fun.ymin=median,
                fun.ymax=median,
                geom='crossbar',
                width=0.4,
                position=position_dodge(0.9))+
            labs(title=paste0(
                gene_name,'(',organ,');\n  Male (r=',r_m%>%round(3),', p =',p_m%>%round(4),')\n',
                '  Female (r=',r_f%>%round(3),', p =',p_f%>%round(4),')'
            ),y='Gene expression')
    }
    if(plot_type=='all'||'dot') {
        p3 = ggplot(dat,aes(x=AGE_db,y=value.s))+theme_bw()+#ylim(-4,4)
            geom_point(alpha=0.5,size=2)+geom_rug()+
            geom_density_2d()+
            geom_smooth(method=lm)+ #
            facet_grid(.~eQTL_f)+
            labs(title=paste0(
                eqtl_id,'-',gene_name,' (',organ,'); r= ',r%>%round(3),', anova p = ',anova_p%>%round(4)
            ),x='AGE',y='Gene expression')
        p4 = ggplot(dat,aes(x=AGE_db,y=value.s))+theme_bw()+#ylim(-4,4)+#colour=SEX
            geom_point(alpha=0.5,size=2)+geom_rug()+
            geom_density_2d()+
            geom_smooth(method=lm)+ #
            facet_grid(SEX~eQTL_f)+
            labs(title=paste0(
                eqtl_id,'-',gene_name,' (',organ,');\n  ',
                'Male (r=',r_m%>%round(3),', p =',p_m%>%round(4),') and ',
                'Female (r=',r_f%>%round(3),', p =',p_f%>%round(4),')'
            ),x='AGE',y='Gene expression')+
            theme(legend.position='bottom',legend.box='horisontal')
        
        p5 = ggplot(dat,aes(x=AGE_db,y=value.s,colour=eQTL_f))+theme_bw()+#ylim(-4,4)+
            geom_point(alpha=0.5,size=2)+geom_rug()+
            geom_smooth(aes(fill=eQTL_f),method=lm,fullragne=T)+ #
            labs(title=paste0(
                eqtl_id,'-',gene_name,' (',organ,');\n  ',
                'r= ',r%>%round(3),', anova p = ',anova_p%>%round(4)
            ),x='AGE',y='Gene expression')+
            theme(legend.position='bottom',legend.box='horisontal')
        p6 = ggplot(dat,aes(x=AGE_db,y=value.s,colour=SEX))+theme_bw()+#ylim(-4,4)+
            geom_point(alpha=0.5,size=2)+geom_rug()+
            geom_smooth(aes(fill=SEX),method=lm,fullragne=T)+ #
            facet_grid(.~eQTL_f)+
            labs(title=paste0(
                eqtl_id,'-',gene_name,' (',organ,');\n  ',
                'Male (r=',r_m%>%round(3),', p =',p_m%>%round(4),') and ',
                'Female (r=',r_f%>%round(3),', p =',p_f%>%round(4),')'
            ),x='AGE',y='Gene expression')+
            theme(legend.position='bottom',legend.box='horisontal')
        
        p7 = ggplot(dat,aes(x=AGE_db,y=value.s,colour=eQTL_f))+theme_bw()+#ylim(-4,4)+
            geom_point(alpha=0.5,size=2)+geom_rug()+
            geom_smooth(aes(fill=eQTL_f),method=lm,fullragne=T)+ #
            facet_grid(.~SEX)+
            labs(title=paste0(
                eqtl_id,'-',gene_name,' (',organ,');\n  ',
                'r= ',r%>%round(3),', anova p = ',anova_p%>%round(4)
            ),x='AGE',y='Gene expression')+
            theme(legend.position='bottom',legend.box='horisontal')
    }
    return(list(p1,p2,p3,p4,p5,p6,p7, tb))
}

gtex_eqtl_plot = function(
    gtex_exp_df      = NULL,
    gtex_genotype_df = NULL,
    eqtl_rsid        = NULL,
    gene_name        = NULL,
    eqtl_id          = NULL,
    f_write          = F
) {
    # Process 1-3.
    gene_df = gtex_eqtl_gene_table(
        gtex_exp_df      = gtex_exp_df,
        gtex_genotype_df = gtex_genotype_df,
        eqtl_rsid        = eqtl_rsid,
        gene_name        = gene_name,
        eqtl_id          = eqtl_id,
        f_write          = F
    )
    if(length(gene_df)==0) return()
    
    '\n-> Genotype significant tissue list is generated.\n'%>%cat
    if(nrow(gene_df[[2]])==0) {
        '\n[STOP] No significant difference between genotype.\n'%>%cat
        return() }
    gene_df_1 = subset(gene_df[[2]],cor_p<0.1)
    gene_df_1%>%dim%>%print
    gene_df_1[,c(1,3,8)]%>%print
    
    cat('\n'); gene_df[[1]]%>%dim%>%print
    eqtl_na_which = which(!is.na(gene_df[[1]]$eQTL))
    gene_df_1 = gene_df[[1]][eqtl_na_which,]
    gene_df_1%>%dim%>%print
    '-> Eliminating rows with genotype NA is done.\n\n'%>%cat
    eQTL_lv = levels(gene_df_1$eQTL)[c(1,3,2)]; #print(eQTL_lv)
    gene_df_1$eQTL_f = factor(gene_df_1$eQTL,levels=eQTL_lv)
    eqtl_f = gene_df_1$eQTL_f%>%summary%>%print
    
    paste0('\n4. Draw plots:','\n')%>%cat
    pl = gtex_eqtl_gene_plot(
        dat_wo_na = gene_df_1,
        dat_list  = gene_df,
        gene_name = gene_name,
        eqtl_id   = eqtl_id,
        plot_type = 'all') #c('all','violin','dot')
    
    return(pl)
}