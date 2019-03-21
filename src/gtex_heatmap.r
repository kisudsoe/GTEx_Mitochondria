# Generate heatmap from GTEx gene expression data
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(dendextend))

gtex_heatmap = function(
    gtex_which, gtex_colids, mito_age, # Mendatory DB files
    organ,                             # Organ name is mendatory.
    median=F,                          # Need median summry?
    category='',                       # GO terms/GeneCards description
    category_filt=NULL,                # Category filter
    organ_beta=NULL,                   # Aging coefficient values as organ matrix
    width,height                       # Figure size
) {
    if(is.na(organ)) stop('Please input organ name')
    
    # 1. Extract expression data
    gtex_age_sex = subset(gtex_colids,SMTSD_SYM%in%organ)[,c(1,3,4)]
    gtex_organ = gtex_age_sex$SAMPID
    gtex_age   = gtex_age_sex$AGE
    gtex_sex   = gtex_age_sex$SEX
    gtex_sex[gtex_sex==1] = 'Male'
    gtex_sex[gtex_sex==2] = 'Female'
    gtex_age_sex$SEX = gtex_sex
    
    mito_age_organ = subset(mito_age,Organs%in%organ)
    organ_Ensgid = mito_age_organ$Ensgid_var
    db_source    = mito_age_organ[,c(14,2:4)]
    colnames(db_source)[1] = 'Ensgid'
    rownames(db_source) = organ_Ensgid
    
    #which_row = which(gtex_which$Name%in%organ_Ensgid)
    gtex_which_row = subset(gtex_which,Name%in%organ_Ensgid)
    which_col = which(colnames(gtex_which)%in%gtex_organ)
    gene_exp_which = gtex_which_row[,which_col]#, with=FALSE] # ERROR?!
    Ensgid = unlist(gtex_which_row[,1])        # 1: Ensgid, 2: Symbol
    row.names(gene_exp_which) = Ensgid
    gene_exp_ann1 = cbind(gene_exp_which,Ensgid=Ensgid)
    gene_exp_ann  = as.data.frame(merge(gene_exp_ann1,db_source,by='Ensgid'))
    category_sub = category[,c(3,4,1)] #<- careful!
    colnames(category_sub) = c('Ensgid','Cat','Symbol')
    gene_exp_cat_ = merge(gene_exp_ann,category_sub,by='Ensgid',all.x=T) #<- debug 190214
    
    gene_exp_cat_ = gene_exp_cat_[order(gene_exp_cat_$Cat),]             #<- not work?!
    
    row_names = unlist(gene_exp_cat_$Symbol)
    rownames(gene_exp_cat_) = row_names
    
    # Filt by category
    if(!is.null(category_filt)) {
        gene_exp_cat = subset(gene_exp_cat_,Cat%in%category_filt)
        row_names = unlist(gene_exp_cat$Symbol)
    } else gene_exp_cat = gene_exp_cat_
    
    # Add beta heatmap
    if(!is.null(organ_beta)) {
        m = ncol(organ_beta)-1
        organ_beta_ = organ_beta[,c(2,5:m)]
        colnames(organ_beta_)[1] = 'Ensgid'
        gene_exp_beta = merge(gene_exp_cat,organ_beta_,by='Ensgid')
        #print(colnames(gene_exp_beta))
        
        k = ncol(organ_beta_)-1; n = ncol(gene_exp_beta)
        m1=n-5-k; m2=n-4-k; m3=n-2-k; m4=n-1-k; m5=n-k+1
        gene_exp_  = gene_exp_beta[,c(2:m1)]
        db_source_ = gene_exp_beta[,c(m2:m3)] #; print(colnames(db_source_))
        category_  = gene_exp_beta[,m4]       #; print(head(category_))
        beta_      = gene_exp_beta[,c(m5:n)]  #; print(colnames(beta_))
    } else {
        k = 0; n = ncol(gene_exp_cat)
        m1=n-5-k; m2=n-4-k; m3=n-2-k; m4=n-1-k
        gene_exp_  = gene_exp_cat[,c(2:m1)]
        db_source_ = gene_exp_cat[,c(m2:m3)]
        category_  = gene_exp_cat[,m4]
        #category_[category_==NA] <- 'Others'
    }
    
    #f_name1 = paste0('data/mat_ann_',paste0(organ,collapse='_'),'.tsv')
    #write.table(gene_exp_cat,f_name1,row.names=F,quote=F,sep='\t')
    #cat(paste0('1) Write table: ',f_name1,'\n'))
        
    if(median==TRUE) {
        cat('Median summarization == TRUE\n')
        df0 = as.data.frame(t(gene_exp_))
        df1 = cbind(df0,SAMPID=colnames(gene_exp_)) #row_names
        df2 = merge(df1,gtex_age_sex,by='SAMPID')
        #write.table(df2,'data/gene_exp_m.tsv',row.names=F,quote=F,sep='\t')
        gene_med   = aggregate(df0,list(SEX=df2$SEX,AGE=df2$AGE),median)
        #write.table(gene_med,'data/gene_med.tsv',row.names=F,quote=F,sep='\t')
        gene_exp = as.data.frame(t(gene_med[,c(3:ncol(gene_med))]))
        sex_age = unlist(paste0(gene_med[,1],' - ',gene_med[,2]))
        colnames(gene_exp) = sex_age
    } else gene_exp = gene_exp_
    
    
    # 2. Prepare heatmap 1
    mat_raw = as.matrix(gene_exp); print(dim(mat_raw))
    if(median==TRUE) mat_raw = mat_raw[,order(sex_age)]
    else mat_raw = mat_raw[,order(gtex_sex,gtex_age)]
    #write.table(mat_raw,'data/mat_raw.tsv',quote=F,sep='\t')
    #mat_raw[mat_raw<=1] <- 1 # To avoid NaN error
    mat_log = log(mat_raw,2) # log2 transform
    
    mat.s = t(scale(t(mat_log)))
    row.names(mat.s) = row_names
    #f_name1 = paste0('data/mat_s_',paste0(organ,collapse='_'),'.tsv')
    #write.table(mat.s,f_name1,quote=F,sep='\t')
    
    mat.s_max = max(mat.s,na.rm=T); mat.s_min = min(mat.s,na.rm=T)
    col.rg = c(mat.s_min,0,mat.s_max)
    #col.rg = c(-4,0,4)
    cell.cols = colorRamp2(col.rg,c("Cyan","black","Yellow"))
    
    
    # 3. Prepare top_annotation: Sex-Age
    rb_col   = rainbow(6,start=0.5,end=0.9)
    ha_color = list(Age=c('20-29'=rb_col[1],
                          '30-39'=rb_col[2],
                          '40-49'=rb_col[3],
                          '50-59'=rb_col[4],
                          '60-69'=rb_col[5],
                          '70-79'=rb_col[6]),
                    Sex=c('Male'='dodgerblue','Female'='Pink'))
    if(median==TRUE) {
        ha_ann = data.frame(Sex=c(rep('Female',6),rep('Male',6)),
                            Age=c('20-29','30-39','40-49','50-59','60-69','70-79',
                                  '20-29','30-39','40-49','50-59','60-69','70-79'))
        ha = HeatmapAnnotation(df=ha_ann,col=ha_color)
    } else {
        col_names = data.frame(SAMPID=colnames(mat.s))
        ha_df = merge(col_names,gtex_age_sex,by='SAMPID')
        ha_df = ha_df[order(gtex_sex,gtex_age),]
        ha_ann = data.frame(Age=ha_df$AGE,Sex=ha_df$SEX)
        ha = HeatmapAnnotation(df=ha_ann,col=ha_color)
    }
    
    
    # 4. Prepare heatmap 2
    db_source_ = as.matrix(db_source_)
    db_source_[db_source_== TRUE] <- 'Exist'
    db_source_[db_source_==FALSE] <- 'Non'
    rownames(db_source_) = row_names
    binary_cols = c(`Exist`='black',`Non`='gray')
        
    
    # 5. Draw heatmaps
    ht1=Heatmap(
        mat.s,
        top_annotation=ha,
        column_title=organ,
        #name='GTEx Gene\nMedian TPM\n(Z-scaled)',
        name='Gene\nexpression',
        cluster_columns=F,
        #show_column_names=F,
        #column_names_side='top',
        #cluster_rows=F,
        split=category_,
        #show_row_names=F,
        row_dend_reorder=F,
        #row_dend_width=unit(1.5,'in'),
        row_names_side='left',
        rect_gp=gpar(col='black'),
        col=cell.cols)
    
    ht2=Heatmap(
        db_source_,
        column_title='DBs',
        name='DB exists',
        cluster_column=F,
        cluster_row=F,
        show_row_names=F,
        row_dend_reorder=F,
        #column_names_side='top',
        rect_gp=gpar(col='white'),
        col=binary_cols)
    
    # 6. Draw heatmap 3
    if(!is.null(organ_beta)) {
        beta_[sapply(beta_,function(x) all(is.na(x)))] <- NULL # remove NA only columns
        beta = as.matrix(beta_)
        col_na = apply(beta,2,function(x) sum(is.na(x)))
        col_na_order = order(col_na)
        beta = beta[,col_na_order]
        
        beta_max = max(beta,na.rm=T); beta_min = min(beta,na.rm=T)
        col.rg  = c(beta_min,0,beta_max)
        cell.cols2 = colorRamp2(col.rg,c('deepskyblue','white','indianred1'))
        
        ht3=Heatmap(
            beta,
            name='Beta',
            cluster_column=F,
            column_title='Tissues',
            cluster_row=F,
            #row_title='Gerontic factor',
            show_row_names=F,
            row_dend_reorder=F,
            rect_gp=gpar(col='white'),
            col=cell.cols2)
        
        draw(ht1+ht2+ht3) #,heatmap_legend_side='left'
        width=ncol(beta)*0.1+10; height=nrow(beta)*0.1+7 # True length version
        #width=ncol(beta)*0.1+10; height=15 # for compact version
    } else {
        draw(ht1+ht2,heatmap_legend_side='left')
        width=15; height=nrow(beta)*0.15+5
    }
    if(!is.null(category_filt)) {
        f_name2 = paste0('fig/',paste0(organ,collapse='_'),'_',
                         paste0(category_filt,collapse='_'),'.png')
    } else f_name2 = paste0('fig/',paste0(organ,collapse='_'),'.png')
    dev.copy(png,f_name2,width=width,height=height,units='in',res=500)
    cat(paste0('2) Draw plot:   ',f_name2,'\n'))
    #while(!is.null(dev.list())) dev.off()
    graphics.off() # killing all devices
    return()
}

#cell.cols = colorRamp2(col.rg,c("black", "cornflowerblue", "yellow", "red"))
#cell.cols = colorRamp2(col.rg,c("black", "cornflowerblue", "yellow"))

#hr = hclust(dist(mat.s),method="average")
#clu = as.data.frame(cutree(hr,h=8))
#clu[clu==3] <- 'C1'
#clu[clu==8] <- 'C2'