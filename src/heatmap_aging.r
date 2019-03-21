# Heatmap using Aging coefficient values
suppressMessages(library(dplyr))
suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap))
heatmap_aging = function(
	organ_ann,
	filter_n=1,
	over=TRUE,
	f_name,
	width=10,
	height=10
) {
	# 1. Count organ_n of genes
	organs = as.character(unique(organ_ann$Organs))
	cat(paste0('Organs_n = ',length(organs),'\n'))
	organ_ann_df = organ_ann[,c(1:4,7:8)] %>%
		spread(Organs,Age.Coef)
	organ_ann_df = organ_ann_df[order(-organ_ann_df$Organ_n),]
	
	# Add organ list tag
	organ_ann_df$organ_list = unlist(apply(organ_ann_df,1,function(row) {
		row_ = row[c(5:length(row))]
		row_ = row_[!is.na(row_)]
		return(paste0(names(row_),collapse=", "))
	}))
	
	# Column sort by hits
	#organ_major = c(
	#	'Artery_Tibial','Artery_Aorta','Nerve_Tibial','Brain_Cortex',
	#	'Uterus','Ovary',
	#	'Adipose_Subcutaneous','Muscle_Skeletal','Esophagus_Mucosa'
	#)
	#organ_major_cols = which(colnames(organ_ann_df)%in%organ_major)
	
	# Row sort by directions
	for(i in 31:5) {
		organ_ann_df = organ_ann_df[order(-organ_ann_df[,i]),]
	}
	organ_ann_df = organ_ann_df[order(-organ_ann_df$Organ_n),]
	
	# Save heatmap data as tsv file
	f_name1 = paste0('data/',f_name,'_heatmap.tsv')
	write.table(organ_ann_df,f_name1,sep='\t',row.names=F)
	cat(paste0('File write: ',f_name1,'\n'))
	
	# 2. Filtering genes by their shared organ numbers
	if(over==TRUE) organ_ann_df_subset = subset(organ_ann_df,Organ_n>=filter_n)
	else organ_ann_df_subset = subset(organ_ann_df,Organ_n==filter_n)
	
	if(filter_n==1) {
		# Column reorder by beta hits
		col_na = apply(organ_ann_df_subset[,c(5:31)],2,function(x) sum(is.na(x)))
		col_na_order = order(col_na)+4
		organ_ann_df_subset = organ_ann_df_subset[,c(1:4,col_na_order,32)]
		
		# Row sort by directions
		for(i in 31:5) {
			organ_ann_df_subset = organ_ann_df_subset[order(-organ_ann_df_subset[,i]),]
		}
		#organ_ann_df_subset = organ_ann_df_subset[order(-organ_ann_df_1$Organ_n),]
	} else {
		# Row reorder by directions
		row_which_as = organ_ann_df_subset[,5]
		row_which_as[which(organ_ann_df_subset[,5]>0)] <- 1
		row_which_as[which(organ_ann_df_subset[,5]<0)] <- 2
		row_which_as[which(organ_ann_df_subset[,5]%>%is.na)] <- 3
		row_which = apply(organ_ann_df_subset[,5:27],1,function(row) {
			up_n = length(which(row>0))
			dn_n = length(which(row<0))
			if(up_n>0&dn_n==0) return(1)
			else if(up_n==0&dn_n>0) return(2)
			else if(up_n>0 &dn_n>0) return(3)
			else if(up_n==0&dn_n==0) return(NA)
		}) %>% unlist
		organ_ann_df_subset = organ_ann_df_subset[order(row_which_as,row_which),]
	}
	# Row split
	if(over==T) {
		organ_split = factor(organ_ann_df_subset$Organ_n)
		organ_split_levels = levels(organ_split)
		organ_split = factor(
			organ_split,
			levels=organ_split_levels[length(organ_split_levels):1])
		m = ncol(organ_ann_df_subset)-1
	} else {
		m = ncol(organ_ann_df_subset)-1
		organ_split = factor(
			organ_ann_df_subset$organ_list,
			levels=colnames(organ_ann_df_subset)[5:m]
		)
	}
	df = organ_ann_df_subset[,c(5:m)]
	df = df[,colSums(is.na(df))<nrow(df)]
	mat = as.matrix(df)
	
	# Column reorder by beta hits
	col_na = apply(mat[,c(2:ncol(mat))],2,function(x) sum(is.na(x)))
	col_na_order = order(col_na)+1
	mat = mat[,c(1,col_na_order)]
	
	# Set rownames
	row.names(mat) = organ_ann_df_subset$GeneSymbol # 1: Ensgid, 2: Ensgid_var, 3: GeneSymbol
	paste0('Heatmap dimension: ',dim(mat)[1],', ',dim(mat)[2])
	
	# Draw Heatmap
	mat_max = max(mat,na.rm=T); mat_min = min(mat,na.rm=T)
	col.rg  = c(mat_min,0,mat_max)
	cell.cols = colorRamp2(col.rg,c('deepskyblue','white','indianred1'))
	h = Heatmap(
		mat,
		name='Beta',
		cluster_column=F,
		cluster_row=F,
		split=organ_split,
		row_title='Genes',
		column_title='Tissues',
		row_names_side='left',
		rect_gp=gpar(col='black'),
		column_names_max_height=unit(5,'in'),
		col=cell.cols
	)
	draw(h,heatmap_legend_side='left')
	if(over==T) f_name1 = paste0('fig/',f_name,'_',filter_n,'_over.svg')
	else f_name1 = paste0('fig/',f_name,'_',filter_n,'.svg')
	dev.copy(svg,f_name1,width=width,height=height)
	paste0('Plot write: ',f_name1,'\n') %>% cat
	graphics.off() # killing all devices
}