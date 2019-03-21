suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
data_summary = function(
	organ_df,
	f_name
) {
	df = organ_df$Organ %>% as.character %>% table %>% data.frame
	colnames(df)=c('Organs','Gene_n')
	df$Organs = factor(
		df$Organs,
		levels=df$Organs[order(-df$Gene_n)] # descending order by Gene_n
	)
	
	age_df_up = subset(organ_df,Age.Coef>0)
	up_age_n = age_df_up$Organs %>% table %>% data.frame
	colnames(up_age_n)=c('Organs','Age_up')
	df_1 = merge(df,up_age_n,by='Organs',all.x=T)

	age_df_dn = subset(organ_df,Age.Coef<0)
	dn_age_n = age_df_dn$Organs %>% table %>% data.frame
	colnames(dn_age_n)=c('Organs','Age_down')
	df_2 = merge(df_1,dn_age_n,by='Organs',all.x=T)
	df_2$Organs = factor(
		df_2$Organs,
		levels=df_2$Organs[order(-df_2$Gene_n)] # descending order by Gene_n
	)
	
	df_3 = melt(df_2[,c(1,3,4)],id='Organs')
	colnames(df_3)=c('Organs','Direction','Gene_n')
	df_3$Direction = factor(
		df_3$Direction,
		levels=c('Age_down','Age_up')
	)
	df_ = ddply(
		df_3,
		c('Organs'),
		transform,
		label_y=cumsum(Gene_n)
	)
	f_name1 = paste0('data/',f_name,'.tsv')
	write.table(df_,f_name1,row.names=F,quote=F,sep='\t')
	paste0('File write: ',f_name1,'\n') %>% cat
	
	# Draw plot as SVG format
	p=ggplot(df_,aes(x=Organs,y=Gene_n,fill=Direction))+theme_bw()+
		geom_bar(stat='identity')+
		scale_fill_manual(values=c('cyan3','coral1'))+
		geom_text(aes(y=label_y,label=Gene_n),
				  hjust=0.3,vjust=.5,angle=45)+
		theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
	f_name2 = paste0('fig/',f_name,'.svg')
	ggsave(f_name2,p,device='svg',width=12,height=8)
	paste0('Plot write: ',f_name2,'\n') %>% cat
	return(p)
}