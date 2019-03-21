# Read file and save as rds file
suppressMessages(library(data.table))
suppressMessages(library(tools))

saveasrds = function(f_paths) {
    n = length(f_paths); df.li = list()
    df.li=lapply(c(1:n),function(i) {
        if(file_ext(f_paths[i])=='gz') { # if file extension is '.gz'
            try(R.utils::gunzip(f_paths[i]))
            f_path = file_path_sans_ext(f_paths[i])
            return(fread(f_path,sep='\t',header=T,stringsAsFactors=F))
        } else {
            return(fread(f_paths[i],sep='\t',header=T,stringsAsFactors=F))
        }
    })
    df = rbindlist(df.li) # data.table
    f_name = paste0(f_paths[1],'.rds')
    saveRDS(df,file=f_name)
    cat(paste0('>> File write: ',f_name,'\n'))
}

suppressMessages(library(rtracklayer))

bedasrds = function(f_path) {
    df = as.data.frame(import(f_path,format='bed'))
    f_name = paste0(f_path,'.rds')
    saveRDS(df,file=f_name)
    cat(paste0('>> File write: ',f_name,'\n'))
}