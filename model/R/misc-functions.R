load.exp.file <- function(exp.file){
  yaml.load(readChar(exp.file, file.info(exp.file)$size))
}
read.dt.csv <- function(csvfile,key=NULL,...){
    data.table(read.csv(csvfile,...),key=key)
}
read.matrix.csv <- function(csvfile){
  df <- read.csv(csvfile,header=T)
  header.names <- substr(tail(names(df),-1),2,nchar(tail(names(df),-1)))
  row.names <- df[,1]
  mat <- as.matrix(df)[,2:ncol(df)]
  rownames(mat) <- row.names
  colnames(mat) <- header.names
  mat
}
parse.results <- function(soln){
  splitted <- str_split(names(soln),"-")
  dt <- data.table(var=unlist(lapply(splitted,function(ll){ ll[1] })),
                   node=unlist(lapply(splitted,function(ll){ ifelse(ll[1]%in%c('simp','simn'),NA,as.numeric(substr(ll[2],2,nchar(ll[2])))) })),
                   inode=unlist(lapply(splitted,function(ll){ ifelse(ll[1]%in%c('simp','simn'),as.numeric(substr(ll[2],1,str_locate(ll[2],ifelse( str_detect(ll[2],'to'),'to','from'))[,'start']-1)),NA) })),
                   dir=unlist(lapply(splitted,function(ll){ ifelse(ll[1]%in%c('simp','simn'),ifelse( str_detect(ll[2],'to'),'to','from'),NA) })),
                   jnode=unlist(lapply(splitted,function(ll){ ifelse(ll[1]%in%c('simp','simn'),as.numeric(substr(ll[2],str_locate(ll[2],ifelse( str_detect(ll[2],'to'),'to','from'))[,'end']+1,nchar(ll[2]))),NA) })),
                   x=unlist(lapply(splitted,function(ll){ as.numeric(substr(ll[3],2,nchar(ll[3]))) })),
                   t=unlist(lapply(splitted,function(ll){ as.numeric(substr(ll[4],2,nchar(ll[4]))) })),
                   val=soln)
  dt
}
status.code <- function(status){
  if(status==0)return('optimal soln found')
  if(status==1)return('suboptimal')
  if(status==2)return('infeasible')
  if(status==3)return('unbounded')
  if(status==4)return('degenerate')
  if(status==5)return('numerical failure')
  if(status==6)return('user abort')
  if(status==7)return('timeout')
  if(status==9)return('presolved')
  if(status==10)return('procfail')
  if(status==11)return('procbreak')
  if(status==12)return('feasible B&B soln found')
  if(status==13)return('no feasible B&B soln found')
}

