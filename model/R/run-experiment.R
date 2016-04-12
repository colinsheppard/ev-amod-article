load.libraries(c('optparse','yaml','lpSolve','lpSolveAPI'))
source(pp(ev.amod.model,'model/R/misc-functions.R'))

option_list <- list(
   make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print extra output [%default]")
)
if(interactive()){
  args<-'C:\\Users\\bsj12_000\\Google Drive\\CE295/ev-amod-outtages/model/inputs/exp.yaml'
  args<-'C:\\Users\\bsj12_000\\Google Drive\\CE295/ev-amod-outtages/model/inputs/tiny/tiny.yaml'
  args <- parse_args(OptionParser(option_list = option_list,usage = "autonomous-sim.R [experiment YAML file]"),positional_arguments=T,args=args)
  exp.dirs <- args$args
}else{
  args <- parse_args(OptionParser(option_list = option_list,usage = "autonomous-sim.R [experiment YAML file]"),positional_arguments=T)
  exp.dirs <- args$args
}
verbose <- args$options$verbose
exp.file <- args$args[1]

exp <- load.exp.file(exp.file)
factors <- unlist(lapply(exp$Factors,function(ll){ streval(pp('list(',ll$Code,'=c("',pp(unlist(lapply(ll$Levels,function(lll){ lll$Code })),collapse='","'),'"))')) }),recursive=F)
factor.grid <- expand.grid(factors)

# arrange parameter values and convert to numeric data types where appropriate
params.base <- data.frame(t(unlist(exp$DefaultParams)),stringsAsFactors=F)
for(coln in names(params.base)){
  if(!is.na(as.numeric(params.base[,coln]))){
    params.base[,coln] <- as.numeric(params.base[,coln])
  }
}

for(i in 1:nrow(factor.grid)){
  params <- params.base
  this.level <- factor.grid[i,] 
  for(factor.i in 1:length(this.level)){
    level.i <- which(unlist(lapply(exp$Factors[[factor.i]]$Levels,function(ll){ ll$Code == this.level[factor.i] })))
    for(param.val in exp$Factors[[factor.i]]$Levels[[level.i]]$Params){
      streval(pp('params$',names(param.val),'<- ',param.val))
    }
  }
  # Replace filenames with full path
  file.inds <- grep('File',names(params))
  params[,file.inds] <- pp(ev.amod.shared,'model/inputs/',params[,file.inds])

  my.cat('Experimental Group: ')
  print(this.level)
  my.cat('Parameter Values: ')
  print(params)

  # Now we can run the model
  #ev.amod.sim(params)
}

save(params,file=pp(ev.amod.shared,'model/inputs/devo-params.Rdata'))

