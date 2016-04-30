
# Recover from crash/freeze
recovery.dir <- 'ExtremeOutages-2016-04-29_17-29-46'
exp$OutputsDirectory <- pp(pp(head(str_split(exp$OutputsDirectory,"/")[[1]],-2),collapse="/"),"/",recovery.dir,"/")
load(file=pp(exp$OutputsDirectory,'params.Rdata'))
final.solution <- data.table(read.csv(pp(exp$OutputsDirectory,'final-solution.csv')))
working.solution <- read.csv(pp(exp$OutputsDirectory,'working-solution.csv'))
working.solution <- array(working.solution[,2],dimnames=list(working.solution[,1]))
#animate.soln(final.solution)

# Do results make sense?
node.meta <- read.dt.csv(params$NodeFile)
setkey(node.meta,id)
odt <- read.matrix.csv(params$TravelTimeFile)
ode <- read.matrix.csv(params$TravelEnergyFile)
odf <- read.matrix.csv(params$TravelFareFile)
dem <- read.dt.csv(params$TripDemand) # demand is specified in trips/min
load <- read.dt.csv(params$PowerDemand) # load is specified in kW
odf.m <- data.table(melt(odf))
odf.m[,':='(inode=X1,jnode=X2,X1=NULL,X2=NULL)]

# What is the max revenue we would expect if all load and demand is served
nodes <- node.meta$id
max.load.revs <- data.table(t=rep(seq(0,params$FullHorizonT,by=params$dt),length(nodes)),node=rep(nodes,each=length(seq(0,params$FullHorizonT,by=params$dt))))
max.load.revs <- join.on(max.load.revs,node.meta,'node','id','price.discharge')
setkey(load,node)
setkey(max.load.revs,node)
for(the.node in nodes){
  max.load.revs[J(the.node),dem.load:=load[J(the.node)]$demand[findInterval(t,load[J(the.node)]$time)]]
}
max.load.revs[,rev.load:=dem.load*params$dt/60*price.discharge]
max.load.revs[,dem.load.veh:=dem.load/params$DischargingRate]
dem.with.price <- join.on(dem,odf.m,c('orig','dest'),c('inode','jnode'),'value')
dem.with.price[,':='(dem.trips=demand,demand=NULL,price.trips=value,value=NULL)]
setkey(dem.with.price,orig,dest)
max.trip.revs <- data.table(expand.grid(list(t=seq(0,params$FullHorizonT,by=params$dt),inode=nodes,jnode=nodes)),key=c('inode','jnode'))
for(i.node in nodes){
  for(j.node in nodes){
    max.trip.revs[J(i.node,j.node),dem.trips:=dem.with.price[J(i.node,j.node)]$dem.trips[findInterval(t,dem.with.price[J(i.node,j.node)]$time)]]
    max.trip.revs[J(i.node,j.node),price.trips:=dem.with.price[J(i.node,j.node)]$price.trips[findInterval(t,dem.with.price[J(i.node,j.node)]$time)]]
  }
}
max.trip.revs[,rev.trips:=dem.trips*params$dt*price.trips]
max.revs <- join.on(max.load.revs,max.trip.revs[,list(rev.trip=sum(rev.trips),dem.trip.per.min=sum(dem.trips)),by=c('inode','t')],c('t','node'),c('t','inode'),c('rev.trip','dem.trip.per.min'))
#ggplot(data.table(melt(max.revs[,list(rev.load=sum(rev.load)/1e3,rev.trip=sum(rev.trip)/1e3),by='t'],id.vars='t',measure.vars=c('rev.load','rev.trip'))),aes(x=t,y=value,colour=variable))+geom_line()+labs(x="Minute of Day",y="Revenue ($1000)",title="Max Revenue Possible per 10 Minute Interval (Extreme Scenario)")

#ggplot(data.table(melt(max.revs[,list(dem.load.veh=sum(dem.load.veh),dem.trip=sum(dem.trip.per.min)*params$dt),by='t'],id.vars='t',measure.vars=c('dem.load.veh','dem.trip'))),aes(x=t,y=value,colour=variable))+geom_line()+labs(x="Minute of Day",y="# Vehicles",title="Vehicles Required to Satisfy Demands (Extreme Scenario)")

ggplot(data.table(melt(max.revs[,list(dem.load.veh=sum(rev.load,na.rm=T)/sum(dem.load.veh,na.rm=T),dem.trip=sum(rev.trip,na.rm=T)/params$dt/sum(dem.trip.per.min,na.rm=T)),by='t'],id.vars='t',measure.vars=c('dem.load.veh','dem.trip'))),aes(x=t,y=value,colour=variable))+geom_line()+labs(x="Minute of Day",y="Revenue / Vehicles / 10-min Interval ($/veh)",title="Avg. Revenue Per Vehicle with Full Capture (Extreme Scenario)")

final.solution <- join.on(final.solution,node.meta,'node','id','price.discharge')
final.solution <- join.on(final.solution,odf.m,c('inode','jnode'),c('inode','jnode'),'value')
final.solution[var=='w',revenue:=val*price.discharge*params$dx*params$DischargingRate*params$dt/60]
final.solution[var=='simp',revenue:=val*value*params$dx*params$dt]
val.sums <- final.solution[var%in%c('w','simp') & (is.na(dir) | dir=='to'),list(revenue=sum(revenue)),by=c('t','var')]
#dev.new()
ggplot(val.sums,aes(x=t,y=revenue,colour=var))+geom_line()
val.sums <- final.solution[var%in%c('w','simp') & (is.na(dir) | dir=='to'),list(revenue=sum(revenue)),by=c('t','var','node')]
ggplot(val.sums,aes(x=t,y=revenue,colour=var))+geom_line()+facet_wrap(~node)

# What does revenue look like as a fraction of max possible
frac.rev <- rbindlist(list(join.on(final.solution[var=='w',list(revenue=sum(revenue)),by=c('t','node')],max.revs,c('t','node'),c('t','node'),'rev.load','max.'),join.on(final.solution[var=='simp'&dir=='to',list(revenue=sum(revenue)),by=c('t','inode')],max.revs,c('t','inode'),c('t','node'),'rev.trip','max.')),fill=T)
frac.rev[is.na(node),':='(node=inode,max.rev=max.rev.trip,type='Outage Load')]
frac.rev[is.na(inode),':='(max.rev=max.rev.load,type='Trip Demand')]
frac.rev[,capture.fraction:=revenue/max.rev]
ggplot(frac.rev,aes(x=t,y=capture.fraction,colour=type))+geom_line()+facet_wrap(~node)+labs(x='Minute',y='Revenue Captured / Max Possible Revenue',title=pp('Fraction of Revenue Captured by Node (Extreme Scenario: ',params$FleetSize,' veh)'))


## Debugging
#parse.results(working.solution)->working.sys
#ggplot(working.sys[t==0],aes(x=x,y=val,colour=dir))+geom_point()+facet_wrap(node~var)
#dev.new()
#ggplot(load,aes(x=time,y=demand))+geom_line()+facet_wrap(~node)
#dev.new()
#ggplot(dem,aes(x=time,y=demand))+geom_line()+facet_grid(orig~dest)

## Check for conservation of vehicles
#setkey(final.solution,t,var,x)
#final.solution[-grep('sic|sid',var)][is.na(dir) | dir=='to',list(count=sum((head(val,-1)+tail(val,-1))/2*params$dx)),by=c('t','var')]
#state.sums <- final.solution[-grep('si',var)][,list(count=sum((head(val,-1)+tail(val,-1))/2*params$dx)),by=c('t','var')][,list(count=sum(count)),by='t']
#net.sums <- join.on(state.sums,final.solution[grep('simp|simn',var)][,list(count=sum(val)*params$dx*params$dt),by=c('t','dir')],'t','t',c('count','dir'),'trans.')
#net.sums <- data.table(as.data.frame(cast(melt(net.sums,id.vars=c('t','trans.dir')),t ~ variable + trans.dir)))
#net.sums[,':='(count=count_from,count_from=NULL,count_to=NULL)]
#net.sums[,count.next:=count+trans.count_from-trans.count_to]
#final.solution[grep('simp|simn',var)][is.na(dir),list(count=sum(val)*params$dx*params$dt),by=c('t','var')]
#final.solution[t<=5,list(round(sum(val*params$dx),1)),by=c('t','var')]
#final.solution[x>=0.9 & t<=5,list(round(sum(val*params$dx),1)),by=c('t','x','var')]

