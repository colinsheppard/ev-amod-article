load.libraries(c('maptools','diagram','TeachingDemos'))
# Recover from crash/freeze
#recovery.dir <- 'ExtremeOutages-2016-04-29_17-29-46' # success 5k veh
#recovery.dir <- 'ExtremeOutages-2016-04-29_19-53-10' # success 10k veh
#recovery.dir <- 'ExtremeOutages-2016-04-30_09-17-10' # success 20k veh
#recovery.dir <- 'ExtremeOutages-2016-04-30_09-37-44' # success 40k veh
#recovery.dir <- 'ModerateOutages-2016-05-01_11-49-41'
#recovery.dir <- 'ExtremeOutages-2016-05-01_12-25-04'

post.process <- function(recovery.dir){
  exp$OutputsDirectory <- pp(pp(head(str_split(exp$OutputsDirectory,"/")[[1]],-2),collapse="/"),"/",recovery.dir,"/")
  load(file=pp(exp$OutputsDirectory,'params.Rdata'))
  final.solution <- data.table(read.csv(pp(exp$OutputsDirectory,'final-solution.csv')))
  working.solution <- read.csv(pp(exp$OutputsDirectory,'working-solution.csv'))
  working.solution <- array(working.solution[,2],dimnames=list(working.solution[,1]))
  #animate.soln(final.solution)

  # Do results make sense?
  node.meta <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$NodeFile,"model")[[1]][2]))
  setkey(node.meta,id)
  odt <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelTimeFile,"model")[[1]][2]))
  ode <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelEnergyFile,"model")[[1]][2]))
  odf <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelFareFile,"model")[[1]][2]))
  dem <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$TripDemand,"model")[[1]][2])) # demand is specified in trips/min
  load <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$PowerDemand,"model")[[1]][2])) # load is specified in kW
  odf.m <- data.table(melt(odf))
  odf.m[,':='(inode=X1,jnode=X2,X1=NULL,X2=NULL)]
  ts <- seq(0,params$MovingHorizonT,by=params$dt)

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
  #ggplot(data.table(melt(max.revs[,list(rev.load=sum(rev.load)/1e3,rev.trip=sum(rev.trip)/1e3),by='t'],id.vars='t',measure.vars=c('rev.load','rev.trip'))),aes(x=t,y=value,colour=variable))+geom_line()+labs(x="Minute of Day",y="Revenue ($1000)",title=pp("Max Revenue Possible per 10 Minute Interval (",params$Title,")"))

  #ggplot(data.table(melt(max.revs[,list(dem.load.veh=sum(dem.load.veh),dem.trip=sum(dem.trip.per.min)*params$dt),by='t'],id.vars='t',measure.vars=c('dem.load.veh','dem.trip'))),aes(x=t,y=value,colour=variable))+geom_line()+labs(x="Minute of Day",y="# Vehicles",title=pp("Vehicles Required to Satisfy Demands (",params$Title,")"))

  #ggplot(data.table(melt(max.revs[,list(dem.load.veh=sum(rev.load,na.rm=T)/sum(dem.load.veh,na.rm=T),dem.trip=sum(rev.trip,na.rm=T)/params$dt/sum(dem.trip.per.min,na.rm=T)),by='t'],id.vars='t',measure.vars=c('dem.load.veh','dem.trip'))),aes(x=t,y=value,colour=variable))+geom_line()+labs(x="Minute of Day",y="Revenue / Vehicles / 10-min Interval ($/veh)",title=pp("Avg. Revenue Per Vehicle with Full Capture (",params$Title,")"))

  final.solution <- join.on(final.solution,node.meta,'node','id','price.discharge')
  final.solution <- join.on(final.solution,odf.m,c('inode','jnode'),c('inode','jnode'),'value')
  final.solution[var=='w',revenue:=val*price.discharge*params$dx*params$DischargingRate*params$dt/60]
  final.solution[var=='simp',revenue:=val*value*params$dx*params$dt]
  #val.sums <- final.solution[var%in%c('w','simp') & (is.na(dir) | dir=='to'),list(revenue=sum(revenue)),by=c('t','var')]
  #dev.new()
  #ggplot(val.sums,aes(x=t,y=revenue,colour=var))+geom_line()
  #val.sums <- final.solution[var%in%c('w','simp') & (is.na(dir) | dir=='to'),list(revenue=sum(revenue)),by=c('t','var','node')]
  #ggplot(val.sums,aes(x=t,y=revenue,colour=var))+geom_line()+facet_wrap(~node)

  # What does revenue look like as a fraction of max possible
  frac.rev <- rbindlist(list(join.on(final.solution[var=='w',list(revenue=sum(revenue)),by=c('t','node')],max.revs,c('t','node'),c('t','node'),'rev.load','max.'),join.on(final.solution[var=='simp'&dir=='to',list(revenue=sum(revenue)),by=c('t','inode')],max.revs,c('t','inode'),c('t','node'),'rev.trip','max.')),fill=T)
  frac.rev[is.na(node),':='(node=inode,max.rev=max.rev.trip,type='Trips')]
  frac.rev[is.na(inode),':='(max.rev=max.rev.load,type='Outtages')]
  frac.rev[,capture.fraction:=revenue/max.rev]
  write.csv(frac.rev,pp(exp$OutputsDirectory,'revenue.csv'))

  p <- ggplot(frac.rev,aes(x=t,y=capture.fraction,colour=type))+geom_line()+facet_wrap(~node)+labs(x='Minute',y='Revenue Captured / Max Possible Revenue',title=pp('Fraction of Revenue Captured by Node (',params$Title,': ',params$FleetSize,' veh)'))
  ggsave(p,file=pp(exp$OutputsDirectory,'revenue-capture-ratio-by-node.pdf'),width=8,height=4)

  p <- ggplot(frac.rev,aes(x=t,y=revenue/1e3,colour=type))+geom_line()+facet_wrap(~node)+labs(x="Minute of Day",y="Revenue ($1000)",title=pp('Revenue per 10 Minute Interval by Node (',params$Title,': ',params$FleetSize," veh)"))
  ggsave(p,file=pp(exp$OutputsDirectory,'revenue-by-node.pdf'),width=8,height=4)

  # What is the time varying fleet-wide SOE of the aggregated batteries
  final.solution[var%in%c('u','v','w'),num.veh:=val*params$dx]
  final.solution[!var%in%c('u','v','w'),num.veh:=val*params$dx*params$dt]
  final.solution <- join.on(final.solution,data.table(melt(round(odt/params$dt)))[,':='(dindex=value,value=NULL)],c('inode','jnode'),c('X1','X2'))
  final.solution[,next.time:=t+params$dt*(dindex-1)]
  final.solution <- join.on(final.solution,final.solution,c('var','dir','x','inode','jnode','t'),c('var','dir','x','inode','jnode','next.time'),'num.veh','next.')
  final.solution[is.na(next.num.veh) | dindex==1,next.num.veh:=0]
  final.solution[,batt.energy:=(num.veh+next.num.veh)*params$BatteryCapacity*x]
  tot.energy <- params$BatteryCapacity * params$FleetSize 
  soe <- final.solution[var%in%c('u','v','w','simp','simn') & (is.na(dir) | dir=='from'),list(soe=sum(batt.energy)/tot.energy),by='t']
  write.csv(soe,pp(exp$OutputsDirectory,'agg-soe.csv'))
  p <- ggplot(soe,aes(x=t,y=soe))+geom_line()+labs(x="Minute",y="SOE",title=pp('Aggregate SOE (',params$Title,': ',params$FleetSize," veh)"))
  ggsave(p,file=pp(exp$OutputsDirectory,'agg-soe.pdf'),width=8,height=4)
  num.veh.sys <- final.solution[var%in%c('u','v','w','simp','simn') & (is.na(dir) | dir=='from'),list(num=sum(num.veh+next.num.veh)),by='t']
  write.csv(num.veh.sys,pp(exp$OutputsDirectory,'num-veh-in-sys.csv'))
  p<-ggplot(num.veh.sys,aes(x=t,y=num))+geom_line()+labs(x="Minute",y="# Vehicles",title=pp('Number of Vehicles in System (',params$Title,': ',params$FleetSize," veh)"))
  ggsave(p,file=pp(exp$OutputsDirectory,'num-vehicles.pdf'),width=8,height=4)
}

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

# recovery.dir<-'Extreme-20k-Veh'
exp$OutputsDirectory <- pp(ev.amod.shared,'model/successful-runs/blah/')
vehicle.states <- function(recovery.dir){
  exp$OutputsDirectory <- pp(pp(head(str_split(exp$OutputsDirectory,"/")[[1]],-2),collapse="/"),"/",recovery.dir,"/")
  load(file=pp(exp$OutputsDirectory,'params.Rdata'))
  final.solution <- data.table(read.csv(pp(exp$OutputsDirectory,'final-solution.csv')))
  working.solution <- read.csv(pp(exp$OutputsDirectory,'working-solution.csv'))
  working.solution <- array(working.solution[,2],dimnames=list(working.solution[,1]))
  #animate.soln(final.solution)

  # Do results make sense?
  node.meta <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$NodeFile,"model")[[1]][2]))
  setkey(node.meta,id)
  odt <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelTimeFile,"model")[[1]][2]))
  ode <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelEnergyFile,"model")[[1]][2]))
  odf <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelFareFile,"model")[[1]][2]))
  dem <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$TripDemand,"model")[[1]][2])) # demand is specified in trips/min
  load <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$PowerDemand,"model")[[1]][2])) # load is specified in kW
  odf.m <- data.table(melt(odf))
  odf.m[,':='(inode=X1,jnode=X2,X1=NULL,X2=NULL)]
  ts <- seq(0,params$MovingHorizonT,by=params$dt)

  final.solution[var%in%c('u','v','w'),num.veh:=val*params$dx]
  final.solution[!var%in%c('u','v','w'),num.veh:=val*params$dx*params$dt]
  final.solution <- join.on(final.solution,data.table(melt(round(odt/params$dt)))[,':='(dindex=value,value=NULL)],c('inode','jnode'),c('X1','X2'))
  final.solution[,next.time:=t+params$dt*(dindex-1)]
  final.solution <- join.on(final.solution,final.solution,c('var','dir','x','inode','jnode','t'),c('var','dir','x','inode','jnode','next.time'),'num.veh','next.')
  final.solution[is.na(next.num.veh) | dindex==1,next.num.veh:=0]
  final.solution[,batt.energy:=(num.veh+next.num.veh)*params$BatteryCapacity*x]
  tot.energy <- params$BatteryCapacity * params$FleetSize 
  soe <- final.solution[var%in%c('u','v','w','simp','simn') & (is.na(dir) | dir=='from'),list(soe=sum(batt.energy)/tot.energy),by='t']

  num.veh.sys <- final.solution[var%in%c('u','v','w','simp','simn') & (is.na(dir) | dir=='to'),list(num=sum(num.veh)),by=c('t','var')]

  num.veh.sys.sum <- num.veh.sys[,list(num.sum=sum(num)),by='t']
  num.veh.sys.sum[,correction:=params$FleetSize/num.sum]
   
  num.veh.sys <- join.on(num.veh.sys,num.veh.sys.sum,'t','t')
  num.veh.sys[,num.corrected:=num*correction]
  num.veh.sys[,var:=revalue(var,c('simp'='Mobile with Passenger','simn'='Mobile no Passenger','u'='Charging','v'='Idle','w'='Discharging'))]
  p <- ggplot(num.veh.sys,aes(x=t,y=num.corrected,fill=var))+geom_area()+geom_line(position='stack')+labs(x="Minute",y="# Vehicles",fill="Vehicle State",colour=NULL,title=params$Title)+theme_bw()
  ggsave(p,file=pp(exp$OutputsDirectory,'states-area-plot-time-series.pdf'),width=8,height=4)

  return(num.veh.sys)
}

# recovery.dir<-'Extreme-20k-Veh'
exp$OutputsDirectory <- pp(ev.amod.shared,'model/successful-runs/blah/')
money.animation <- function(recovery.dir){
  exp$OutputsDirectory <- pp(pp(head(str_split(exp$OutputsDirectory,"/")[[1]],-2),collapse="/"),"/",recovery.dir,"/")
  load(file=pp(exp$OutputsDirectory,'params.Rdata'))
  final.solution <- data.table(read.csv(pp(exp$OutputsDirectory,'final-solution.csv')))
  working.solution <- read.csv(pp(exp$OutputsDirectory,'working-solution.csv'))
  working.solution <- array(working.solution[,2],dimnames=list(working.solution[,1]))
  #animate.soln(final.solution)

  # Do results make sense?
  node.meta <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$NodeFile,"model")[[1]][2]))
  setkey(node.meta,id)
  odt <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelTimeFile,"model")[[1]][2]))
  ode <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelEnergyFile,"model")[[1]][2]))
  odf <- read.matrix.csv(pp(ev.amod.shared,'model/',str_split(params$TravelFareFile,"model")[[1]][2]))
  dem <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$TripDemand,"model")[[1]][2])) # demand is specified in trips/min
  load <- read.dt.csv(pp(ev.amod.shared,'model/',str_split(params$PowerDemand,"model")[[1]][2])) # load is specified in kW
  odf.m <- data.table(melt(odf))
  odf.m[,':='(inode=X1,jnode=X2,X1=NULL,X2=NULL)]
  ts <- seq(0,params$MovingHorizonT,by=params$dt)

  final.solution[var%in%c('u','v','w'),num.veh:=val*params$dx]
  final.solution[!var%in%c('u','v','w'),num.veh:=val*params$dx*params$dt]
  final.solution <- join.on(final.solution,data.table(melt(round(odt/params$dt)))[,':='(dindex=value,value=NULL)],c('inode','jnode'),c('X1','X2'))
  final.solution[,next.time:=t+params$dt*(dindex-1)]
  final.solution <- join.on(final.solution,final.solution,c('var','dir','x','inode','jnode','t'),c('var','dir','x','inode','jnode','next.time'),'num.veh','next.')
  final.solution[is.na(next.num.veh) | dindex==1,next.num.veh:=0]
  final.solution[,batt.energy:=(num.veh+next.num.veh)*params$BatteryCapacity*x]
  tot.energy <- params$BatteryCapacity * params$FleetSize 
  soe <- final.solution[var%in%c('u','v','w','simp','simn') & (is.na(dir) | dir=='from'),list(soe=sum(batt.energy)/tot.energy),by='t']

  num.veh.sys <- final.solution[var%in%c('u','v','w','simp','simn') & (is.na(dir) | dir=='to'),list(num=sum(num.veh)),by=c('t','var')]

  num.veh.sys.sum <- num.veh.sys[,list(num.sum=sum(num)),by='t']
  num.veh.sys.sum[,correction:=params$FleetSize/num.sum]
   
  num.veh.sys <- join.on(num.veh.sys,num.veh.sys.sum,'t','t')
  num.veh.sys[,num.corrected:=num*correction]
  num.veh.sys[,var:=revalue(var,c('simp'='Mobile with Passenger','simn'='Mobile no Passenger','u'='Charging','v'='Idle','w'='Discharging'))]

  sf <- readShapePoly(pp(ev.amod.shared,'/data/sf_boundary/sf_boundary'))
  sf.bb <- sf@bbox
  node.pts <- data.table(node=1:4,
                      x=c(sf.bb[1,1]+diff(sf.bb[1,])/4-.025*diff(sf.bb[1,]),
                          sf.bb[1,1]+3*diff(sf.bb[1,])/4-.025*diff(sf.bb[1,])-7e3,
                          sf.bb[1,1]+diff(sf.bb[1,])/4-.025*diff(sf.bb[1,]),
                          sf.bb[1,1]+3*diff(sf.bb[1,])/4-.025*diff(sf.bb[1,])-3e3),
                      y=c(sf.bb[2,1]+3*diff(sf.bb[2,])/4-3e3,
                          sf.bb[2,1]+3*diff(sf.bb[2,])/4,
                          sf.bb[2,1]+diff(sf.bb[2,])/4,
                          sf.bb[2,1]+diff(sf.bb[2,])/4),
                      x2=c(sf.bb[1,1]+diff(sf.bb[1,])/8-.025*diff(sf.bb[1,]),
                          sf.bb[1,1]+6*diff(sf.bb[1,])/8-.025*diff(sf.bb[1,])+500,
                          sf.bb[1,1]+diff(sf.bb[1,])/4-.025*diff(sf.bb[1,]),
                          sf.bb[1,1]+6*diff(sf.bb[1,])/8-.025*diff(sf.bb[1,])+500),
                      y2=c(sf.bb[2,1]+5*diff(sf.bb[2,])/8-3e3,
                          sf.bb[2,1]+5*diff(sf.bb[2,])/8-3e3,
                          sf.bb[2,1]+diff(sf.bb[2,])/6,
                          sf.bb[2,1]+diff(sf.bb[2,])/12),key='node')
  setkey(final.solution,t)
  max.mob <- max(final.solution[dir=='to',sum(num.veh),by=c('t','var','inode','jnode')]$V1)
  max.charge <- max(final.solution[var=='u',sum(num.veh),by=c('t','node')]$V1)
  the.inset <- function(){
    plot(load[node!=3,list(dem=sum(demand)),by='time']$time, load[node!=3,list(dem=sum(demand)),by='time']$dem,type='l',xlab='',ylab='',main='',xaxt='n',yaxt='n')
    abline(v=the.t,col='red4')
  }

  make.dir(pp(exp$OutputsDirectory,'money'))
  image.i <- 1
  for(the.t in u(num.veh.sys$t)){
    png(pp(exp$OutputsDirectory,'money/money-img',sprintf('%05d',image.i),'.png'),width = 10, height = 10, units = "in",res=150)
    plot(sf)
    abline(v=sf.bb[1,1]+diff(sf.bb[1,])/2-.05*diff(sf.bb[1,]))
    abline(h=sf.bb[2,1]+diff(sf.bb[2,])/2)
    for(i.node in node.meta$id){
      for(j.node in node.meta$id){
        end.pt <- unlist(node.pts[J(j.node),list(x,y)])
        if(i.node==j.node){
          end.pt[1] <- end.pt[1] + ifelse(i.node==1,-8000,ifelse(i.node==2,9000,0))
          end.pt[2] <- end.pt[2] + ifelse(i.node==1,-4000,ifelse(i.node==2,6000,-10000))
          arr.pos <- 0.4
        }else{
          end.pt <- unlist(node.pts[J(j.node),list(x,y)])
          arr.pos <- 0.8
        }
        lwd.p <- 15*sum(final.solution[the.t==t & var=='simp' & inode==i.node & jnode==j.node & dir=='to']$num.veh)/max.mob
        lwd.n <- 15*sum(final.solution[the.t==t & var=='simn' & inode==i.node & jnode==j.node & dir=='to']$num.veh)/max.mob
        if(lwd.p>0.01){
          curvedarrow(unlist(node.pts[J(i.node),list(x,y)]), end.pt, lcol = 'blue', arr.col='blue', lwd = lwd.p, dr = 0.1, curve=0.25, 
                    endhead = T, segment = c(0, 1),arr.type = "triangle",arr.pos=arr.pos)
        }
        if(lwd.n>0.01){
          curvedarrow(unlist(node.pts[J(i.node),list(x,y)]), end.pt,lcol = 'grey', arr.col='grey',  lwd = lwd.n, dr = 0.01, curve=0.15, endhead = T, segment = c(0, 1),
                    arr.type = "triangle",arr.pos=arr.pos)
        }
      }
      ch.ht <- 8e3*sum(final.solution[the.t==t & var=='u' & node==i.node]$num.veh)/max.charge
      disch.ht <- 8e3*sum(final.solution[the.t==t & var=='w' & node==i.node]$num.veh)/max.charge
      segments(node.pts[J(i.node),x2]-1e3,node.pts[J(i.node),y2],y1=node.pts[J(i.node),y2]+ch.ht,lend=1,lwd=15,col='red3')
      segments(node.pts[J(i.node),x2]+1e3,node.pts[J(i.node),y2],y1=node.pts[J(i.node),y2]+disch.ht,lend=1,lwd=15,col='green4')
    }
    soe.ht.max <- 8e3
    soe.ht <- soe.ht.max*soe[t==the.t,soe]
    segments(node.pts[J(3),x2],node.pts[J(3),y2],y1=node.pts[J(3),y2]+soe.ht.max*1.003,lend=1,lwd=15,col='black')
    segments(node.pts[J(3),x2],node.pts[J(3),y2],y1=node.pts[J(3),y2]+soe.ht.max,lend=1,lwd=14.5,col='white')
    segments(node.pts[J(3),x2],node.pts[J(3),y2],y1=node.pts[J(3),y2]+soe.ht,lend=1,lwd=15,col='black')
    text(node.pts[J(3),x2],node.pts[J(3),y2],label=pp('SOE: ',roundC(soe[t==the.t,soe]*100,0),'%'),cex=1.25,pos=1)
    text(node.pts[1,x]-7e3,node.pts[1,y]+10e3,label=pp('Hour: ',roundC(the.t/60,1)),cex=1.25)
    title(main=params$Title,cex=2.5)
    legend('bottomleft',legend=c('With Passengers','Without Passengers','Charging','Discharging','Agg. Battery SOE'),col=c('blue','grey','red3','green4','black'),bg='white',lwd=4)
    subplot(the.inset(),x=c(node.pts[J(2),x],node.pts[J(2),x]+10e3)+10e3,y=c(node.pts[J(2),y],node.pts[J(2),y]+6e3)+6e3)
    dev.off()
    image.i <- image.i + 1
  }
  system(pp('ffmpeg -framerate 6 -i ',exp$OutputsDirectory,'money/money-img%05d.png -vf "scale=640:-1" -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -r 30 ',exp$OutputsDirectory,'money-animation.mp4'),intern=T,input='Y')
}
