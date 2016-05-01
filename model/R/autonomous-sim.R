### TODO
source(pp(ev.amod.model,'model/R/run-experiment.R'))
source(pp(ev.amod.model,'model/R/misc-functions.R'))

animate.soln <- function(the.sys){
  my.cat('animating')
  the.sys[,var:=factor(var,c('u','sic','v','sid','w','simp','simn'))]
  for(the.node in na.omit(u(the.sys$node))){
    image.i <- 1
    for(the.t in u(the.sys$t)){
      png(filename = pp(exp$OutputsDirectory,'soln-node',the.node,'-img',sprintf('%05d',image.i),".png"),width=800,height=500)
      p<-ggplot(the.sys[!var%in%c('simp','simn') & node==the.node & t==the.t],aes(x=x,y=val,colour=var))+geom_line()+facet_grid(var~.,scales='free_y')+labs(title=pp('t=',the.t))+scale_color_manual(values=c('black','red','black','red','black'))+scale_y_continuous(limits=range(the.sys[!var%in%c('simp','simn')]$val))
      #p<-ggplot(the.sys[t==the.t],aes(x=x,y=val,colour=var))+geom_line()+facet_grid(var~.)+labs(title=pp('t=',the.t))+scale_y_continuous(limits=range(the.sys$val))
      print(p)
      dev.off()
      image.i <- image.i + 1
    }
    # system(pp('ffmpeg -framerate 30 -i ',exp$OutputsDirectory,'soln-node',the.node,'-img%05d.png -vf "scale=640:-1" -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -r 30 ',exp$OutputsDirectory,'soln-node',the.node,'.mp4'),intern=T,input='Y')
  }
  image.i <- 1
  for(the.t in u(the.sys$t)){
    png(filename = pp(exp$OutputsDirectory,'soln-transp-img',sprintf('%05d',image.i),".png"),width=500,height=500)
    p<-ggplot(the.sys[var%in%c('simp','simn') & dir=='to' & t==the.t],aes(x=x,y=val,colour=var))+geom_line()+facet_grid(inode~jnode,scales='free_y')+labs(title=pp('t=',the.t))+scale_color_manual(values=c('black','red','black','red','black'))+scale_y_continuous(limits=range(the.sys[var%in%c('simp','simn') & dir=='to']$val))
    print(p)
    dev.off()
    image.i <- image.i + 1
  }
  # system(pp('ffmpeg -framerate 30 -i ',exp$OutputsDirectory,'soln-transp-img%05d.png -vf "scale=640:-1" -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -r 30 ',exp$OutputsDirectory,'soln-transp.mp4'),intern=T,input='Y')
}

ev.amod.sim <- function(params,prev.solution=NULL,abs.t=0,prev.params=NULL,the.timeout=500,fleet.correct=F){
  
  my.cat('initializing')
  node.meta <- read.dt.csv(params$NodeFile)
  setkey(node.meta,id)
  odt <- read.matrix.csv(params$TravelTimeFile)
  ode <- read.matrix.csv(params$TravelEnergyFile)
  odf <- read.matrix.csv(params$TravelFareFile)
  dem <- read.dt.csv(params$TripDemand) # demand is specified in trips/min
  load <- read.dt.csv(params$PowerDemand) # load is specified in kW

  n.nodes <- nrow(node.meta)
  nodes <- node.meta$id
  xs <- seq(0,1,by=params$dx)
  n.x <- length(xs)
  ts <- seq(0,params$MovingHorizonT,by=params$dt)
  n.t <- length(ts)
  moving.horizon.dindex <- round(params$MovingHorizonDT/params$dt)
  prev.zero.index <- findInterval(params$MovingHorizonDT,ts)
  
  Qc <- function(x){ # SOE/min
    if(x==1)return(0)
    params$ChargingRate*params$ChargingEfficiency/params$BatteryCapacity/60 # 60 converts from hour to min
  }
  uLW1 <- function(x){
    Qc(x)*params$dt/params$dx/2
  }
  uLW2 <- function(x){
   (Qc(x)*params$dt/params$dx)^2/2
  }
  Qd <- function(x){ # SOE/min
    if(x==0)return(0)
    -params$DischargingRate/params$BatteryCapacity/60 
  }
  wLW1 <- function(x){
    Qd(x)*params$dt/params$dx/2
  }
  wLW2 <- function(x){
    (Qd(x)*params$dt/params$dx)^2/2
  }
  my.cat(pp('Courant #: ',roundC(Qc(0.5)*params$dt/params$dx,3),' (charging), ',roundC(-Qd(0.5)*params$dt/params$dx,3),' (discharging)'))

  #######################
  # OBJECTIVE
  #######################
  # First get a named vector for indexing the correct location in the decision/state vars
  state.combs <- expand.grid(xs,ts,nodes)
  names(state.combs) <- c('x','t','i')
  state.combs$row <- 1:nrow(state.combs)
  obj <- array(0,5*n.t*n.nodes*n.x)
  state.suffix <- ddply(state.combs,.(row),function(xx){ pp(pp(names(xx)[c(3,1,2)],xx[c(3,1,2)]),collapse='-')})$V1
  mob.combs <- expand.grid(xs,ts,nodes,nodes)
  names(mob.combs) <- c('x','t','i','j')
  mob.combs$row <- 1:nrow(mob.combs)

  obj <- array(0,5*n.t*n.nodes*n.x + 4*n.t*n.nodes^2*n.x)
  state.suffix <- ddply(state.combs,.(row),function(xx){ pp(pp(names(xx)[c(3,1,2)],xx[c(3,1,2)]),collapse='-')})$V1
  mob.suffix <- c(ddply(mob.combs,.(row),function(xx){ pp(xx[3],'to',xx[4],'-',pp(names(xx)[c(1,2)],xx[c(1,2)],collapse='-'))})$V1,ddply(mob.combs,.(row),function(xx){ pp(xx[3],'from',xx[4],'-',pp(names(xx)[c(1,2)],xx[c(1,2)],collapse='-'))})$V1)
  names(obj) <- c(pp('u-',state.suffix),pp('v-',state.suffix),pp('w-',state.suffix),pp('sic-',state.suffix),pp('sid-',state.suffix),pp('simp-',mob.suffix),pp('simn-',mob.suffix))

  if(!exists('prev.solution') || is.null(prev.solution)){
    boundary.correction.factor <- (1/params$dx)/(1/params$dx+1)
    prev.solution<-obj
    for(ix in 1:length(xs)){
      for(inode in 1:length(nodes)){
        prev.solution[pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[prev.zero.index])] <- params$FleetSize/n.nodes*boundary.correction.factor
      }
    }
    prev.params <- params
  }

  for(it in 1:length(ts)){
    for(ix in 1:length(xs)){
      for(inode in 1:length(nodes)){
        obj[pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- -params$CostCharge*params$dx*params$ChargingRate*params$dt/60
        obj[pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- node.meta[J(nodes[inode]),price.discharge]*params$dx*params$DischargingRate*params$dt/60
        for(jnode in 1:length(nodes)){
          obj[pp('simp-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- odf[inode,jnode]*params$dx*params$dt
        }
      }
    }
  }
  
  #######################
  # CONSTRAINTS
  #######################
  
  #------------------------------------ Equality constraints ------------------------------------------#
  
  #######################
  # Initial Conditions
  #######################

  # First we sum the total vehicles in the previous solution for the time step corresponding to time=0 in this run
  # Then we scale back to the fleet size to manage numerical leakage.
  prev.num.veh <- 0
  for(ix in 1:length(xs)){
    for(inode in 1:length(nodes)){
      prev.num.veh <- prev.num.veh + params$dx*(prev.solution[pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[prev.zero.index])] + 
                                                prev.solution[pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[prev.zero.index])] +
                                                prev.solution[pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[prev.zero.index])])
      for(jnode in 1:length(nodes)){
        trip.time.dindex <- round(odt[nodes[inode],nodes[jnode]]/params$dt)
        if(trip.time.dindex>0){
          for(it in 1:trip.time.dindex){
            prev.num.veh <- prev.num.veh + params$dx*params$dt*(
                              prev.solution[pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it + moving.horizon.dindex])] + 
                              prev.solution[pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it + moving.horizon.dindex])]) 
          }
        }
      }
    }
  }
  fleet.correction.fact <- ifelse(fleet.correct,mean(params$FleetSize / prev.num.veh,1),1)
  my.cat(pp('Fleet Correct by ',fleet.correction.fact))

  n.constr <- 3*n.nodes*n.x
  constr <- array(0,c(n.constr,length(obj)),dimnames=list(pp('init',1:n.constr),names(obj)))
  rhs <- array(0,n.constr)
  i.constr <- 1
  for(ix in 1:length(xs)){
    for(inode in 1:length(nodes)){
      for(state.var in c('u','v','w')){
        constr[i.constr,pp(state.var,'-i',nodes[inode],'-x',xs[ix],'-t0')] <- 1
        rhs[i.constr] <- prev.solution[pp(state.var,'-i',nodes[inode],'-x',xs[ix],'-t',ts[prev.zero.index])] * fleet.correction.fact
        i.constr <- i.constr + 1
      }
    }
  }
  #######################
  # Equations of State
  #######################
  if(params$Scheme=='upwind1'){ # Will use this one for coupled nodes 
    # First order upwind
    n.constr <- 3*(n.t-1)*n.nodes*n.x 
    constr <- rbind(constr,array(0,c(n.constr,length(obj)),dimnames=list(pp('state',1:n.constr),names(obj))))
    rhs <- c(rhs,array(0,n.constr))
    for(it in 1:(length(ts)-1)){
      for(ix in 1:length(xs)){
        for(inode in 1:length(nodes)){
          # u: charging state 
          if(ix>1)constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix-1],'-t',ts[it])] <- uLW1(xs[ix])*2
          constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1-uLW1(xs[ix])*2
          constr[i.constr,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- params$dt
          constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it+1])] <- -1
          i.constr <- i.constr + 1
          # v: idle state
          constr[i.constr,pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1
          constr[i.constr,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- -params$dt
          constr[i.constr,pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- -params$dt
          for(jnode in 1:length(nodes)){
            constr[i.constr,pp('simp-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- -params$dt
            constr[i.constr,pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- params$dt
            constr[i.constr,pp('simn-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- -params$dt
            constr[i.constr,pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- params$dt
          }
          constr[i.constr,pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[it+1])] <- -1
          i.constr <- i.constr + 1
          # w: discharging state
          constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1+wLW1(xs[ix])*2
          if(ix<length(xs))constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix+1],'-t',ts[it])] <- -wLW1(xs[ix])*2
          constr[i.constr,pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- params$dt
          constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it+1])] <- -1
          i.constr <- i.constr + 1
        }
      }
    }
  }else{
    my.cat(pp('Error: unrecognized discritization scheme: ',params$Scheme))
    return()
  }

  #######################
  # Boundary Conditions
  #######################
  n.constr <- 4*n.nodes*n.t
  constr <- rbind(constr,array(0,c(n.constr,length(obj)),dimnames=list(pp('bound',1:n.constr),names(obj))))
  rhs <- c(rhs,array(0,n.constr))
  for(it in 1:length(ts)){
    for(inode in 1:length(nodes)){
      # The flow rate from charging state to idle at full charge is equal to the total number of vehicles with full charge
      constr[i.constr,pp('u-i',nodes[inode],'-x1-t',ts[it])] <- 1/params$dt
      constr[i.constr,pp('sic-i',nodes[inode],'-x1-t',ts[it])] <- 1
      i.constr <- i.constr + 1
      # The flow rate from discharging state to idle at zero SOE is equal to the total number of vehicles with zero charged
      constr[i.constr,pp('w-i',nodes[inode],'-x0-t',ts[it])] <- 1/params$dt
      constr[i.constr,pp('sid-i',nodes[inode],'-x0-t',ts[it])] <- 1
      i.constr <- i.constr + 1
    }
  }
  #######################
  # Transport Conservation
  #######################
  n.constr <- 2*n.nodes^2*n.t*n.x # note that this is the max number of new constraints but it's complicated to calc 
                                  # the real number ahead of time so we figure out the final size dynamically
  new.constr <- array(0,c(n.constr,length(obj)),dimnames=list(pp('transport',1:n.constr),names(obj)))
  new.rhs <- array(0,n.constr)
  new.i.constr <- 1
  for(inode in 1:length(nodes)){
    for(jnode in 1:length(nodes)){
      trip.soe.dindex <- ceiling(ode[nodes[inode],nodes[jnode]]/params$BatteryCapacity/params$dx)
      trip.time.dindex <- round(odt[nodes[inode],nodes[jnode]]/params$dt)
      for(it in 1:(length(ts)-trip.time.dindex)){
        for(ix in (1+trip.soe.dindex):length(xs)){
          new.constr[new.i.constr,pp('simp-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
          new.constr[new.i.constr,pp('simp-',nodes[jnode],'from',nodes[inode],'-x',xs[ix-trip.soe.dindex],'-t',ts[it+trip.time.dindex])] <- -1 
          new.i.constr <- new.i.constr + 1
          new.constr[new.i.constr,pp('simn-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
          new.constr[new.i.constr,pp('simn-',nodes[jnode],'from',nodes[inode],'-x',xs[ix-trip.soe.dindex],'-t',ts[it+trip.time.dindex])] <- -1 
          new.i.constr <- new.i.constr + 1
        }
      }
    }
  }
  constr <- rbind(constr,new.constr[1:new.i.constr,])
  rhs <- c(rhs,new.rhs[1:new.i.constr])
  i.constr <- i.constr + new.i.constr  

  ##############################################
  # Arrivals from previous window
  # We need to schedule these arrivals that come
  # from decisions made in the recent past
  ##############################################
  n.constr <- 2*n.nodes^2*n.t*n.x # note that this is the max number of new constraints but it's complicated to calc 
                                  # the real number ahead of time so we figure out the final size dynamically
  new.constr <- array(0,c(n.constr,length(obj)),dimnames=list(pp('transport',1:n.constr),names(obj)))
  new.rhs <- array(0,n.constr)
  new.i.constr <- 1
  for(inode in 1:length(nodes)){
    for(jnode in 1:length(nodes)){
      trip.soe.dindex <- ceiling(ode[nodes[inode],nodes[jnode]]/params$BatteryCapacity/params$dx)
      trip.time.dindex <- round(odt[nodes[inode],nodes[jnode]]/params$dt)
      if(trip.time.dindex>0){
        for(it in 1:trip.time.dindex){
          for(ix in 1:(length(xs)-trip.soe.dindex)){
            new.constr[new.i.constr,pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
            new.rhs[new.i.constr] <- prev.solution[pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it+moving.horizon.dindex])] * fleet.correction.fact
            #new.rhs[new.i.constr] <- prev.solution[pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it+moving.horizon.dindex])]
            new.i.constr <- new.i.constr + 1
            new.constr[new.i.constr,pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
            new.rhs[new.i.constr] <- prev.solution[pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it+moving.horizon.dindex])] * fleet.correction.fact
            #new.rhs[new.i.constr] <- prev.solution[pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it+moving.horizon.dindex])]
            new.i.constr <- new.i.constr + 1
          }
        }
      }
      # set arrivals to zero if they are infeasible from an energy perspective
      if(trip.soe.dindex>0){
        for(it in 1:length(ts)){
          for(ix in (length(xs)-trip.soe.dindex+1):length(xs)){
            new.constr[new.i.constr,pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
            new.i.constr <- new.i.constr + 1
            new.constr[new.i.constr,pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
            new.i.constr <- new.i.constr + 1
          }
        }
      }
    }
  }
  constr <- rbind(constr,new.constr[1:new.i.constr,])
  rhs <- c(rhs,new.rhs[1:new.i.constr])
  i.constr <- i.constr + new.i.constr  

  ########################################################################
  ## Transport energy limit -- must have minium SOC to make a trip 
  ## the number of departures at SOC less than the minimum SOC to make a trip should be 0
  ## - load mobility data from i to j node at each time
  ## - rearrange mobility data in order of nodes and by ascending time
  ## - make the number of vehicles of less than the minimum SOC zero
  ########################################################################
  n.constr <- 2*n.nodes^2*n.t*n.x # note that this is the max number of new constraints but it's complicated to calc 
                                  # the real number ahead of time so we figure out the final size dynamically
  new.constr <- array(0,c(n.constr,length(obj)),dimnames=list(pp('trans.energy',1:n.constr),names(obj)))
  new.rhs <- array(0,n.constr)
  new.i.constr <- 1
  for(inode in 1:length(nodes)){ # from i
    for(jnode in 1:length(nodes)){ # to j
      for(it in 1:length(ts)){
        # Get required SOE from inode to jnode
        soe.min <-  ceiling(ode[nodes[inode],nodes[jnode]]/params$BatteryCapacity/params$dx)*params$dx
        # Set up the constraints -- flow of vehicles with states less than the minimum SOE to make a trip is zero
        for(ix in which(xs < soe.min)){
          new.constr[new.i.constr,pp('simp-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1 # with passenger
          new.i.constr <- new.i.constr + 1
          new.constr[new.i.constr,pp('simn-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1 # without passenger
          new.i.constr <- new.i.constr + 1
        }
      }
    }
  }
  constr <- rbind(constr,new.constr[1:new.i.constr,])
  rhs <- c(rhs,new.rhs[1:new.i.constr])
  i.constr <- i.constr + new.i.constr  
  
  #------------------------------------ Inequality Contraints ------------------------------------------#
  
  #######################
  # Flow Limits
  #######################
  n.constr <- 3*n.nodes*n.x*n.t
  constr.ineq <- array(0,c(n.constr,length(obj)),dimnames=list(pp('flow.lim',1:n.constr),names(obj)))
  rhs.ineq <- array(0,n.constr)
  i.constr.ineq <- 1
  for(it in 1:length(ts)){
    for(ix in 1:length(xs)){
      for(inode in 1:length(nodes)){
        # u >= sci: flow from charging to idle state should be less or equal to the total number of vehicles of given SOE in charging state 
        constr.ineq[i.constr.ineq,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- -1/params$dt
        constr.ineq[i.constr.ineq,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- -1
        i.constr.ineq <- i.constr.ineq + 1
        # w >= sdi: flow from discharging to idle state should be less or equal to the total number of vehicles of given SOE in discharging state 
        constr.ineq[i.constr.ineq,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- -1/params$dt
        constr.ineq[i.constr.ineq,pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- -1
        i.constr.ineq <- i.constr.ineq + 1
        # v >= sic + sid: flow from idle to charging/discharging states plus mobility should be less or equal to the total number of vehicles of given SOC
        constr.ineq[i.constr.ineq,pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- -1/params$dt
        constr.ineq[i.constr.ineq,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- 1
        constr.ineq[i.constr.ineq,pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- 1
        for(jnode in 1:length(nodes)){
          trip.time.dindex <- round(odt[nodes[inode],nodes[jnode]]/params$dt)
          constr.ineq[i.constr.ineq,pp('simp-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
          constr.ineq[i.constr.ineq,pp('simn-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
          constr.ineq[i.constr.ineq,pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- -1
          constr.ineq[i.constr.ineq,pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- -1
        }
        i.constr.ineq <- i.constr.ineq + 1
      }
    }
  }
  
  ########################################################################
  ## Discharging limit -- should be less or equal to the load demand
  ########################################################################
  n.constr <- 1*n.nodes*n.t
  constr.ineq <- rbind(constr.ineq, array(0,c(n.constr,length(obj)),dimnames=list(pp('discharge.lim',1:n.constr),names(obj))))
  rhs.ineq <- c(rhs.ineq,array(0,n.constr))
  for(inode in 1:length(nodes)){
    the.load <- load[node==nodes[inode]]
    for(it in 1:length(ts)){
      for(ix in 1:length(xs)){ 
        constr.ineq[i.constr.ineq,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- 1
      }
      #rhs.ineq[i.constr.ineq] <- min(the.load[findInterval(ts[it],the.load[,time]-abs.t)]$demand/params$DischargingRate/params$dx,
                                     #params$FleetSize/params$dx)*fleet.correction.fact*1.0000001
      if(it==1){
        rhs.ineq[i.constr.ineq] <- the.load[findInterval(ts[it],the.load[,time]-abs.t)]$demand/params$DischargingRate/params$dx*fleet.correction.fact
      }else{
        rhs.ineq[i.constr.ineq] <- the.load[findInterval(ts[it],the.load[,time]-abs.t)]$demand/params$DischargingRate/params$dx
      }
      i.constr.ineq <- i.constr.ineq + 1
    }
  }
  ########################################################################
  ## Transporting with passengers -- should be less or equal to the mobility demand *********************
  ########################################################################
  n.constr <- 1*n.nodes^2*n.t
  constr.ineq <- rbind(constr.ineq, array(0,c(n.constr,length(obj)),dimnames=list(pp('transp.lim',1:n.constr),names(obj))))
  rhs.ineq <- c(rhs.ineq,array(0,n.constr))
  for(inode in 1:length(nodes)){
    for(jnode in 1:length(nodes)){
      the.dem <- dem[orig==nodes[inode] & dest==nodes[jnode]]
      for(it in 1:length(ts)){
        for(ix in 1:length(xs)){ 
          constr.ineq[i.constr.ineq,pp('simp-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
        }
        rhs.ineq[i.constr.ineq] <- the.dem[findInterval(ts[it],the.dem[,time]-abs.t)]$demand/params$dx
        #rhs.ineq[i.constr.ineq] <- the.dem[findInterval(ts[it],the.dem[,time]-abs.t)]$demand/params$dx*fleet.correction.fact
        i.constr.ineq <- i.constr.ineq + 1
      }
    }
  }
  
  #------------------------------------------ End Contraints --------------------------------------------------#
  
  write.csv(obj,file=pp(exp$OutputsDirectory,'objective.csv'))
  write.csv(constr,file=pp(exp$OutputsDirectory,'constraints-eq.csv'))
  write.csv(constr.ineq,file=pp(exp$OutputsDirectory,'constraints-ineq.csv'))
  write.csv(rhs,file=pp(exp$OutputsDirectory,'rhs-eq.csv'))
  write.csv(rhs.ineq,file=pp(exp$OutputsDirectory,'rhs-ineq.csv'))
  
  #######################
  # Solve
  #######################
  lprec <- make.lp(0,ncol(constr),verbose='important') # verb can be normal detailed full
  lp.control(lprec,sense='max')
  lp.control(lprec,epslevel=params$EpsLevel)
  #lp.control(lprec,basis.crash='leastdegenerate')
  #lp.control(lprec,maxpivot=100)
  lp.control(lprec,pivoting=c('steepestedge','primalfallback'))
  lp.control(lprec,timeout=the.timeout)
  #lp.control(lprec,infinite=1e20)
  set.objfn(lprec, obj)
  for(i in 1:nrow(constr)){
    add.constraint(lprec, constr[i,], "=", rhs[i])
  }
  for(i in 1:nrow(constr.ineq)){
    add.constraint(lprec, constr.ineq[i,], "<=", rhs.ineq[i])
  }
  lower.b <- rep(0,length(obj))
  lower.b[grep('sic|sid',names(obj))] <- -Inf
  set.bounds(lprec, lower = lower.b, upper = rep(Inf,length(obj)))
  dimnames(lprec) <- list(c(rownames(constr), rownames(constr.ineq)),names(obj))
  set.basis(lprec,guess.basis(lprec,rep(0,length(obj))))

  my.cat('solving')
  status <- solve(lprec)
  sol <- list(status=status,solution=get.variables(lprec),obj=get.objective(lprec))
  my.cat(pp("result: ",status.code(status)))

  if(sol$status==0){
    names(sol$solution) <- names(obj)
    sol$sys <- parse.results(sol$solution)
  }else{
    write.csv(get.objective(lprec),file=pp(exp$OutputsDirectory,'debug-objective.csv'))
    write.csv(get.variables(lprec),file=pp(exp$OutputsDirectory,'debug-variables.csv'))
  }
  return(sol)
}

the.timeout <- 60
ev.amod.sim.horizon <- function(params,the.timeout=100,t.initial=0,final.solution=NULL,prev.solution=NULL,prev.params=NULL,fleet.correct=F){
  t.initials <- seq(t.initial,params$FullHorizonT,by=params$MovingHorizonDT)
  solution.list <- list()
  if(!is.null(final.solution))solution.list[[1]] <- final.solution
  t.initial <- t.initials[1]
  for(t.initial in t.initials){
    my.cat('##############################################################')
    my.cat(pp('Solving Window Starting at T=',t.initial))
    ptm <- proc.time()
    if(is.null(prev.solution)){
      sol <- tryCatch(ev.amod.sim(params,the.timeout=the.timeout,fleet.correct=fleet.correct),error=function(e){ print(e); return(NA) })
    }else{
      sol <- tryCatch(ev.amod.sim(params,prev.solution,the.timeout=the.timeout,abs.t=t.initial,prev.params=prev.params,fleet.correct=fleet.correct),error=function(e){ print(e); return(NA) })
    }
    if(length(sol)==1 | is.null(sol$sys)){
      stop('Error')
    }
    my.cat(pp('Solved in ',(proc.time() - ptm)[['elapsed']],' seconds'))
    prev.solution <- sol$solution
    solution.list[[length(solution.list)+1]] <- sol$sys[t<params$MovingHorizonDT]
    solution.list[[length(solution.list)]][,t:=t+t.initial]
    final.solution <- rbindlist(solution.list)
    # Write results as we go since you never know when it will hang
    write.csv(prev.solution,pp(exp$OutputsDirectory,'working-solution.csv'))
    write.csv(final.solution,pp(exp$OutputsDirectory,'final-solution.csv'),row.names=F)
  }
  final.solution
}

final.solution <- ev.amod.sim.horizon(params,fleet.correct=T)
#animate.soln(final.solution)

recovery.dir <- 'ModerateOutages-2016-05-01_11-49-41'
exp$OutputsDirectory <- pp(pp(head(str_split(exp$OutputsDirectory,"/")[[1]],-2),collapse="/"),"/",recovery.dir,"/")
load(file=pp(exp$OutputsDirectory,'params.Rdata'))
final.solution <- data.table(read.csv(pp(exp$OutputsDirectory,'final-solution.csv')))
working.solution <- read.csv(pp(exp$OutputsDirectory,'working-solution.csv'))
working.solution <- array(working.solution[,2],dimnames=list(working.solution[,1]))
#animate.soln(final.solution)

### Adjust params and reboot
prev.params <- params
#params$EpsLevel <- 'baggy'
params$MovingHorizonDT <- 20
params$MovingHorizonT <- 40
#params$dx <- 0.125
final.solution <- ev.amod.sim.horizon(params,t.initial=max(final.solution$t)+params$dt,
                                      final.solution=final.solution,prev.solution=working.solution,
                                      prev.params=prev.params,the.timeout=60,fleet.correct=T)
