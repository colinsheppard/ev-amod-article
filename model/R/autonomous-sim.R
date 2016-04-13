### TODO

source(pp(ev.amod.model,'model/R/run-experiment.R'))
load(pp(ev.amod.shared,'model/inputs/devo-params.Rdata'))
source(pp(ev.amod.model,'model/R/misc-functions.R'))
params$TravelDistanceFile = pp(ev.amod.shared,'model/inputs/travel-energy.csv');

animate.soln <- function(the.sys){
  my.cat('animating')
  the.sys[,var:=factor(var,c('u','sic','v','sid','w','simp','simn'))]
  for(the.node in u(the.sys$node)){
    image.i <- 1
    for(the.t in u(the.sys$t)){
      png(filename = pp(ev.amod.shared,'model/results/soln-node',the.node,'-img',sprintf('%05d',image.i),".png"),width=800,height=500)
      p<-ggplot(the.sys[!var%in%c('simp','simn') & node==the.node & t==the.t],aes(x=x,y=val,colour=var))+geom_line()+facet_grid(var~.,scales='free_y')+labs(title=pp('t=',the.t))+scale_color_manual(values=c('black','red','black','red','black'))+scale_y_continuous(limits=range(the.sys[!var%in%c('simp','simn')]$val))
      #p<-ggplot(the.sys[t==the.t],aes(x=x,y=val,colour=var))+geom_line()+facet_grid(var~.)+labs(title=pp('t=',the.t))+scale_y_continuous(limits=range(the.sys$val))
      print(p)
      dev.off()
      image.i <- image.i + 1
    }
    # system(pp('ffmpeg -framerate 30 -i ',ev.amod.shared,'model/results/soln-node',the.node,'-img%05d.png -vf "scale=640:-1" -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -r 30 ',ev.amod.shared,'model/results/soln-node',the.node,'.mp4'),intern=T,input='Y')
  }
  image.i <- 1
  for(the.t in u(the.sys$t)){
    png(filename = pp(ev.amod.shared,'model/results/soln-transp-img',sprintf('%05d',image.i),".png"),width=500,height=500)
    p<-ggplot(the.sys[var%in%c('simp','simn') & dir=='to' & t==the.t],aes(x=x,y=val,colour=var))+geom_line()+facet_grid(inode~jnode,scales='free_y')+labs(title=pp('t=',the.t))+scale_color_manual(values=c('black','red','black','red','black'))+scale_y_continuous(limits=range(the.sys$val))
    print(p)
    dev.off()
    image.i <- image.i + 1
  }
  # system(pp('ffmpeg -framerate 30 -i ',ev.amod.shared,'model/results/soln-transp-img%05d.png -vf "scale=640:-1" -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p -r 30 ',ev.amod.shared,'model/results/soln-transp.mp4'),intern=T,input='Y')
}

ev.amod.sim <- function(params,prev.solution=NULL){
  
  my.cat('initializing')
  node.meta <- read.dt.csv(params$NodeFile)
  odt <- read.matrix.csv(params$TravelTimeFile)
  ode <- read.matrix.csv(params$TravelEnergyFile)
  odf <- read.matrix.csv(params$TravelFareFile)
  dem <- read.dt.csv(params$TripDemand)
  load <- read.dt.csv(params$PowerDemand)

  n.nodes <- nrow(node.meta)
  nodes <- node.meta$id
  xs <- seq(0,1,by=params$dx)
  n.x <- length(xs)
  ts <- seq(0,params$T,by=params$dt)
  n.t <- length(ts)
  
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
  getTimeStepNum <- function(dt, inode, jnode){
    return(round(odt[inode,jnode]/dt))
  }
  getStateStepNum <- function(dx, inode, jnode){
    return(round(odd[inode,jnode]/params$BatteryCapacity/dx))
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

  if(!exists('prev.solution') || is.null(prev.solution))prev.solution<-obj

  for(it in 1:length(ts)){
    for(ix in 1:length(xs)){
      for(inode in 1:length(nodes)){
        obj[pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- -params$CostCharge*params$dx
        obj[pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- params$PriceDischarge*params$dx
        for(jnode in 1:length(nodes)){
          obj[pp('simp-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- params$PriceTrip*params$dx
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
  n.constr <- 3*n.nodes*n.x
  constr <- array(0,c(n.constr,length(obj)),dimnames=list(pp('init',1:n.constr),names(obj)))
  rhs <- array(0,n.constr)
  i.constr <- 1
  for(ix in 1:length(xs)){
    for(inode in 1:length(nodes)){
      constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t0')] <- 1
      rhs[i.constr] <- 0
      ## Constraints from Caroline's paper which I don't fully understand / agree with
      #rhs[i.constr] <- ifelse(ix==1,0,params$FleetSize/n.nodes/3)
      i.constr <- i.constr + 1
      constr[i.constr,pp('v-i',nodes[inode],'-x',xs[ix],'-t0')] <- 1
      rhs[i.constr] <- params$FleetSize/n.nodes
      i.constr <- i.constr + 1
      constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix],'-t0')] <- 1
      rhs[i.constr] <- 0
      ## Constraints from Caroline's paper which I don't fully understand / agree with
      #rhs[i.constr] <- ifelse(ix==length(xs),0,params$FleetSize/n.nodes/3)
      i.constr <- i.constr + 1
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
  }else if(params$Scheme=='upwind2'){
    # Second order upwind
    n.constr <- 3*(n.t-1)*n.nodes*n.x 
    constr <- rbind(constr,array(0,c(n.constr,length(obj)),dimnames=list(pp('state',1:n.constr),names(obj))))
    rhs <- c(rhs,array(0,n.constr))
    for(it in 1:(length(ts)-1)){
      for(ix in 1:length(xs)){
        for(inode in 1:length(nodes)){
          # u: charging state 
          constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1-3*uLW1(xs[ix])
          if(ix>1){
            constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix-1],'-t',ts[it])] <- 4*uLW1(xs[ix])
            if(ix>2){
              constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix-2],'-t',ts[it])] <- -uLW1(xs[ix])
            }else{
              # use first order upwind 
              constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1-uLW1(xs[ix])*2 # 2 accounts for not dividing by 2dx
              constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix-1],'-t',ts[it])] <- uLW1(xs[ix])*2
            }
          }else{
            # use first order upwind with no upwind party
            constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1-uLW1(xs[ix])*2 # 2 accounts for not dividing by 2dx
          }
          constr[i.constr,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- params$dt
          constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it+1])] <- -1
          i.constr <- i.constr + 1
          # v: idle state
          constr[i.constr,pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1
          constr[i.constr,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- -params$dt
          constr[i.constr,pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- -params$dt
          constr[i.constr,pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[it+1])] <- -1
          i.constr <- i.constr + 1
          # w: discharging state
          constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1+3*wLW1(xs[ix])
          if(ix<length(xs)){
            constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix+1],'-t',ts[it])] <- -4*wLW1(xs[ix])
            if(ix<(length(xs)-1)){
              constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix+2],'-t',ts[it])] <- wLW1(xs[ix])
            }else{
              # use first order downwind 
              constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- 1+wLW1(xs[ix])*2
              constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix+1],'-t',ts[it])] <- -wLW1(xs[ix])*2
            }
          }else{
            # use first order downwind with no downwind party
            constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- 1+wLW1(xs[ix])*2
          }
          constr[i.constr,pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- params$dt
          constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it+1])] <- -1
          i.constr <- i.constr + 1
        }
      }
    }
  }else if(params$Scheme=='lax-wendroff'){
    # lax wendroff
    n.constr <- 3*(n.t-1)*n.nodes*n.x 
    constr <- rbind(constr,array(0,c(n.constr,length(obj)),dimnames=list(pp('state',1:n.constr),names(obj))))
    rhs <- c(rhs,array(0,n.constr))
    for(it in 1:(length(ts)-1)){
      for(ix in 1:length(xs)){
        for(inode in 1:length(nodes)){
          # u: charging state 
          if(ix>1)constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix-1],'-t',ts[it])] <- uLW1(xs[ix])+uLW2(xs[ix])
          constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1-2*uLW2(xs[ix])
          if(ix<length(xs))constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix+1],'-t',ts[it])] <- uLW2(xs[ix])-uLW1(xs[ix])
          constr[i.constr,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- params$dt
          constr[i.constr,pp('u-i',nodes[inode],'-x',xs[ix],'-t',ts[it+1])] <- -1
          i.constr <- i.constr + 1
          # v: idle state
          constr[i.constr,pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1
          constr[i.constr,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- -params$dt
          constr[i.constr,pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- -params$dt
          constr[i.constr,pp('v-i',nodes[inode],'-x',xs[ix],'-t',ts[it+1])] <- -1
          i.constr <- i.constr + 1
          # w: discharging state
          if(ix>1)constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix-1],'-t',ts[it])] <- wLW1(xs[ix])+wLW2(xs[ix])
          constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])]   <- 1-2*wLW2(xs[ix])
          if(ix<length(xs))constr[i.constr,pp('w-i',nodes[inode],'-x',xs[ix+1],'-t',ts[it])] <- wLW2(xs[ix])-wLW1(xs[ix])
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
  new.i.constr <- 0
  for(inode in 1:length(nodes)){
    for(jnode in 1:length(nodes)){
      trip.soe.dindex <- round(ode[nodes[inode],nodes[jnode]]/params$BatteryCapacity/params$dx)
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
  new.i.constr <- 0
  for(inode in 1:length(nodes)){
    for(jnode in 1:length(nodes)){
      trip.soe.dindex <- round(ode[nodes[inode],nodes[jnode]]/params$BatteryCapacity/params$dx)
      trip.time.dindex <- round(odt[nodes[inode],nodes[jnode]]/params$dt)
      for(it in 1:trip.time.dindex){
        for(ix in 1:(length(xs)-trip.soe.dindex)){
          # new.constr[new.i.constr,pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
          # new.rhs[new.i.constr] <- prev.solution[pp('simp-',nodes[jnode],'from',nodes[inode],'-x',xs[trip.soe.dindex+ix],'-t',ts[length(ts)-trip.time.dindex+it])]
          # new.i.constr <- new.i.constr + 1
          # new.constr[new.i.constr,pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
          # new.rhs[new.i.constr] <- prev.solution[pp('simn-',nodes[jnode],'from',nodes[inode],'-x',xs[trip.soe.dindex+ix],'-t',ts[length(ts)-trip.time.dindex+it])]
          # new.i.constr <- new.i.constr + 1
        
          new.constr[new.i.constr,pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
          new.rhs[new.i.constr] <- prev.solution[pp('simp-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it+1])]
          new.i.constr <- new.i.constr + 1
          new.constr[new.i.constr,pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
          new.rhs[new.i.constr] <- prev.solution[pp('simn-',nodes[inode],'from',nodes[jnode],'-x',xs[ix],'-t',ts[it+1])]
          new.i.constr <- new.i.constr + 1
        }
      }
      # set arrivals to zero if they are infeasible from an energy perspective
      if(trip.time.dindex>0){
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
  new.i.constr <- 0
  for(inode in 1:length(nodes)){ # from i
    for(jnode in 1:length(nodes)){ # to j
      for(it in 1:length(ts)){
        # Get required SOE from inode to jnode
        soe.min <-  ode[nodes[inode],nodes[jnode]]/params$BatteryCapacity
        # Set up the constraints -- flow of vehicles with states less than the minimum SOE to make a trip is zero
        for(ix in which(xs <= soe.min)){
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

  #######################
  # FOR DEBUGGING SET SIGMAS TO ZERO
  #######################
  #n.constr <- 2*n.nodes*n.t*(n.x-1)
  ##n.constr <- 2*n.nodes*(n.x-1)*3
  #constr <- rbind(constr,array(0,c(n.constr,length(obj)),dimnames=list(pp('noflow',1:n.constr),names(obj))))
  #rhs <- c(rhs,array(0,n.constr))
  #for(it in 1:length(ts)){
  ##for(it in 1:3){
    #for(inode in 1:length(nodes)){
      #for(ix in 2:length(xs)){
        ##print(pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it]))
        #constr[i.constr,pp('sid-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- 1
        #i.constr <- i.constr + 1
      #}
      #for(ix in 1:(length(xs)-1)){
        ##print(pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it]))
        #constr[i.constr,pp('sic-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- 1
        #i.constr <- i.constr + 1
      #}
    #}
  #}
  
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
  constr.ineq <- rbind(constr.ineq, array(0,c(n.constr,length(obj)),dimnames=list(pp('dc.lim',1:n.constr),names(obj))))
  rhs.ineq <- c(rhs.ineq,array(0,n.constr))
  for(inode in 1:length(nodes)){
    the.load <- load[node==nodes[inode]]
    for(it in 1:length(ts)){
      for(ix in 1:length(xs)){ 
        constr.ineq[i.constr.ineq,pp('w-i',nodes[inode],'-x',xs[ix],'-t',ts[it])] <- 1
      }
      rhs.ineq[i.constr.ineq] <- the.load[findInterval(ts[it],the.load[,time])]$demand/params$DischargingRate/params$dx
      #my.cat(rhs.ineq[i.constr.ineq])
      i.constr.ineq <- i.constr.ineq + 1
    }
  }
  ########################################################################
  ## Transporting with passengers -- should be less or equal to the mobility demand *********************
  ########################################################################
  n.constr <- 1*n.nodes^2*n.t
  constr.ineq <- rbind(constr.ineq, array(0,c(n.constr,length(obj)),dimnames=list(pp('dc.lim',1:n.constr),names(obj))))
  rhs.ineq <- c(rhs.ineq,array(0,n.constr))
  for(inode in 1:length(nodes)){
    for(jnode in 1:length(nodes)){
      the.dem <- dem[orig==nodes[inode] & dest==nodes[jnode]]
      for(it in 1:length(ts)){
        for(ix in 1:length(xs)){ 
          constr.ineq[i.constr.ineq,pp('simp-',nodes[inode],'to',nodes[jnode],'-x',xs[ix],'-t',ts[it])] <- 1
        }
        rhs.ineq[i.constr.ineq] <- the.dem[findInterval(ts[it],the.dem[,time])]$demand/params$dx
        #my.cat(rhs.ineq[i.constr.ineq])
        i.constr.ineq <- i.constr.ineq + 1
      }
    }
  }
  
  #------------------------------------------ End Contraints --------------------------------------------------#
  
  #write.csv(rbind(obj,constr,constr.ineq),file=pp(ev.amod.shared,'model/debugging/tiny.csv'))
  #write.csv(c(rhs,rhs.ineq),file=pp(ev.amod.shared,'model/debugging/tiny-rhs.csv'))
  
  #######################
  # Solve
  #######################
  # lprec <- make.lp(nrow(constr),ncol(constr),verbose = "normal")
  lprec <- make.lp((nrow(constr)+nrow(constr.ineq)),ncol(constr),verbose='important') # row of inequality condition added
  lp.control(lprec,sense='max')
  lp.control(lprec,epslevel=params$EpsLevel)
  set.objfn(lprec, obj)
  for(i in 1:nrow(constr)){
    add.constraint(lprec, constr[i,], "=", rhs[i])
  }
  for(i in 1:nrow(constr.ineq)){
    add.constraint(lprec, constr.ineq[i,], "<=", rhs.ineq[i])
  }
  lower.b <- rep(0,length(obj))
  lower.b[grep('sic|sid',names(obj))] <- -Inf
  set.bounds(lprec, lower = lower.b)
  set.bounds(lprec, upper = rep(Inf,length(obj)))

  my.cat('solving')
  status <- solve(lprec)
  sol <- list(status=status,solution=get.variables(lprec),obj=get.objective(lprec))
  my.cat(pp("result: ",status.code(status)))

  if(sol$status==0){
    soln <- sol$solution
    names(soln) <- names(obj)
    sol$sys <- parse.results(soln)

    #ggplot(sys[t%%5==0],aes(x=x,y=val,colour=var))+geom_point()+facet_wrap(~t)
    #ggplot(sys[var%in%c('sic','sid')&t%%5==0],aes(x=x,y=val,colour=var))+geom_line()+facet_wrap(~t)
  }
  return(sol)
}

params$T <- 20
params$dx <- 0.02
params$dt <- 3
params$FleetSize <- 300
params$NodeFile <- pp(ev.amod.shared,'model/inputs/tiny/nodes-2.csv')
params$EpsLevel <- 'baggy' # 'tight','medium','loose','baggy'
params$Scheme <- 'upwind1' # 'lax-wendroff', 'upwind2'
ptm <- proc.time()
sol <- ev.amod.sim(params)
print(proc.time() - ptm)
sys <- sol$sys

# animate.soln(sol$sys)

# Debugging

# Check for conservation of vehicles
setkey(sys,t,var,x)
sys[-grep('si',var)][,list(count=sum((head(val,-1)+tail(val,-1))/2*params$dx)),by=c('t','var')][,list(count=sum(count)),by='t']

#sys[t<=5,list(round(sum(val*params$dx),1)),by=c('t','var')]
#sys[x>=0.9 & t<=5,list(round(sum(val*params$dx),1)),by=c('t','x','var')]

