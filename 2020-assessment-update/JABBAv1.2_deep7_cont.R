##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## Stock Assessment execution File for JABBA
## Developed by Henning Winker & Felipe Carvalho (Cape Town/Hawaii)  
##><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><

cat(paste0("\n","><>><>><>><>><>><>><>><>"))
cat(paste0("\n","><> Run Model ",Mod.names,"<><"))
cat(paste0("\n","><>><>><>><>><>><>><>><>","\n","\n"))
# setwd(paste(File))
dir.create(paste0(File,"/",assessment),showWarnings = FALSE)
dir.create(paste0(File,"/",assessment,"/",Scenario,"_",Mod.names),showWarnings = FALSE)
dir.create(paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Input"),showWarnings = FALSE)
input.dir = paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Input")

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Define objects to make sure they exist if not included in Prime file
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
if(exists("igamma")==FALSE) igamma = c(4,0.01)  # Generic process error prior
if(exists("BmsyK")==FALSE) BmsyK = 0.4  # JABBA default for Pella model
if(exists("Model")==FALSE){ model = 1; Mod.names = c("Schaefer")} # Run Schaefer if model is not specified
if(exists("proc.dev.all")==FALSE) proc.dev.all = FALSE # process error deviation if only catch is available  
if(exists("Plim")==FALSE) Plim = 0  # Standard non-compound model
if(exists("P_bound")==FALSE) P_bound = c(0.02,1)  # Soft penalty bounds for P 
if(exists("KOBE.plot")==FALSE) KOBE.plot = TRUE # Produces JABBA Kobe plot 
if(exists("KOBE.type")==FALSE) KOBE.type = c("ICCAT","IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange) 
if(exists("SP.plot")==FALSE) SP.plot = c("standard","phase")[2] # Produces standard or 'Kobe phase' SP plot  
if(exists("Biplot")==FALSE) Biplot= TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
if(exists("save.trajectories")==FALSE) save.trajectories =FALSE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object 
if(exists("harvest.label")==FALSE) harvest.label = c("Hmsy","Fmsy")[2] # choose label preference H/Hmsy versus Fmsy
if(exists("CPUE.plot")==FALSE) CPUE.plot= TRUE # Runs state-tool to produce "alligned" multi-CPUE plot  
if(exists("catch.metric")==FALSE) catch.metric = "(t)" # Runs state-tool to produce "alligned" multi-CPUE plot  
if(exists("meanCPUE")==FALSE) meanCPUE = FALSE # Uses averaged CPUE from state-space tool instead of individual indices  
if(exists("Projection")==FALSE) Projection = FALSE # Use Projections: requires to define TACs vectors 
if(exists("save.projections")==FALSE) save.projections = FALSE# saves projection posteriors as .RData object 
if(exists("Reproduce.seed")==FALSE) Reproduce.seed = FALSE # If FALSE a random seed assigned to each run (default)
if(exists("TACint")==FALSE) TACint = mean(catch[nrow(catch)-3,2]:catch[nrow(catch),2]) # use mean catch from last years
if(exists("imp.yr")==FALSE) imp.yr = as.numeric(format(Sys.Date(), "%Y"))+1 # use next year from now
if(exists("init.values")==FALSE) init.values =FALSE # Allows to add manual starting values for K, r, q
if(exists("sigmaobs_bound")==FALSE) sigmaobs_bound = 1 # Adds an upper bound to the observation variance  
if(exists("sigmaproc_bound")==FALSE) sigmaproc_bound = 0.2 # Adds an upper bound to the process variance  
if(exists("q_bounds")==FALSE) q_bounds= c(10^-30,100) # Defines lower and upper bounds for q 
if(exists("K_bounds")==FALSE) K_bounds= c(10,10^8) # Defines lower and upper bounds for q 
# Save entire posterior as .RData object
if(exists("save.all")==FALSE) save.all = FALSE #  
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

## Prepare input data ----

cat(paste0("\n","><> Prepare input data <><","\n"))
indices = names(cpue)[2:ncol(cpue)]
n.indices = max(length(indices),1)
catches = names(catch)[2:ncol(catch)]
n.catches = length(catches)

years=catch[,1]
styr = min(years)
endyr = max(years)
n.years = length(years)
styr.cpue = min(cpue[,1])
styr.I = styr.cpue-styr+1 
syrs = which(!is.na(cpue$CPUE_3)) #added by Meg years of survey data


# Convert input data to matrices for JAGS input
conv.cpue = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
CPUE=matrix(conv.cpue,nrow=n.years,ncol=n.indices)


if(SE.I==FALSE){
  se = cpue  
  conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(cpue[,-1])))
  se2 = matrix(ifelse(fixed.obsE>0,fixed.obsE^2,10^-10),n.years,n.indices)#/2
} else{
  conv.se = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.indices),styr.I-1,n.indices),as.matrix(se[,-1])))
  #conv.se = sqrt(conv.se^2+fixed.obsE^2) 
  se2 = matrix(ifelse(is.na(conv.se),0.3^2,conv.se)^2,n.years,n.indices)+fixed.obsE^2#/2
}

conv.catch = as.numeric(rbind(matrix(rep(NA,(styr.I-1)*n.catches),styr.I-1,n.catches),as.matrix(catch[,-1])))
Catch=matrix(conv.catch,nrow=n.years,ncol=n.catches)
Catch[is.na(Catch)] = 0 # Replace any NA by zero 

# Total Catch
TC = apply(Catch,1,sum)

# Plot Catch 

cat(paste0("\n","><> Plot Catch in Input subfolder <><","\n"))

Par = list(mfrow=c(1,1),mar = c(5, 5, 1, 1), mgp =c(3,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/Catches_",assessment,".png"), width = 7, height = 5, 
    res = 200, units = "in")
par(Par)
plot(catch[,1],catch[,1],ylim=c(0,max(catch[,2:ncol(catch)],na.rm=TRUE)),ylab=paste0("Catch ",catch.metric),xlab="Year",type="n")
for(i in 2:ncol(catch)) lines(catch[,1],catch[,i],lty=(i-1),lwd=2)
legend("topright",paste(names(catch)[2:ncol(catch)]),lty=1:(ncol(catch)-1),bty="n")
dev.off()

#---------------------
# Index color palette
#---------------------
jabba.colors = as.character(rep(c('#e6194b', "#3cb44b", "#ffe119",
                                  "#0082c8","#f58231", "#911eb4",
                                  "#46f0f0", "#f032e6", "#d2f53c",
                                  "#fabebe", "#008080","#e6beff", "#aa6e28"),2))
#--------------------
# Set seed
#--------------------
if(Reproduce.seed==FALSE){ set.seed(ceiling(runif(1,min=0,max=1e6))) } else {set.seed(123)}

#---------------------------------------------------------------------------
# CPUE run State-Space model for averaging CPUE
#---------------------------------------------------------------------------
if(CPUE.plot==TRUE){ 
  cat(paste0("\n","><> Run State-Space CPUE averaging tool","\n"))
  #find first time-series with first CPUE
  q1.y = c(1:n.years)[is.na(apply(CPUE,1,mean,na.rm=TRUE))==FALSE][1] #first year with CPUE
  q1.I = which.max(CPUE[q1.y,])
  
  qs = c(q1.I,c(1:(ncol(cpue)-1))[-q1.I])
  
  
  sink("cpueAVG.jags")
  cat("
      model {
      
      # Prior specifications  
      eps <- 0.0000000000001 # small constant    
      
      iq[1] ~ dgamma(1000,1000)
      q[1] <-  pow(iq[1],-1)
      logq[1] <- log(1)
      for(i in 2:nI){
      iq[i] ~ dgamma(0.001,0.001)
      q[i] <- pow(iq[i],-1)
      logq[i] <-  log(q[i])
      }
      
      
      ")
  
  if(sigma.proc==TRUE){
    cat("
        # Process variance
        isigma2 <- isigma2.est 
        sigma2 <- pow(isigma2,-1)
        sigma <- sqrt(sigma2)
        fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
        ",append=TRUE)  
  }else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
             sigma2 <- pow(isigma2,-1)
             sigma <- sqrt(sigma2)
             
             ",append=TRUE)}
  
  if(sigma.est==TRUE){
    cat("
        # Obsevation variance
        # Observation error
        itau2~ dgamma(0.001,0.001)
        tau2 <- 1/itau2
        
        
        for(i in 1:nI)
        {
        for(t in 1:N)
        {
        var.obs[t,i] <- SE2[t,i]+tau2
        ivar.obs[t,i] <- 1/var.obs[t,i]
        # note total observation error (TOE)     
        TOE[t,i] <- sqrt(var.obs[t,i])
        
        }}
        ",append=TRUE)  
  }else{ cat(" 
      # Obsevation variance
             # Observation error
             itau2~ dgamma(2,2)
             tau2 <- 1/itau2
             
             
             for(i in 1:nI)
             {
             for(t in 1:N)
             {
             var.obs[t,i] <- SE2[t,i] # drop tau2
             fake.tau[t,i] <- tau2
             
             ivar.obs[t,i] <- 1/var.obs[t,i]
             # note total observation error (TOE)     
             TOE[t,i] <- sqrt(var.obs[t,i])
             
             }}
             
             ",append=TRUE)}
  
  # Run rest of code  
  cat("  
      # Process variance prior
      isigma2.est ~ dgamma(0.2,0.1)
      
      
      # Priors and constraints
      logY.est[1] ~ dnorm(logY1, 1)       # Prior for initial population size
      
      mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
      
      # Likelihood
      # State process
      for (t in 1:(N-1)){
      r[t] ~ dnorm(mean.r, isigma2)
      logY.est[t+1] <- logY.est[t] + r[t] }
      
      # Observation process
      for (t in 1:N) {
      for(i in 1:nI){
      y[t,i] ~ dnorm(logY.est[t]+logq[i], ivar.obs[t,i])
      }}
      
      # Population sizes on real scale
      for (t in 1:N) {
      Y.est[t] <- exp(logY.est[t])
      }
      
  } 
      ",fill = TRUE)
  sink()
  
  
  
  q.init = 1
  mCPUE = as.matrix(CPUE[q1.y:n.years,qs])
  mSE2 = as.matrix(se2[q1.y:n.years,qs])
  if(n.indices>1) for(i in 2:n.indices){q.init[i] = mean(mCPUE[,i],na.rm=TRUE)/mean(mCPUE[,1],na.rm=TRUE)}
  # Bundle data
  jags.data <- list(
    y = log(mCPUE),
    SE2=mSE2, 
    logY1 = log(mCPUE[1,1]),
    N = length(q1.y:n.years),
    nI=n.indices,
    sigma.fixed=ifelse(sigma.proc==TRUE,0,sigma.proc)
    )
  
  # Initial values
  inits <- function(){list(isigma2.est=runif(1,20,100), itau2=runif(1,80,200), mean.r = rnorm(1),iq = 1/q.init)}
  
  # Parameters monitored
  parameters <- c("mean.r", "sigma","r", "Y.est","q")
  
  
  # Call JAGS from R (BRT 3 min)
  mod.cpue <- jags(jags.data, inits, parameters, "cpueAVG.jags", n.chains = nc, n.thin = max(nt,2), n.iter = max(ni/5,10000), n.burnin = nb/10)
  
  
  cat(paste0("\n","><> Plot State-Space CPUE fits in Input subfolder <><","\n"))
  # get individual trends
  fitted <- lower <- upper <- NULL
  cpue.yrs = years[q1.y:n.years]
  
  for (t in 1:nrow(mCPUE)){
    fitted[t] <- median(mod.cpue$BUGSoutput$sims.list$Y.est[,t])
    lower[t] <- quantile(mod.cpue$BUGSoutput$sims.list$Y.est[,t], 0.025)
    upper[t] <- quantile(mod.cpue$BUGSoutput$sims.list$Y.est[,t], 0.975)}
  
  
  q.adj = apply(mod.cpue$BUGSoutput$sims.list$q,2,median)
  
  
  Par = list(mfrow=c(1,1),mar = c(3.5, 3.5, 0.1, 0.1), mgp =c(2.,0.5,0), tck = -0.02,cex=0.8)
  png(file = paste0(input.dir,"/CPUE_",assessment,"_",Scenario,".png"), width = 5, height = 3.5, 
      res = 200, units = "in")
  par(Par)
  u.ylim = NULL
  for(i in 1:n.indices){ u.ylim = c(u.ylim,exp(log(mCPUE[,i]/q.adj[i])+1.96*sqrt(mSE2[,i])))}  
  ylim = c(0,max(u.ylim,na.rm=TRUE))
  plot(0, 0, ylim = ylim, xlim = range(cpue.yrs), ylab = "Expected CPUE", xlab = "Year", col = "black", type = "n")
  legend("topright",paste(indices),lwd=2,col=(jabba.colors)[1:n.indices],bty="n")
  polygon(x = c(cpue.yrs,rev(cpue.yrs)), y = c(lower,rev(upper)), col = "gray", border = "gray90")
  
  for(i in 1:n.indices)
  {
    shift = runif(1,-0.1,0.1)
    cols=jabba.colors[qs[i]]
    plotCI(cpue.yrs+shift,mCPUE[,i]/q.adj[i],ui=exp(log(mCPUE[,i]/q.adj[i])+1.96*sqrt(mSE2[,i])),li=exp(log(mCPUE[,i]/q.adj[i])-1.96*sqrt(mSE2[,i])),add=TRUE,col= cols,pt.bg = cols,pch=21,gap=0)
    lines(cpue.yrs+shift,mCPUE[,i]/q.adj[i], col = cols,lwd=2)
    points(cpue.yrs+shift,mCPUE[,i]/q.adj[i], bg = cols,pch=21)
  }
  lines(cpue.yrs,fitted,lwd=2)
  
  dev.off()
  
  logSE = apply(log(mod.cpue$BUGSoutput$sims.list$Y.est),2,sd)
  
  
  if(nrow(mCPUE)<n.years) {
    fitted = c(rep(NA,q1.y-1),fitted)
    logSE = c(rep(0.2,q1.y-1),logSE)
  }    
  avgCPUE = data.frame(Year=years,CPUE= fitted,logSE=logSE)
  
  write.csv(avgCPUE,paste0(input.dir,"/avgCPUE_",assessment,"_",Scenario,".csv"))
  
  if(meanCPUE==TRUE){
    cat(paste0("\n","><> Use average CPUE as input for JABBA <><","\n"))
    
    CPUE = as.matrix(avgCPUE[,2]) 
    cpue.check = cpue[,-1]
    cpue.check[is.na(cpue[,-1])]=0
    CPUE[,1] = ifelse(apply(cpue.check,1,sum)==0,rep(NA,length(CPUE[,1])),CPUE[,1])
    se2 =  as.matrix(avgCPUE[,3]^2)     
    n.indices=1
    indices = "All"
    sets.q =1
    sets.var =1
  }
  
  }

#--------------------------------------------------------------------------
# END of CPUE State-Space tool
#--------------------------------------------------------------------------

#-----------
# FUNCTIONS
#-----------
cat(paste0("\n","><> Prepare JABBA prior inputs <><","\n"))

#--------------------------------------------------
# Function to get beta prior parameters
#--------------------------------------------------
get_beta <- function(mu,CV,Min=0,Prior="x"){
  a = seq(0.0001,1000,0.001)
  b= (a-mu*a)/mu
  s2 = a*b/((a+b)^2*(a+b+1))
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = (a-mu*a)/mu
  x = seq(Min,1,0.001)  
  pdf = dbeta(x,a,b)  
  plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(a,b))
}

#--------------------------------------------------
# Function to get gamma prior parameters
#--------------------------------------------------

get_gamma <- function(mu,CV,Prior="x"){
  a = seq(0.0001,1000,0.0001)
  b = a/mu
  s2 = (a/b^2)
  sdev = sqrt(s2)
  # find beta )parameter a
  CV.check = (sdev/mu-CV)^2
  a = a[CV.check==min(CV.check)]
  #find beta parameter b
  b = a/mu
  x = sort(rgamma(1000,a,b))  
  pdf = dgamma(x,a,b)  
  plot(x,pdf,type="l",xlim=range(x[pdf>0.01]),xlab=paste(Prior),ylab="",yaxt="n")
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(a,b))
}


#--------------------------------------------------
# Function to get lognormal prior parameters
#--------------------------------------------------
plot_lnorm <- function(mu,CV,Prior="x"){
  sdev= sqrt(log(CV^2+1))
  rand.pr = rlnorm(1000,log(mu),sdev)
  x = seq(min(rand.pr),quantile(rand.pr,0.995),max(rand.pr/500))  
  pdf = dlnorm(x,log(mu),sdev)  
  plot(x,pdf,type="l",xlim=range(x),xlab=paste(Prior),ylab="",yaxt="n")
  polygon(c(x,rev(x)),c(rep(0,length(x)),rev(ifelse(pdf==Inf,100000,pdf))),col="grey")
  return(c(mu,sdev))
}

#------------------------------------
# Function kobeJabba for FLR
#------------------------------------
kobeJabba<-function(x,minyear=1){
  
  out=cbind(melt(x[,,2]),c(x[,,3]))
  names(out)=c("iter","year","stock","harvest")
  out$year=out$year+minyear-1
  out}

#-------------------------------------------------
# Function kobeJabbaProj for projections with FLR
#-------------------------------------------------
kobeJabbaProj<-function(x,minyear=1,tac=NULL){
  
  out=cbind(melt(x[,,,2]),c(x[,,,3]))
  names(out)=c("iter","year","tac","stock","harvest")
  out$year=out$year+minyear-1
  
  out}

#----------------------------------------------------
# Determine initial ranges for r
#----------------------------------------------------
if(r.dist=="range"){
  # initial range of r based on resilience (FishBase.org)
  if(length(r.prior)>1){ start.r = r.prior} else
    if(r.prior == "High") {
      start.r <- c(0.6,1.5)} else if(r.prior == "Medium") {
        start.r <- c(0.2,0.8)}    else if(r.prior == "Low") {
          start.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
            start.r <- c(0.015,0.1)}  
  
  log.r = mean(log(start.r))
  sd.r = abs(log.r - log(start.r[1]))/2
  r.prior = c(exp(log.r),sd.r)  
  CV.r = sqrt(exp(sd.r^2)-1) 
} else {
  log.r = log(r.prior[1])
  sd.r = r.prior[2]
  CV.r = sqrt(exp(sd.r^2)-1) 
}

#----------------------------------------------------
# Prepare K prior
#----------------------------------------------------
if(K.dist=="range"){
  log.K = mean(log(K.prior))
  sd.K= abs(log.K - log(K.prior[1]))/2
  CV.K = sqrt(exp(sd.K^2)-1)
} else {
  
  log.K = log(K.prior[1])
  CV.K = K.prior[2]
  sd.K=sqrt(log(CV.K^2+1))
}


#----------------------------------------------------
# Prepare radius prior
#----------------------------------------------------

if(!is.null(target_rad_mean)){
  log.rad = log(target_rad_mean)
  CV.rad = CV_rad
  sd.rad = sqrt(log(CV.rad^2+1))
}

#----------------------------------------------------------
# Get JABBA parameterization and suplus production function
#----------------------------------------------------------
# For Pella-Tomlinson
if (Model == 3 | Model == 4) {
  ## run thru sensitivities if given shape
  # find inflection point
  ishape = NULL
  # Find shape for  SBmsytoK
  ishape = seq(0.1, 10, 0.001)
  check.shape = ((ishape) ^ (-1 / (ishape - 1)) - BmsyK) ^ 2
  #  Set shape (> 0, with 1.001 ~ Fox and 2 = Schaefer)
  shape =  ishape[check.shape == min(check.shape)]
  
  if (exists("sensname")) {
    if (sensname == 'M') {
      m = shape = sensmean[s] ## becomes m.mu
      shape.CV = sensvar[s]
      cat(paste0(shape, " ", shape.CV), "\n")
    } ## end sensname == M
  } ## end sensname exists
} else {shape = FALSE}


#------------------------------------------------


# Set shape m for Fox and Schaefer: Fox m ~1; Schaefer m =2
if(shape==FALSE){
  if(Model == 1){m=2} else {m = 1.001}}else{m=shape}

cat(paste0("\n","><> Plot Prior distributions in Input subfolder  <><","\n"))

Par = list(mfrow=c(2,2),mai=c(0.5,0.1,0,.1),omi = c(0.1,0.2,0.1,0) + 0.1,mgp=c(2,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/Priors_",assessment,"_",Scenario,".png"), width = 9, height = 9, 
    res = 200, units = "in")
par(Par)
K.pr = plot_lnorm(exp(log.K),CV.K,Prior="K")

if(psi.dist=="beta"){
  psi.pr = get_beta(mu=psi.prior[1],CV=psi.prior[2],Min=0,Prior=paste0("Prior B(",years[1],")/K"))} else {
    psi.pr = plot_lnorm(mu=psi.prior[1],CV=psi.prior[2],Prior=paste0("Prior B(",years[1],")/K"))  
  }


r.pr = plot_lnorm(mu=exp(log.r),CV=CV.r,Prior="r")
mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)

### Added for calculating effective area radius priors ###################
rad.pr = plot_lnorm(mu = log(target_rad_mean), CV_rad, Prior = "Radius")
dev.off() 
##########################################################################

cat(paste0("\n","><> Plot assumed Surplus Production shape in Input subfolder  <><","\n"))


# Plot MSY
Par = list(mfrow=c(1,1),mai=c(0.6,0.3,0,.15),omi = c(0.1,0.2,0.2,0) + 0.1,mgp=c(2,1,0), tck = -0.02,cex=0.8)
png(file = paste0(input.dir,"/Production",assessment,"_",Scenario,".png"), width = 6, height = 5, 
    res = 200, units = "in")
par(Par)

# Get Bmsy/B0 as a fucntion of M 
Bmsy=(m)^(-1/(m-1))
P = seq(0.0001,1,0.001) 
SP = ifelse(P>Plim,r.pr[1]/(m-1)*P*(1-P^(m-1)),r.pr[1]/(m-1)*P*(1-P^(m-1))*4*P)
#if(is.null(refBmsy)==TRUE) refBmsy = Bmsy
plot(P,SP/max(SP),type="l",ylab="Relative Yield",xlab="B/B0",lwd=2)
mtext(paste("Relative Yield"), side=2, outer=TRUE, at=0.6,line=1,cex=0.9)
legend("topright",c("SPM"),col=c(1),lwd=2,bty="n")  


if(Model==4){
  # shape density
  #dm = dgamma(seq(0.001,5,0.1),5,5)*m
  dm = dlnorm((seq(0.001,5,0.1)),log(m),shape.CV)
  dm = dm/max(dm)
  bmsyk  = (seq(0.001,5,0.1))^(-1/(seq(0.001,5,0.1)-1))
  
  polygon(c(bmsyk,rev(bmsyk)),c(dm,rep(0,length(dm))),col="grey",border=0)  
}
abline(v=Bmsy,lty=2)
mtext(paste("Relative Yield"), side=2, outer=TRUE, at=0.6,line=1,cex=0.9)
legend("topright",c("SPM"),col=c(1),lwd=2,bty="n")  
abline(v=Bmsy,lty=2)
dev.off()


# Note PRIORS and save input subfolder
Priors =rbind(K.pr,psi.prior,c(r.pr[1],CV.r))
row.names(Priors) = c("K","Psi","r")
colnames(Priors) = c("Mean","CV")                          
write.csv(Priors,paste0(input.dir,"/Priors",assessment,"_",Scenario,".csv"))


#----------------------------------------------------------
# Set up JABBA 
#----------------------------------------------------------
cat(paste0("\n","><> Set up JAGS input <><","\n"))
# Plot MSY
# remove scientific numbers
options(scipen=999)
#----------------------------------------------------------


# starting values
nq = length(unique(sets.q))
nvar = length(unique(sets.var))

## MDFitchett change to TAC setup
#----------------------------------------------------------
# Setup TAC projection
#---------------------------------------------------------
if(Projection==TRUE) {
  nTAC = length(TACs)
  TAC = mat.or.vec(pyrs, nTAC)
  yr.now = as.numeric(format(Sys.Date(), "%Y")) + 1
  yr.last = max(years) # assessment year
  
  for (i in 1:nTAC) {
    TAC[, i] = c(rep(TC[n.years], yr.now - yr.last), rep(TACs[i], pyrs - (yr.now -yr.last)))
  }
  
} else if(Projection == FALSE){
  nTAC = 1
  TAC = TC[n.years]
  pyrs = 1
}

#---------------------------------------------------------------
# JABBA Schaefer/Fox Models 1-2, Pella 3
#---------------------------------------------------------------

# Slope of hockey-stick
slope.HS = ifelse(Plim==0,1/10^-10,1/Plim)


nSel = 1 # setup for JABBA-SELECT version (in prep)
nI = ncol(CPUE) # number of CPUE series
stI = ifelse(proc.dev.all == TRUE, 1, c(1:n.years)[is.na(apply(CPUE,1,mean,na.rm = TRUE)) == FALSE][1]) #first year with CPUE


# Initial starting values
inits <- function(){list(K= rlnorm(1,log.K,0.3),
                         q = c(runif(nq-1,0.005,0.5),NA),   # added NA for 3rd q which is the BFISH survey catchability since it is not stochastic
                         isigma2.est=runif(1,20,100), 
                         itau2=runif(nvar,80,200),
                         rad = runif(1, 10, 60))}
# starting value option
if(init.values==TRUE){
inits <- function(){list(K= K.init,
                         r=r.init,
                         q = q.init, 
                         isigma2.est=runif(1,20,100), 
                         itau2=runif(nvar,80,200))}
}

# JABBA input data 
surplus.dat <- list(
  N = n.years,
  TC = TC,
  I = CPUE,
  SE2 = se2,
  mu.m = m,
  r.pr = r.pr,
  psi.pr = psi.pr,
  K.pr = K.pr,
  rad.pr = rad.pr, #added for effective area radius prior
  n.grid = n.grid, #added for number of sampling grids in sampling domain
  a.grid = a.grid, #added for area of sampling grid
  # nq=nq,   #removed because not using a for-loop for q priors anymore
  s_lambda = s_lambda, #added for BFISH survey weighting
  nI = nI,
  nvar = nvar,
  sigma.fixed = ifelse(sigma.proc == TRUE, 0, sigma.proc),
  sets.var = sets.var,
  sets.q = sets.q,
  pen.bk = rep(0, n.years),
  Plim = Plim,
  slope.HS = slope.HS,
  nTAC = nTAC,
  pyrs = pyrs,
  #Syears = syrs, #added for loop to fit BFISH biomass, survey years
  TAC = TAC,
  igamma = igamma,
  stI = stI,
  TACint = TACint,
  P_bound = P_bound,
  proc.pen = 0,
  K.pen = 0,
  obs.pen = rep(0, nvar),
  q_bounds = q_bounds,
  sigmaobs_bound = sigmaobs_bound,
  sigmaproc_bound = sigmaproc_bound,
  K_bounds = K_bounds
)
# If shape parameter is estimated (Model =4)
if(Model==4){
  surplus.dat$m.CV = shape.CV }


# JAGS model file
JABBA = "JABBA.jags"


# PARAMETERS TO MONITOR
params <- c("K",
            "r", 
            "q", 
            "psi",
            "sigma2", 
            "tau2",
            "m",
            "Hmsy",
            "SBmsy", 
            "MSY", 
            "BtoBmsy",
            "HtoHmsy",
            "CPUE",
            "Proc.Dev",
            "P",
            "SB",
            "prP",
            "prBtoBmsy",
            "prHtoHmsy",
            "TOE", 
            "rad")


cat(paste0("\n","><> RUN ",Mod.names," model for ",assessment," ",Scenario," in JAGS <><","\n","\n"))

# JAGS MODEL Standard
sink("JABBA.jags")
cat("
    
    model {
    
    # Prior specifications  
    eps <- 0.0000000000000000000000000000000001 # small constant    
    
    
    ### BFISH Effective Sampling Area Radius Prior ################################
    
    #rad_precision <- 1.0/log(1.0 + CV_rad * CV_rad)
    #rad_mean <- log(target_rad_mean) - (0.5/rad_precision)
    #rad ~ dlnorm(rad_mean, rad_precision)
    ### rad.pr was calculated using the plot_lnorm function that was used to generate lognormal priors for K and r 
    ### see line 502 that is why code above (which was taken from BUGS code) was omitted and replaced with line below
    rad ~ dlnorm(rad.pr[1], pow(rad.pr[2],-2))
    
    
    #Catchability coefficients
   
    q[1] ~ dunif(q_bounds[1],q_bounds[2]) 
    q[2] ~ dunif(q_bounds[1],q_bounds[2]) 
    q[3] <- a.grid/(rad*rad*3.14159)
   
   ## removed for loop for q priors because having a deterministic q was causing problems 
    # for(i in 1:nq-1)
    # {   
    # q[i] ~ dunif(q_bounds[1],q_bounds[2])    
    # }  
    
    
    ")
   
if(Model==4){
  cat("
      # Shape m prior
      m ~ dlnorm(log(mu.m),pow(m.CV,-2))
      ",append=TRUE)  
  }else{ cat(" 
      m <- mu.m     
    ",append=TRUE)}

if(psi.dist =="beta"){
  cat("
      # Beta Prior for Biomass depletion at the start (deteministic)
      psi ~ dbeta(psi.pr[1],psi.pr[2])
      ",append=TRUE)
} else {
  cat("
      # Lognormal for Biomass depletion at the start (deteministic)
      psi ~ dlnorm(log(psi.pr[1]),pow(psi.pr[2],-2)) #I(0.1,1.1)    
      ",append=TRUE)  
}

if(sigma.proc==TRUE){
  cat("
      # Process variance
      isigma2 <- isigma2.est 
      sigma2 <- pow(isigma2,-1)
      sigma <- sqrt(sigma2)
      fakesigma.fixed <- sigma.fixed # Prevent unused variable error msg    
      ",append=TRUE)  
}else{ cat(" 
      isigma2 <- pow(sigma.fixed+eps,-2) 
           sigma2 <- pow(isigma2,-1)
           sigma <- sqrt(sigma2)
           
           ",append=TRUE)}

if(sigma.est==TRUE){
  cat("
      # Obsevation variance
      for(i in 1:nvar)
      {
      # Observation error
      itau2[i]~ dgamma(0.2,1.0)
      tau2[i] <- 1/itau2[i]
      }
      
      for(i in 1:nI)
      {
      for(t in 1:N)
      {
      var.obs[t,i] <- SE2[t,i]+tau2[sets.var[i]] 
      ivar.obs[t,i] <- 1/var.obs[t,i]
      # note total observation error (TOE)     
      TOE[t,i] <- sqrt(var.obs[t,i]) # Total observation variance
      
      }}
      ",append=TRUE)  
}else{ cat(" 
      # Obsevation variance
           for(i in 1:nvar)
           {
           # Observation error
           itau2[i]~ dgamma(4,0.01)
           tau2[i] <- 1/itau2[i]
           }
           
           for(i in 1:nI)
           {
           for(t in 1:N)
           {
           var.obs[t,i] <- SE2[t,i] # drop tau2
           fake.tau[t,i] <- tau2[sets.var[i]]
           
           ivar.obs[t,i] <- 1/var.obs[t,i]
           # note total observation error (TOE)     
           TOE[t,i] <- sqrt(var.obs[t,i])
           
           }}
           
           ",append=TRUE)}

# Run rest of code  
cat("  
    # Process variance prior
    isigma2.est ~ dgamma(igamma[1],igamma[2])
    
    # Carrying Capacity SB0
    K ~ dlnorm(log(K.pr[1]),pow(K.pr[2], -2))
    
    # informative priors for Hmsy as a function of r
    r ~ dlnorm(log(r.pr[1]),pow(r.pr[2],-2))
    
    
    #Process equation
    Pmean[1] <- log(psi)
    iPV[1] <- ifelse(1<(stI),10000,isigma2) # inverse process variance
    P[1] ~ dlnorm(Pmean[1],iPV[1]) # set to small noise instead of isigma2
    
    penB[1]  <- ifelse(P[1]<P_bound[1],log(K*P[1])-log(K*P_bound[1]),ifelse(P[1]>P_bound[2],log(K*P[1])-log(K*P_bound[2]),0)) # penalty if Pmean is outside viable biomass
    
    # Process equation
    for (t in 2:N) 
    {
    Pmean[t] <- ifelse(P[t-1] > Plim,
    log(max(P[t-1] +  r/(m-1)*P[t-1]*(1-pow(P[t-1],m-1)) - TC[t-1]/K,0.005)),
    log(max(P[t-1] +  r/(m-1)*P[t-1]*(1-pow(P[t-1],m-1))*P[t-1]*slope.HS - TC[t-1]/K,0.005)))
    iPV[t] <- ifelse(t<(stI),10000,isigma2) # inverse process variance
    P[t] ~ dlnorm(Pmean[t],iPV[t])
    penB[t]  <- ifelse(P[t]<(P_bound[1]),log(K*P[t])-log(K*(P_bound[1])),ifelse(P[t]>P_bound[2],log(K*P[t])-log(K*(P_bound[2])),0)) # penalty if Pmean is outside viable biomass
    }
    

    # Process error deviation 
    for(t in 1:N){
    Proc.Dev[t] <- P[t]-exp(Pmean[t])} 
    
    # Enforce soft penalties on bounds for P
    for(t in 1:N){
    pen.bk[t] ~ dnorm(penB[t],1000) # enforce penalty with CV = 0.1
    }
    
    Hmsy <- r*pow(m-1,-1)*(1-1/m) 
    
    for (t in 1:N) 
    { 
    SB[t] <- K*P[t]    
    H[t] <- TC[t]/SB[t] 
    }
    
    # Observation equation in related to EB
    
    for(i in 1:(nI-1))
    {
    for (t in 1:N) 
    { 
    Imean[t,i] <- log(q[sets.q[i]]*P[t]*K);
    I[t,i] ~ dlnorm(Imean[t,i],(ivar.obs[t,i]));
    CPUE[t,i] <- q[sets.q[i]]*P[t]*K
    }}
    
     for(i in nI){ ## added by Meg for BFISH index
      for (t in 1:N){ 

        Imean[t,i] <- log(P[t] * K/(q[sets.q[i]])*n.grid)
        survey_precision[t] <- (s_lambda*s_lambda)/( ivar.obs[t,i])
        I[t, i] ~ dlnorm(Imean[t,i], survey_precision[t])
        CPUE[t,i] <- P[t]*K/(q[sets.q[i]]*n.grid)

    }
    }
    
    #Management quantities
    SBmsy_K <- (m)^(-1/(m-1))
    SBmsy <- SBmsy_K*K
    
    MSY <- SBmsy*Hmsy
    for (t in 1:N)
    {
    # use x y to put them towards the end of the alphabetically sorted  mcmc object
    #SP[t] <- pow(r.pella,-(m-1))*SB[t]*(1-pow(P[t],m-1))
    BtoBmsy[t] <- SB[t]/SBmsy
    HtoHmsy[t] <- H[t]/(Hmsy) 
    }
    
    
    # Enforce soft penalty on K if < K_bounds >  
    K.pen ~ dnorm(penK,1000) # enforce penalty 
    penK  <- ifelse(K<(K_bounds[1]),log(K)-log(K_bounds[1]),ifelse(K>K_bounds[2],log(K)-log(K_bounds[2]),0)) # penalty if Pmean is outside viable biomass
    

    # Enforce soft penalty on process deviance if sigma.proc > 0.2 
    proc.pen ~ dnorm(penProc,1000) # enforce penalty 
    penProc  <- ifelse(sigma>sigmaproc_bound,log(sigma)-log(sigmaproc_bound),0) 
     
    
    # Enforce soft penalty on observation error if sigma.obs > sigma_bound 
    for(i in 1:nvar){
    obs.pen[i] ~ dnorm(penObs[i],1000) # enforce penalty 
    penObs[i]  <- ifelse(pow(tau2[i],0.5)>sigmaobs_bound,log(pow(tau2[i],0.5))-log(sigmaobs_bound),0) 
    }
    
    ", append=TRUE)

# PROJECTION
if(Projection==TRUE){
  cat("
      for(i in 1:nTAC){
      # Project first year into the future
      prPmean[1,i] <- ifelse(P[N] > Plim,
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1)) - TACint/K,0.005)),
      log(max(P[N] +  Hmsy/(1-1/m)*P[N]*(1-pow(P[N],m-1))*4*P[N] - TACint/K,0.005)))
      prP[1,i] ~ dlnorm(prPmean[1,i],isigma2) 
      # Project all following years
      for(t in 2:pyrs){
      prPmean[t,i] <- ifelse(prP[t-1,i] > Plim,
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1)) - TAC[t-1,i]/K,0.001)),
      log(max(prP[t-1,i] +  Hmsy/(1-1/m)*prP[t-1,i]*(1-pow(prP[t-1,i],m-1))*slope.HS*prP[t-1,i] - TAC[t-1,i]/K,0.005)))
      # process error (as monte-carlo simular)
      prP[t,i] ~ dlnorm(prPmean[t,i],isigma2)}
      for(t in 1:pyrs){
      prB[t,i] <- prP[t,i]*K
      prH[t,i] <- TAC[t,i]/prB[t,i]
      prHtoHmsy[t,i] <- prH[t,i]/Hmsy
      prBtoBmsy[t,i] <- prB[t,i]/SBmsy
      }}  
      ",append=TRUE)} else {
        cat("
            #Prevent error for unused input    
            fakeTAC <-  TAC
            fakepyrs <- pyrs 
            fakenTAC <- nTAC    
            fakeTACint <- TACint
            prHtoHmsy <- 1
            prP <- 1 
            prBtoBmsy <- 1    
            ", append=TRUE)}  

cat("
    
} # END OF MODEL
    ",append=TRUE,fill = TRUE)
sink()


ptm <- proc.time()

mod <- jags(surplus.dat, 
            inits, 
            params, 
            paste(JABBA), 
            n.chains = nc, 
            n.thin = nt, 
            n.iter = ni, 
            n.burnin = nb)  # adapt is burn-in

proc.time() - ptm
save.time = proc.time() - ptm

cat(paste0("\n",paste0("><> Scenario ",Scenario,"_",Mod.names," completed in ",as.integer(save.time[3]/60)," min and ",round((save.time[3]/60-as.integer(save.time[3]/60))*100)," sec <><","\n")))

cat(paste0("\n","><> Produce results output of ",Mod.names," model for ",assessment," ",Scenario," <><","\n"))

# if run with library(rjags)
posteriors = mod$BUGSoutput$sims.list

#-----------------------------------------------------------
# <><<><<><<><<><<><<>< Outputs ><>><>><>><>><>><>><>><>><>
#-----------------------------------------------------------
output.dir = paste0(File,"/",assessment,"/",Scenario,"_",Mod.names,"/Output")
dir.create(output.dir, showWarnings = FALSE)

# run some mcmc convergence tests
par.dat <- data.frame(posteriors[params[c(1:7,21)]])  #added radius paramater to outputs
geweke <- geweke.diag(data.frame(par.dat))
pvalues <- 2*pnorm(-abs(geweke$z))

heidle <- heidel.diag(data.frame(par.dat))

# postrior means + 95% BCIs
#Model  parameter
apply(par.dat,2,quantile,c(0.025,0.5,0.975))

man.dat <- data.frame(posteriors[params[8:10]])
#Management quantaties
apply(man.dat,2,quantile,c(0.025,0.5,0.975))

# Depletion
Depletion <- posteriors$P[,c(1,n.years)]
colnames(Depletion) <- c(paste0("P", years[1]), paste0("P", years[n.years]))

# Current stock status (Kobe posterior)
H_Hmsy.cur <- posteriors$HtoHmsy[, c(n.years)]
B_Bmsy.cur <- posteriors$BtoBmsy[, c(n.years)]


# Prepare posterior quantaties
man.dat <- data.frame(man.dat, Depletion, B_Bmsy.cur, H_Hmsy.cur)

results <- round(t(cbind(apply(par.dat, 2, quantile, c(0.025, 0.5, 0.975)))), 6)

results.mean <- round(cbind(apply(par.dat, 2, mean)), 6)

results <- data.frame(Median = results[, 2], 
                      LCI = results[, 1], 
                      UCI = results[, 3], 
                      Geweke.p = round(pvalues, 3), 
                      Heidel.p = round(heidle[, 3], 3), 
                      Mean = results.mean[, 1])

ref.points <- round(t(cbind(apply(man.dat, 2, quantile, c(0.025, 0.5, 0.975)))), 3)

ref.points <- data.frame(Median = ref.points[, 2], 
                         LCI = ref.points[, 1], 
                         UCI = ref.points[, 3])

# get number of parameters
npar <- length(par.dat)
# number of years
N <- n.years

# Save posteriors (Produces large object!)
if (save.all == TRUE){
  
  save(posteriors, file = paste0(output.dir, "/", Scenario, "_posteriors"))

  }

#-------------------------------------------------------------------------
# Save parameters, results table and current status posterior in csv files
#-------------------------------------------------------------------------

# Save model estimates and convergence p-values
write.csv(data.frame(results), paste0(output.dir, "/Estimates_", assessment, "_", Scenario, ".csv"))


############## THIS SECTION SEEMS SUSPECT. ITS CHANGING THE NAME OF q.2 TO SIGMA.PROC WHICH DOESN'T SEEM CORRECT. COMMENTING OUT FOR NOW. WILL CHECK WITH FELILPE OR JOHN. 3/18/2022 ##############################################

# Make standard results table with parameter estimates and reference points
# Table <- rbind(data.frame(results)[, 1:3], data.frame(ref.points))
# 
# Table[4, ] <- round(sqrt((Table[4, ])), 3)
# 
# rownames(Table)[4] <- "sigma.proc"
# 
# write.csv(Table, paste0(output.dir, "/Results_", assessment, "_", Scenario, ".csv"))
################################################################################################################

# Save posterior of recent assessment year (KOBE posterior)
write.csv(data.frame(BtoBmsy = B_Bmsy.cur, 
                     FtoFmsy = H_Hmsy.cur), 
          paste0(output.dir, "/Status_posterior", assessment, ".csv"))

## source all plotting scripts
source(paste0(JABBA.file, "/plot_JABBA.R"))


if (save.trajectories == TRUE) {
  cat(paste0("\n", "><> Saving Posteriors of FRP trajectories <><", "\n"))

  # FRP trajectories
  trajectories <- array(NA, c(nsaved, n.years, 3))
  trajectories[, , 1] <- posteriors$P
  trajectories[, , 2] <- posteriors$BtoBmsy
  trajectories[, , 3] <- posteriors$HtoHmsy


  kb <- kobeJabba(trajectories, years[1])
  save(kb, file = paste0(output.dir, "/", Scenario, "_trajectories"))
}

cat(paste0("\n", "><> Scenario ", Mod.names, "_", Scenario, " for ", assessment, " - DONE! <><", "\n"))




