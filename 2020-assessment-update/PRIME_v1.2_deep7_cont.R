## ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><><
## JABBA: Just Another Bayesian Biomass Assessment
## PRIME ('Input' File)
## Developed by Henning Winker & Felipe Carvalho (Cape Town/Hawaii)
##
## Modified from Guam_base_FINAL files by megumi.oshima@noaa.gov
## for use in the 2023 Deep 7 Bottomfish assessment
## Now includes:
##  variability around catch
##
## Last updated Feb. 24, 2022
## ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>


# required packages
library(gplots)
library(coda)
library(rjags)
library(R2jags)
library(fitdistrplus)
library(reshape)
library(here) # use with Rprojects
library(nmfspalette)

# Set Working directory file, where assessments are stored
File <- here("2020-assessment-update")
# Set working directory for JABBA R source code
JABBA.file <- here("2020-assessment-update")
# JABBA version
version <- "v1.2_deep7_cont"
# Set Assessment file: assement folder within File that includes .csv input files (suffix included in data file names)
assessment <- "deep7_continuity"


#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Graphic, Output, Saving (.RData) settings
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
KOBE.plot <- TRUE # Produces JABBA Kobe plot
KOBE.type <- c("ICCAT", "IOTC")[2] # ICCAT uses 3 colors; IOTC 4 (incl. orange)
Biplot <- TRUE # Produces a "post-modern" biplot with buffer and target zones (Quinn & Collie 2005)
SP.plot <- c("standard", "phase")[2] # Produces standard or 'Kobe phase' SP plot
save.trajectories <- TRUE # saves posteriors of P=B/K, B/Bmsy and H/Hmsy as .RData object
harvest.label <- c("Hmsy", "Fmsy")[1] # choose label preference H/Hmsy versus Fmsy
CPUE.plot <- TRUE # Runs state-tool to produce "aligned" multi-CPUE plot
meanCPUE <- FALSE # Uses averaged CPUE from state-space tool instead of individual indices
Projection <- FALSE # Use Projections: requires to define TACs vectors
save.projections <- FALSE # saves projection posteriors as .RData object
catch.metric <- "000 t" # Define catch input metric e.g. (tons) "000 t" etc
Reproduce.seed <- TRUE # If FALSE a random seed assigned to each run, if TRUE set.seed(123)
# P_bound = c(0.02,1.2)  # Soft penalty bounds for P
# Save entire posterior as .RData object
save.all <- TRUE # (if TRUE, a very large R object of entire posterior is saved)
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Optional: Note Scenarios
#><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
# Specify Scenario name for output file names
Scenarios <- c("Cont", NA, NA, NA)

# Execute multiple JABBA runs in loop

for (s in 1:1) {
  Scenario <- Scenarios[s]

  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Suplus Production model specifications
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

  # Choose model type:
  # 1: Schaefer
  # 2: Fox
  # 3: Pella-Tomlinsson
  # 4: Pella-Tomlinsson M

  Model <- c(1, 1, 1, 1)[s]
  Mod.names <- c("Schaefer", "Fox", "Pella", "Pella_m")[Model]

  # Depensation opiton:
  # Set Plim = Blim/K where recruitment may become impaired (e.g. Plim = 0.25)
  # Choose Plim = 0 to reduce to conventional Schaefer, Fox, Pella models
  Plim <- 0

  # Required specification for Pella-Tomlinson (Model = 3/4)
  BmsyK <- 0.5 # Set Surplus Production curve inflection point
  shape.CV <- sqrt(log(1 + 0.5^2)) # Need sd on log-scale, not CV # Must be defined if Model = 4!
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

  #--------------------------------------------------
  # Read csv files
  #--------------------------------------------------

  # Use SE from csv file for survey indices, and CVs from CPUE indices
  SE.I <- TRUE

  # Load assessment data
  ### Catch data (contains reported and unreported catch)
  catch <- read.csv(paste0(File, "/", assessment, "/catch", assessment, ".csv"))

  ### CPUE data
  cpue <- read.csv(paste0(File, "/", assessment, "/cpue", assessment, ".csv"))
  

  if(SE.I == TRUE) {
    
    se <- read.csv(paste0(File, "/", assessment, "/se", assessment, ".csv"))
    
    ## convert mean CV to log(se)
    se$CV_1 <- sqrt(log(1 + se$CV_1^2))
    se$CV_2 <- sqrt(log(1 + se$CV_2^2))
    
  }
  #######################################################
  #######################################################
  
  
  indices2 <- names(cpue)[-1]
  wink.colors <- data.frame(
    idx = indices2,
    cols = nmfs_palette("regional web")(5)[1:length(indices2)]
  )


  #------------------------------------------------------
  # Option use mean CPUE from state-space cpue averaging
  #-----------------------------------------------------
  meanCPUE <- FALSE

  #-----------------------------------------------------
  # Starting value option for r and K
  #-----------------------------------------------------

  #------------------------------------------------
  # Prior for unfished biomass K
  #------------------------------------------------
  # The option are:
  # a) Specify as a lognormal prior with mean and CV
  # b) Specify as range to be converted into lognormal prior

  # ><> new objective K prior
  K.dist <- c("lnorm", "range")[1] # ><> to range
  #  K.prior
  K.prior <- c(29, 0.5) # actual CV needed here (mean = 29, CV = 50%)

  #-----------------------------------------------------------
  # mean and CV and sd for Initial depletion level P1= SB/SB0
  #-----------------------------------------------------------
  # Set the initial depletion prior B1/K
  # To be converted into a lognormal prior (with upper bound at 1.1)

  psi.dist <- c("lnorm", "beta")[1]
  # specify as mean and CV
  # psi.prior = c(psi.prior.mean,0.5)
  psi.prior <- c(0.53, 0.2) # actual CV needed here
  #--------------------------------------------------------------
  # Determine estimation for catchability q and observation error
  #--------------------------------------------------------------
  # Assign q to CPUE
  sets.q <- 1:(ncol(cpue) - 1)

  #----------------------------------------------------
  # Determine r prior
  #----------------------------------------------------
  # The option are:
  # a) Specifying a lognormal prior
  # b) Specifying a resiliance category after Froese et al. (2017; CMSY)
  # Resilience: "Very low", "Low", "Medium", High" (requires r.range = TRUE)

  # use [1] lognormal(mean,stdev) or [2] range (min,max) or
  r.dist <- c("lnorm", "range")[1]
  #  r.prior = c(0.2735,0.3)
  r.prior <- c(0.10, sqrt(log(1 + 0.25^2))) # needs to be sd on log-scale
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Observation  Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>

  # To Estimate additional observation variance set sigma.add = TRUE
  sigma.est <- TRUE

  # Series
  sets.var <- 1:(ncol(cpue) - 1) # estimate individual additional variance

  # As option for data-weighing
  # minimum fixed observation error for each variance set (optional choose 1 value for both)
  fixed.obsE <- c(0.0) # Important if SE.I is not available

  # Total observation error: TOE = sqrt(SE^2+sigma.est^2+fixed.obsE^2)

  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Process Error
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Estimate set sigma.proc == True
  sigma.proc <- TRUE
  # Determines if process error deviation are estimated for all years (TRUE)
  # or only from the point the first abundance index becomes available (FALSE)
  proc.dev.all <- TRUE
  #------------------------------------------
  if (sigma.proc == TRUE) {
    igamma <- c(2, 0.1) # specify inv-gamma parameters  WAS 4,0.01

    # Process error check
    gamma.check <- 1 / rgamma(1000, igamma[1], igamma[2])
    # check mean process error + CV
    # mu.proc = sqrt(mean(gamma.check)); CV.proc = sd(sqrt(gamma.check))/mean(sqrt(gamma.check))

    # check CV
    # round(c(mu.proc,CV.proc),3)
    # quantile(sqrt(gamma.check),c(0.1,0.9))
  } else {
    sigma.proc <- 0.1 # IF Fixed: typicallly 0.05-0.15 (see Ono et al. 2012)
  }
  #--------------------------------------------

  #><><><><><><><><><><><><><><><><><><>>
  # Include variation in catch estimates
  #><><><><><><><><><><><><><><><><><><>>
  Catch.CV <- FALSE

  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  # Optional: Do TAC Projections
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>>
  Projection <- FALSE # Switch on by Projection = TRUE

  # Set range for alternative TAC projections
  TACs <- seq(10000, 160000, 15000)

  # TACint = mean(catch[nrow(catch)-3,2]:catch[nrow(catch),2]) # avg last 3 years
  TACint <- catch[nrow(catch), 2] ## Catch 2016 is initital

  # Set number of projections years
  pyrs <- 10

  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
  # Execute model and produce output
  #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

  # MCMC settings
  ni <- 55500 # Number of iterations
  nt <- 20 # steps saved
  nb <- 25500 # Burn-in
  nc <- 2 # number of chains
  nsaved <- (ni - nb) / nt * nc

  source(paste0(JABBA.file, "/JABBA", version, ".R"))
} # THE END
