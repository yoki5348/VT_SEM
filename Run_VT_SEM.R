
# Parameters for simulation dataset generation
num.cvs <- 100 #The number of causal variants
am <-  0.25 # assortative mating parameter, mu
herit0 <- 0.5 #Inital heritability
vf0 <- 0.15 #Initial variance of vertical transmitted environment factor
######ve0=1-vf0-herit0, so herit0+vf0 must be less than 1
Unequal_AM=TRUE # TRUE: Generating disequilibrium AM dataset, FALSE: Generating equilibrium AM dataset
prop.h2.latent <- 0.5 #Proportion of heritability 

seed <- 1234567
max.cores <- 1 ##Number of used cores for OpenMx
num.gen <- 20 #Number of generating generation
pop.size <- 1000 #The number of independent trio dataset
num.its <- 1000
RUN.MARKERS <- FALSE  #whether to only consider GRMs built from CVs (FALSE) or both CVs and SNPs (TRUE)
avoid.inb <- FALSE
save.covariances <- TRUE
save.history <- TRUE

#AM Simulation wildcards
MIN.MAF <- .1
MAX.MAF <- .50

#CVs and their effect sizes
MAF.VECTOR <- runif(num.cvs,MIN.MAF,MAX.MAF)  #Can change the distribution of MAFs here
GENTP.VAR <- MAF.VECTOR*(1-MAF.VECTOR)*2
#ALPHA.VECTOR <- rnorm(CVS,0,sqrt(1/(CVS*GENTP.VAR))) #Can change the distribution of effect sizes here - random effects
ALPHA.VECTOR <- sample(c(-1,1),num.cvs,replace=TRUE)*sqrt(1/(num.cvs*GENTP.VAR)) #Can change the distribution of effect sizes here - fixed f'n of MAF
CV.INFO <- data.frame(MAF=MAF.VECTOR,alpha=ALPHA.VECTOR) #we'll use this for both the observed and latent

####### Call packages ##########
library(MASS)
library(data.table)
require(OpenMx)
library(stringr)


#################
#1 Functions & source functions

#Assortative mating simulation functions
source("./SIMULATE_DAT_GEN.R")
source("./VT_SEM_functions.R")


AM.DATA <- AM.SIMULATE(CV.INFO=CV.INFO, H2.T0=herit0, NUM.GENERATIONS=num.gen, POP.SIZE=pop.size*3, MATE.COR=am, AVOID.INB=avoid.inb, SAVE.EACH.GEN=save.history, SAVE.COVS=save.covariances, SEED=seed, VF.T0=vf0,PROP.H2.LATENT=prop.h2.latent,Unequal_AM=Unequal_AM) #Dataset generation

##########Run models#######

  #Add important options to openmx
  mxOption(NULL, 'Number of Threads', max.cores) #alternatively, you specify the number of cores yourself
  mxOption(NULL,"Default optimizer","NPSOL")
  mxOption(NULL,"Calculate Hessian","Yes")
  mxOption(NULL,"Standard Errors","Yes")
  
  options()$mxOptions$'Number of Threads'
  options()$mxOptions$'Default optimizer'
  
  #Optimizer issues - make tolerance smaller to increase NPSOL precision - see ?mxOption
  mxOption(NULL,"Feasibility tolerance","1e-5")
  #mxOption(NULL,"Analytic Gradients","No")
  
  options()$mxOptions$'Feasibility tolerance'
  #options()$mxOptions$'Analytic Gradients'
  options()$mxOptions$'Gradient step size'  #1e-7
  options()$mxOptions$'Optimality tolerance'  #1e-7
  
  dat_init=as.data.frame(AM.DATA$PHEN)
  dat_init=dat_init[sample(1:nrow(dat_init),pop.size),] ##independent dataset check
  
  model0=perform_SEM_model0(dat_init,"NTMO","TMO","NTPO","TPO","Y")
  model1_eq=perform_SEM_model1_eq(dat_init,"NTMO","TMO","NTPO","TPO","Y")
  model1_dis=perform_SEM_model1_dis(dat_init,"NTMO","TMO","NTPO","TPO","Y")
  model2_eq=perform_SEM_model2_eq(dat_init,"NTMO","TMO","NTPO","TPO","YM","YP","Y")
  model2_dis=perform_SEM_model2_dis(dat_init,"NTMO","TMO","NTPO","TPO","YM","YP","Y")
  model2_eq_NP=perform_SEM_model2_eq_NP(dat_init,"NTMO","TMO","NTPO","TPO","Y",herit0,prop.h2.latent)
  model2_dis_NP=perform_SEM_model2_dis_NP(dat_init,"NTMO","TMO","NTPO","TPO","Y",herit0,prop.h2.latent)
  
  results=rbind(model0$result_summary,model1_eq$result_summary,model1_dis$result_summary,model2_eq$result_summary,model2_dis$result_summary,
                model2_eq_NP$result_summary,model2_dis_NP$result_summary)
  print(results) ##Print results
  
  summary(model0$AFE.Fit) #Check the fitted results 
  
  