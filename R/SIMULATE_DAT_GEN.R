

#The purpose of this script is to simulate SNP data under assortative mating WITHOUT any LD
#This script started with Simulate.AM5.R, but has been modified to be simpler, faster, and not require LD
#Start date: May 1, 2017

#v8 makes the entire thing a single function

##############################################################################
#LOAD PACKAGES REQUIRED
##############################################################################

#Load packages & functions needed
require(MASS)

#Global options for this script - important to run this or script breaks
op <-  options(stringsAsFactors=FALSE)

##############################################################################
#END LOAD PACKAGES REQUIRED
##############################################################################

##############################################################################
#FUNCTIONS
##############################################################################


##########################
#F1 - SHORT FUNCTIONS
is.even <- function(x) {x%%2 == 0}
is.odd <- function(x) {x%%2 != 0}

corn <- function(x,nrow=10){
    x[1:nrow,1:nrow]}

popvar <- function(x){sum((x-mean(x))^2)/length(x)}


# Constructing delta
ddd<-function(a,b,c){
    b^2-4*a*c
}

# Constructing Quadratic Formula
quadformula <- function(a,b,c){
    if(ddd(a,b,c) > 0){ # first case D>0
        x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
        x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
        result = c(x_1,x_2)
    }
    else if(ddd(a,b,c) == 0){ # second case D=0
        x = -b/(2*a)
    }
    else {"There are no real roots."} # third case D<0
    x = NA
}


next.smallest <- function(x,value,return.index=TRUE){
x2 <- sort(x,decreasing=FALSE)
if(return.index){
    ans <- max(which(x2<value))
} else {
    ans <- x2[max(which(x2<value))]
}
return(ans)}

#########################################################
#F2 - ASSORTATIVE MATING FUNCTION

#PHENDATA=PHEN
#MATCOR=PHENO.MATE.CUR
#avoid.inbreeding=AVOID.INB
#POPSIZE=POP.CUR

assort.mate <- function(PHENDATA, MATCOR, POPSIZE, avoid.inbreeding=TRUE){

#probability that each individual marries for generation currun; capped at .95
size.mate <- nrow(PHENDATA)

### MARRIAGEABLE PEOPLE - this section creates an equal number of males and females who will be paired off below
#We need to make the number of marriable people the same in females & males
num.males.mate  <- sum(PHENDATA[,"MALE"])
num.females.mate <- size.mate-num.males.mate

### MATING VIA PRIMARY ASSORTMENT - ANCESTORS
#split up the PHENDATA matrix
males.PHENDATA <- PHENDATA[PHENDATA[,"MALE"]==1,]
females.PHENDATA <- PHENDATA[PHENDATA[,"MALE"]==0,]

#make exactly equal numbers of males & females
if(num.males.mate>num.females.mate) {males.PHENDATA <- males.PHENDATA[- sample(1:nrow(males.PHENDATA),1),]}
if(num.males.mate<num.females.mate) {females.PHENDATA <- females.PHENDATA[- sample(1:nrow(females.PHENDATA),1),]}

#order the males & females by their mating phenotypic value, lowest to highest
males.PHENDATA <- males.PHENDATA[order(males.PHENDATA[,"Y"]),]
females.PHENDATA <- females.PHENDATA[order(females.PHENDATA[,"Y"]),]

#making the template correlation matrix and ordering the two variables
#as it is (empirical=TRUE), the correlation is *exactly* the AM coefficient each generation; may change to FALSE to make it more realistic
#nevertheless, there is a stochastic element to it due to the sorting
template.AM.dist <- mvrnorm(n=nrow(males.PHENDATA),mu=c(0,0),Sigma=matrix(c(1,MATCOR,MATCOR,1),nrow=2),empirical=TRUE)
#rank.template.males <- rank(template.AM.dist[,1])
#rank.template.females <- rank(template.AM.dist[,2])

#marry the males and females according to the assortative mating coefficient earlier inputed
males.PHENDATA <-  males.PHENDATA[rank(template.AM.dist[,1]),]
females.PHENDATA <-  females.PHENDATA[rank(template.AM.dist[,2]),]


###AVOID INBREEDING
if (avoid.inbreeding==TRUE){

#for now, we just drop the inbreeding pairs; an alternative is to remate them
no.sib.inbreeding <- males.PHENDATA[,"Father.ID"]!=females.PHENDATA[,"Father.ID"]
no.cousin.inbreeding <- {males.PHENDATA[,"Fathers.Father.ID"]!=females.PHENDATA[,"Fathers.Father.ID"]&
                      males.PHENDATA[,"Fathers.Father.ID"]!=females.PHENDATA[,"Mothers.Father.ID"]&
                      males.PHENDATA[,"Mothers.Father.ID"]!=females.PHENDATA[,"Fathers.Father.ID"]&
                      males.PHENDATA[,"Mothers.Father.ID"]!=females.PHENDATA[,"Mothers.Father.ID"]&
                      males.PHENDATA[,"Fathers.Mother.ID"]!=females.PHENDATA[,"Fathers.Mother.ID"]&
                      males.PHENDATA[,"Fathers.Mother.ID"]!=females.PHENDATA[,"Mothers.Mother.ID"]&
                      males.PHENDATA[,"Mothers.Mother.ID"]!=females.PHENDATA[,"Fathers.Mother.ID"]&
                      males.PHENDATA[,"Mothers.Mother.ID"]!=females.PHENDATA[,"Mothers.Mother.ID"]}
no.inbreeding <- (no.sib.inbreeding & no.cousin.inbreeding)


### CREATE FINAL MARRIED PAIRS WHO AREN'T INBRED
#for now, we just drop inbred pairs & create the male & female files (both in order of married pairs)
     males.PHENDATA <- males.PHENDATA[no.inbreeding,]
     females.PHENDATA <- females.PHENDATA[no.inbreeding,]

} ###End avoid inbreeding


### CREATE COR.SPOUSES & "SPOUSE.ID" VARIABLE
males.Spouse.ID <- females.PHENDATA[,'ID']
females.Spouse.ID <- males.PHENDATA[,'ID']

#Create EXACTLY number of offspring desired (note it's possible that, rarely, number offspring won't be exactly that desired, if the # it is off initially is greater than or less than the number of rows in males.PHENDATA)
num.offspring <- rpois(nrow(males.PHENDATA),lambda=POPSIZE/nrow(males.PHENDATA))
(nof <- sum(num.offspring))
if (nof < POPSIZE){num.offspring <- num.offspring + sample(c(rep(1,POPSIZE-nof),rep(0,nrow(males.PHENDATA)-POPSIZE+nof)),nrow(males.PHENDATA),replace=FALSE)}
if (nof > POPSIZE){
subtractor <- rep(0,nrow(males.PHENDATA))
(which.pos <- which(num.offspring>0))
(which.subt <- sample(which.pos,nof-POPSIZE,replace=FALSE))
subtractor[which.subt] <- -1
num.offspring <- num.offspring+subtractor}
#summary(num.offspring) #CHECK

#Also put the spousal info into PHENDATA
males.PHENDATA <- cbind(males.PHENDATA,males.Spouse.ID,num.offspring)
females.PHENDATA <- cbind(females.PHENDATA,females.Spouse.ID,num.offspring)

#make a list of what is returned from the function
list(males.PHENDATA=males.PHENDATA,females.PHENDATA=females.PHENDATA)
}
#########################################################


#########################################################
#F3 Reproduce Function

#MATEDATA <- MATES
#GENDATA <- X
#PHENDATA <- PHEN
#CVDATA <- CV.INFO
#OBS.VAR.BV.GEN0 <- obs.var.bv.gen0
#OBS.VAR.E.GEN0 <- obs.var.e1.gen0

reproduce <- function(MATES,XO,XL,PHEN,CV.INFO,CUR.GEN, DELTA.PATH,A.PATH,F.PATH,E.PATH){

#Create male & female data in mating order for only those who had offspring
males.PHEN <- MATES$males.PHENDATA
males.PHEN <- males.PHEN[males.PHEN[,"num.offspring"]>0,]
females.PHEN <- MATES$females.PHENDATA
females.PHEN <- females.PHEN[females.PHEN[,"num.offspring"]>0,]

#Put male & female genotypes in order of those above
males.index <- match(males.PHEN[,"ID"],PHEN[,"ID"])
males.gentp1.obs <- XO[males.index,]
males.gentp1.lat <- XL[males.index,]
females.index <- match(females.PHEN[,"ID"],PHEN[,"ID"])
females.gentp1.obs <- XO[females.index,]
females.gentp1.lat <- XL[females.index,]

#CHECK - should be 1; it is
#bvm <- males.gentp1 %*% CV.INFO$alpha
#cor(bvm,males.PHEN[,"BV"])
#bvf <- females.gentp1 %*% CV.INFO$alpha
#cor(bvf,females.PHEN[,"BV"])

#Create 1 row per new offspring
males.rep <- rep(1:nrow(males.PHEN),times=males.PHEN[,"num.offspring"])
males.gentp.obs <- males.gentp1.obs[males.rep,]
males.gentp.lat <- males.gentp1.lat[males.rep,]
females.rep <- rep(1:nrow(females.PHEN),times=females.PHEN[,"num.offspring"])
females.gentp.obs <- females.gentp1.obs[females.rep,]
females.gentp.lat <- females.gentp1.lat[females.rep,]

#Create matrices for which of the 2 alleles are passed on at random (male haplotypes) - TRANSMITTED, OBSERVED
male.adder.obs <- matrix(sample(x=c(0,1),size=nrow(males.gentp.obs)*ncol(males.gentp.obs),replace=TRUE),nrow=nrow(males.gentp.obs),ncol=ncol(males.gentp.obs))
XM1.obs <- (males.gentp.obs==1)*1
XM1.adder.obs <- male.adder.obs*XM1.obs #give one or other allele IF you are heterozygous
XM2.obs <- (males.gentp.obs==2)*1
males.haps.obs <- XM1.adder.obs+XM2.obs

#Create matrices for which of the 2 alleles are passed on at random (male haplotypes) - TRANSMITTED, LATENT
male.adder.lat <- matrix(sample(x=c(0,1),size=nrow(males.gentp.lat)*ncol(males.gentp.lat),replace=TRUE),nrow=nrow(males.gentp.lat),ncol=ncol(males.gentp.lat))
XM1.lat <- (males.gentp.lat==1)*1
XM1.adder.lat <- male.adder.lat*XM1.lat #give one or other allele IF you are heterozygous
XM2.lat <- (males.gentp.lat==2)*1
males.haps.lat <- XM1.adder.lat+XM2.lat


#Create matrices for the NON-TRANSMITTED alleles (male haplotypes) - NONTRANSMITTED, OBSERVED
male.NT.adder.obs <- (male.adder.obs*-1)+1
XM1.NT.adder.obs <- male.NT.adder.obs*XM1.obs #give one or other allele IF you are heterozygous
males.NT.haps.obs <- XM1.NT.adder.obs+XM2.obs

#Create matrices for the NON-TRANSMITTED alleles (male haplotypes) - NONTRANSMITTED, LATENT
male.NT.adder.lat <- (male.adder.lat*-1)+1
XM1.NT.adder.lat <- male.NT.adder.lat*XM1.lat #give one or other allele IF you are heterozygous
males.NT.haps.lat <- XM1.NT.adder.lat+XM2.lat



#Create matrices for which of the 2 alleles are passed on at random (female haplotypes) - TRANSMITTED, OBSERVED
female.adder.obs <- matrix(sample(x=c(0,1),size=nrow(females.gentp.obs)*ncol(females.gentp.obs),replace=TRUE),nrow=nrow(females.gentp.obs),ncol=ncol(females.gentp.obs))
XF1.obs <- (females.gentp.obs==1)*1
XF1.adder.obs <- female.adder.obs*XF1.obs
XF2.obs <- (females.gentp.obs==2)*1
females.haps.obs <- XF1.adder.obs+XF2.obs

#Create matrices for which of the 2 alleles are passed on at random (female haplotypes) - TRANSMITTED, LATENT
female.adder.lat <- matrix(sample(x=c(0,1),size=nrow(females.gentp.lat)*ncol(females.gentp.lat),replace=TRUE),nrow=nrow(females.gentp.lat),ncol=ncol(females.gentp.lat))
XF1.lat <- (females.gentp.lat==1)*1
XF1.adder.lat <- female.adder.lat*XF1.lat
XF2.lat <- (females.gentp.lat==2)*1
females.haps.lat <- XF1.adder.lat+XF2.lat


#Create matrices for the NON-TRANSMITTED alleles (female haplotypes) - NONTRANSMITTED, OBSERVED
female.NT.adder.obs <- (female.adder.obs*-1)+1
XF1.NT.adder.obs <- female.NT.adder.obs*XF1.obs #give one or other allele IF you are heterozygous
females.NT.haps.obs <- XF1.NT.adder.obs+XF2.obs

#Create matrices for the NON-TRANSMITTED alleles (female haplotypes) - NONTRANSMITTED, LATENT
female.NT.adder.lat <- (female.adder.lat*-1)+1
XF1.NT.adder.lat <- female.NT.adder.lat*XF1.lat #give one or other allele IF you are heterozygous
females.NT.haps.lat <- XF1.NT.adder.lat+XF2.lat


#Create new genotypes for offspring
XO <- males.haps.obs + females.haps.obs
XL <- males.haps.lat + females.haps.lat

#Create 1 row per new offspring from the phenotype matrix
males.PHEN2 <- males.PHEN[males.rep,]
females.PHEN2 <- females.PHEN[females.rep,]

#Create Phenotypes & haplotype PRS's - OBSERVED
AO <- as.numeric(XO %*% CV.INFO$alpha)
TPO <- as.numeric(males.haps.obs %*% CV.INFO$alpha)
TMO <- as.numeric(females.haps.obs %*% CV.INFO$alpha)
NTPO <- as.numeric(males.NT.haps.obs %*% CV.INFO$alpha)
NTMO <- as.numeric(females.NT.haps.obs %*% CV.INFO$alpha)
BV.NT.O <- NTPO+NTMO

#Create Phenotypes & haplotype breeding values - OBSERVED
AL <- as.numeric(XL %*% CV.INFO$alpha)
TPL <- as.numeric(males.haps.lat %*% CV.INFO$alpha)
TML <- as.numeric(females.haps.lat %*% CV.INFO$alpha)
NTPL <- as.numeric(males.NT.haps.lat %*% CV.INFO$alpha)
NTML <- as.numeric(females.NT.haps.lat %*% CV.INFO$alpha)
BV.NT.L <- NTPL+NTML


CURRENT.N <- length(BV.NT.O)

#Calculate latent variables for offspring
F <- F.PATH*males.PHEN2[,'Y'] + F.PATH*females.PHEN2[,'Y']
E <- scale(rnorm(CURRENT.N))

AOy <- DELTA.PATH*AO
ALy <- A.PATH*AL
Fy <- 1*F
Ey <- E.PATH*E
Y <- AOy + ALy + Fy + Ey # AFE Model

#CHECK
#BV2 <- TP+TM #exactly same as BV
#cov(cbind(TP,NTP,TM,NTM,BV,BV.NT))
#var(Ay);var(Fy);var(Ey);var(Y)

#Create relative information for this generation
ID <- sample(1000000:9999999,size=CURRENT.N,replace=FALSE)
Father.ID <- males.PHEN[,"ID"][males.rep]
Mother.ID <- females.PHEN[,"ID"][females.rep]
Fathers.Father.ID <- males.PHEN[,"Father.ID"][males.rep]
Fathers.Mother.ID <- males.PHEN[,"Mother.ID"][males.rep]
Mothers.Father.ID <- females.PHEN[,"Father.ID"][females.rep]
Mothers.Mother.ID <- females.PHEN[,"Mother.ID"][females.rep]
TOT.ID <- cbind(ID,Father.ID,Mother.ID,Fathers.Father.ID,Fathers.Mother.ID,Mothers.Father.ID,Mothers.Mother.ID)
Father.Y <- males.PHEN[,"Y"][males.rep]
Father.F <- males.PHEN[,"F"][males.rep]
Mother.Y <- females.PHEN[,"Y"][females.rep]
Mother.F <- females.PHEN[,"F"][females.rep]



#Create gender for this gen, making an equal number of males & females
if (is.even(CURRENT.N)) {SEXVEC <- c(rep(0,CURRENT.N/2),rep(1,CURRENT.N/2))}
if (is.odd(CURRENT.N)) {SEXVEC <- c(rep(0,floor(CURRENT.N/2)),rep(1,floor(CURRENT.N/2)),0)}
MALE <- sample(SEXVEC,size=CURRENT.N,replace=FALSE)

#Create Phenotype Matrix
PHEN <- cbind(TOT.ID,MALE,
AO,AL,F,E,
AOy,ALy,Fy,Ey,Y,
BV.NT.O,TPO,TMO,NTPO,NTMO,
BV.NT.L,TPL,TML,NTPL,NTML,
Father.Y,Father.F,Mother.Y,Mother.F)

colnames(PHEN) <- c('ID','Father.ID','Mother.ID','Fathers.Father.ID','Fathers.Mother.ID','Mothers.Father.ID','Mothers.Mother.ID','MALE',
'AO','AL','F','E',
'AOy','ALy','Fy','Ey','Y',
'BV.NT.O','TPO','TMO','NTPO','NTMO',
'BV.NT.L','TPL','TML','NTPL','NTML',
'YP','FP','YM','FM')

#return PHEN & genotypes
list(PHEN=PHEN,XO=XO,XL=XL)
}
#########################################################





#########################################################
#F4 - Assortative mating function

#User supplied to function, for checking only (comment out)
#H2.T0 <- .5 #heritability at initial (time 0) generation
#NUM.GENERATIONS <- 10    #number of generations
#POP.SIZE <- 14000
#MATE.COR <- .4
#AVOID.INB <- TRUE #whether to avoid inbreeding
#SAVE.EACH.GEN <- TRUE #whether to save X & PHEN & MATES from each generation
#SEED <- 1234 #default (0) is that seed is random

#CV.INFO=CV.INFO; H2.T0=HERIT0; NUM.GENERATIONS=NUM.GEN; POP.SIZE=POP.SIZE; MATE.COR=AM; AVOID.INB=AVOID.INB; SAVE.EACH.GEN=SAVE.GEN; SEED=SEED

#CV.INFO=CV.INFO; H2.T0=HERIT0; NUM.GENERATIONS=NUM.GEN; POP.SIZE=POP.SIZE; MATE.COR=MATE.COR; AVOID.INB=AVOID.INB; SAVE.EACH.GEN=SAVE.GEN; SAVE.COVS=FALSE; SEED=SEED

# AVOID.INB=TRUE;SAVE.EACH.GEN=FALSE;SAVE.COVS=FALSE;SEED=0

AM.SIMULATE <- function(CV.INFO, H2.T0, NUM.GENERATIONS, POP.SIZE, MATE.COR, AVOID.INB=TRUE, SAVE.EACH.GEN=FALSE, SAVE.COVS=FALSE,SEED=0, VF.T0,PROP.H2.LATENT,Unequal_AM=FALSE){

###################
#A1 IMPLIED VARIABLES
NUM.CVs <- nrow(CV.INFO)
POP.VECTOR <- rep(POP.SIZE,NUM.GENERATIONS)  #population size over time
PHENO.MATE.VECTOR <- rep(MATE.COR,NUM.GENERATIONS)   #phenotypic correlation between mates at each generation
POP.CUR <- POP.SIZE
###################


###################
#A2 Create Genotypes and Phenotypes for GEN0

#Set seed if needed
if (SEED != 0) {
	set.seed(SEED)
}

#Create path coefficients for A, E, and F
(delta <- sqrt(H2.T0*(1-PROP.H2.LATENT)))
(a <- sqrt(H2.T0*PROP.H2.LATENT))
(f <- sqrt(VF.T0/2))
(e <- sqrt(1 - VF.T0 - H2.T0))

#Create Genotypes of GEN0
NUM.GENTPS <- NUM.CVs*POP.CUR
XO <- matrix(rbinom(NUM.GENTPS,size=2,prob=CV.INFO$MAF),nrow=POP.CUR,ncol=NUM.CVs,byrow=TRUE)
XL <- matrix(rbinom(NUM.GENTPS,size=2,prob=CV.INFO$MAF),nrow=POP.CUR,ncol=NUM.CVs,byrow=TRUE)

#Create Phenotypes of GEN0; VA will be 1 at T0; other components will be relative to VA @t0
AO <- XO %*% CV.INFO$alpha  #this is the observed (PRS) breeding value and should have var=1 at t0
AL <- XL %*% CV.INFO$alpha  #this is the latent breeding value and should have var=1 at t0
(var.bvo.gen0 <- as.numeric(var(AO))) #this should be exactly 1 as n -> infinity
(var.bvl.gen0 <- as.numeric(var(AL))) #this should be exactly 1 as n -> infinity

#Create variance components at T0. These follow exactly the path model Jared sent to me on May 26
F <- f*scale(rnorm(POP.CUR)) + f*scale(rnorm(POP.CUR))
E <- scale(rnorm(POP.CUR))

#Create components of Y (the above multiplied by their path coefficients) and Y
AOy <- delta*AO
ALy <- a*AL
Fy <- 1*F
Ey <- e*E
Y <- (AOy + ALy + Fy + Ey)

#Observed variances
(var.e <- var(Ey))
(var.f <- var(Fy))
(var.ao <- var(AOy))
(var.al <- var(ALy))
(var.y <- var(Y))

#Create relative information for GEN0
TOT.ID.VEC <- matrix(sample(10000000:99999999,size=POP.CUR*7,replace=FALSE),nrow=POP.CUR,ncol=7)

#Create gender for GEN0, making an equal number of males & females
if (is.even(POP.CUR)) {SEXVEC <- c(rep(0,POP.CUR/2),rep(1,POP.CUR/2))}
if (is.odd(POP.CUR)) {SEXVEC <- c(rep(0,floor(POP.CUR/2)),rep(1,floor(POP.CUR/2)),0)}
MALE <- sample(SEXVEC,size=POP.CUR,replace=FALSE)

#Create Phenotype Matrix
PHEN <- cbind(TOT.ID.VEC,MALE,AO,AL,F,E,AOy,ALy,Fy,Ey,Y)
colnames(PHEN) <- c('ID','Father.ID','Mother.ID','Fathers.Father.ID','Fathers.Mother.ID','Mothers.Father.ID','Mothers.Mother.ID','MALE','AO','AL','F','E','AOy','ALy','Fy','Ey','Y')
#X and PHEN are what you want to create each generation

#Potentially save GEN0
HISTORY <- list()
if (SAVE.EACH.GEN){
HISTORY$MATES[[1]] <- PHEN
HISTORY$PHEN[[1]] <- PHEN
HISTORY$XO[[1]] <- XO
HISTORY$XL[[1]] <- XL
}

COVS <- list()
if (SAVE.COVS){
COVS[[1]] <- NULL}

NNN <- rep(NA,NUM.GENERATIONS+1)
#Save basic information
SUMMARY.RES <- data.frame(GEN=seq(0,NUM.GENERATIONS),NUM.CVs=NUM.CVs,
                          MATE.COR=MATE.COR,
                          POPSIZE=POP.CUR,
                          VAO=var.ao, #new0 - 
                          VAL=var.al, #new1  - 
                          VF=var.f,
                          VE=var.e,
                          VP=var.y,
                          h2=(var.ao+var.al)/var.y, #new0 - 
                          h2.obs=var.ao/var.y, #new1 - 
                          h2.lat=var.al/var.y, #new1 - 
                          var.hap.obs=NNN, #new0 - 
                          cov.hap.obs=NNN, #new0 - 
                          var.T.obs=NNN,  #new0 - 
                          var.NT.obs=NNN, #new0 - 
                          cov.T.NT.obs=NNN, #new0 - 
                          w=NNN,
                          w.parental.gen=NNN,
                          y.parental.gen=NNN,
                          cov.y.yp=NNN,
                          cov.y.ym=NNN,
                          omega=NNN,
                          theta.NT=NNN,
                          theta.T=NNN,
                          cov.yp.ym=NNN,
                          var.yp=NNN,
                          var.ym=NNN,
                          var.hap.lat=NNN, #new1
                          cov.hap.lat=NNN, #new1
                          var.T.lat=NNN, #new1
                          var.NT.lat=NNN, #new1
                          cov.T.NT.lat=NNN, #new1
                          y=NNN, #new1
                          gamma=NNN, #new1
                          cov.hap.obs.hap.lat=NNN #new1
                          )

###################



##########################
#A3 - Loop through the generations, beginning with generation 1
for (CUR.GEN in  1:NUM.GENERATIONS){

#WILDCARDS THIS GEN
POP.CUR <- POP.VECTOR[CUR.GEN] #scalar - the pop size for the children of this generation
PHENO.MATE.CUR <- PHENO.MATE.VECTOR[CUR.GEN] #scalar - r_mates this generation
##########################


###################
#A4 Assortatively mate the people
#MATES FOR EACH GROUP THIS GENERATION
#PHENDATA=PHEN;  MATCOR=PHENO.MATE.CUR; POPSIZE=POP.CUR; avoid.inbreeding=AVOID.INB
if(Unequal_AM& CUR.GEN!=NUM.GENERATIONS){
    MATES <- assort.mate(PHENDATA=PHEN,  MATCOR=0, POPSIZE=POP.CUR, avoid.inbreeding=AVOID.INB)
}else{
    MATES <- assort.mate(PHENDATA=PHEN,  MATCOR=PHENO.MATE.CUR, POPSIZE=POP.CUR, avoid.inbreeding=AVOID.INB)
}
#cor(MATES$males.PHENDATA[,'Y'],MATES$females.PHENDATA[,'Y']) #CHECK
###################


###################
#A5 Have mates reproduce and create new genotypes and phenotypes of offspring for this generation
#Reproduce
# MATES=MATES;XO=XO;XL=XL;PHEN=PHEN;CV.INFO=CV.INFO;SUMMARY.RES=SUMMARY.RES;DELTA.PATH=delta;A.PATH=a;F.PATH=f;E.PATH=e

OFFSPRING <- reproduce(MATES=MATES,XO=XO,XL=XL,PHEN=PHEN,CV.INFO=CV.INFO,CUR.GEN=CUR.GEN,DELTA.PATH=delta,A.PATH=a,F.PATH=f,E.PATH=e)

XO <- OFFSPRING$XO
XL <- OFFSPRING$XL
PHEN <- OFFSPRING$PHEN
###################


###################
#A6 Finish loop

covs <- cov(PHEN[,c('TPO','TMO','NTPO','NTMO','AO','BV.NT.O','F','Y','YP','YM','FP','FM','TPL','TML','NTPL','NTML','AL','BV.NT.L')])
#STOPPED HERE
haps.covs <- covs[1:4,1:4]
(var.hap.obs <- mean(diag(haps.covs)))
(cov.hap.obs <- mean(haps.covs[lower.tri(haps.covs)]))
(cov.bv.bvnt.obs <- covs[5,6])
(w <- sum(covs[7,1:4])/2)
(w.parental.gen <- (sum(covs[11,c(1,3)]) + sum(covs[12,c(2,4)]))/2)
(y.parental.gen <- (sum(covs[11,c(13,15)]) + sum(covs[12,c(14,16)]))/2)
(cov.y.yp <- covs[8,9])
(cov.y.ym <- covs[8,10])
(omega <- mean(c(covs[1,9], covs[3,9], covs[2,10], covs[4,10])))
(theta.NT <- sum(covs[8,3:4]))
(theta.T <- sum(covs[8,1:2]))
(cov.yp.ym <- covs[9,10])
(var.yp <- covs[9,9])
(var.ym <- covs[10,10])

haps.covs.lat <- covs[13:16,13:16]
(var.hap.lat <- mean(diag(haps.covs.lat)))
(cov.hap.lat <- mean(haps.covs.lat[lower.tri(haps.covs.lat)]))
(var.T.lat <- var(PHEN[,'AL']))
(var.NT.lat <- var(PHEN[,'BV.NT.L']))
(cov.T.NT.lat <- covs[17,18])
(y <- sum(covs[7,13:16])/2)
(gamma <- mean(c(covs[13,9], covs[15,9], covs[14,10], covs[16,10])))
cross.covs <- covs[c(1:4,13:16),c(1:4,13:16)]
(cov.hap.obs.hap.lat <- mean(cross.covs[5:8,1:4]))


#SUMMARY.RES WITH NAMES OF THE COLUMNS NEXT TO THEM
SUMMARY.RES[CUR.GEN+1,3:ncol(SUMMARY.RES)] <- c(cor(MATES$males.PHENDATA[,"Y"],MATES$females.PHENDATA[,"Y"]), #MATE.COR
                                 nrow(PHEN), #POP.CUR
                                 var(PHEN[,'AOy']), #VAO
                                 var(PHEN[,'ALy']), #VAL
                                 var(PHEN[,'Fy']), #VF
                                 var(PHEN[,'Ey']), #VE
                                 var(PHEN[,'Y']),  #VP
                                 var(PHEN[,'AOy']+PHEN[,'ALy']) /var(PHEN[,'Y']), #h2 tot
                                 var(PHEN[,'AOy']) /var(PHEN[,'Y']), #h2 obs
                                 var(PHEN[,'ALy']) /var(PHEN[,'Y']), #h2 lat
                                 var.hap.obs, #var.hap.prs - should be 1/2 + g
                                 cov.hap.obs, #cov.hap.rs - should be g
                                 var(PHEN[,'AO']), #var.T.BV
                                 var(PHEN[,'BV.NT.O']), #var.NT.BV - should be same as above
                                 cov.bv.bvnt.obs, #cov.BV.BVNT - should be 4g
                                 w, #w - mean of 4 covs bw F & T, NT multipled by 2 - should be w
                                 w.parental.gen, #w.parental.gen - w from parental gen
                                 y.parental.gen, #y.parental.gen - y from parental gen
                                 cov.y.yp, #cov.y.yp - cov(Yo,Yp)
                                 cov.y.ym,  #cov.y.ym - cov(Yo,Ym); should be same as cov.y.yp
                                 omega, #omega - mean of 4 covs bw parental T & NT and Yp, Ym - should be 1/2 theta.T
                                 theta.NT, #theta.NT - mean of 2 covs bw Yo and NTP, NTM
                                 theta.T, #theta.T - sum of 2 covs bw Yo and TP,TM
                                 cov.yp.ym, #cov.yp.ym - spousal covariance
                                 var.yp, #var.yp - var(Yp)
                                 var.ym, #var.ym - var(Ym), should be same as VP and var.yp
                                 var.hap.lat, #new1 - 
                                 cov.hap.lat, #new1 - 
                                 var.T.lat, #new1
                                 var.NT.lat, #new1
                                 cov.T.NT.lat, #new1
                                 y, #new1
                                 gamma, #new1
                                 cov.hap.obs.hap.lat #new1
                                 )
if (SAVE.EACH.GEN){
HISTORY$MATES[[CUR.GEN+1]] <- MATES
HISTORY$PHEN[[CUR.GEN+1]] <- PHEN
HISTORY$XO[[CUR.GEN+1]] <- XO
HISTORY$XL[[CUR.GEN+1]] <- XL
}

if (SAVE.COVS){
    COVS[[CUR.GEN+1]] <- round(covs,3)}

} #End for loop
###################


###################
#A7 Return variables of interest
return(list(SUMMARY.RES=SUMMARY.RES,XO=XO,XL=XL,PHEN=PHEN,HISTORY=HISTORY,COVARIANCES=COVS))

} #End function
###################
#########################################################




##############################################################################
#END FUNCTIONS
##############################################################################



























