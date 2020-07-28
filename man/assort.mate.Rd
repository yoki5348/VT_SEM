\name{assort.mate }
\alias{assort.mate }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assort.mate
%%  ~~function to do ... ~~
}
\description{
VTSEM pakcage is developed to distribute the structural equation model (SEM) for estimating the vertical transmitted environment factors via polygenic risk score (PRS) of Trio data
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SIMULATE_DAT_GEN}{
  assort.mate <- function(PHENDATA, MATCOR, POPSIZE, avoid.inbreeding=TRUE)
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yongkang Kim (yongkangkim87@gmail.com)
Jared Balbona
Matthew C Keller
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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

assort.mate(PHENDATA=PHEN,  MATCOR=PHENO.MATE.CUR, POPSIZE=POP.CUR, avoid.inbreeding=AVOID.INB)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "Assortative Mating (AM)"" }% use one of  RShowDoc("KEYWORDS")
