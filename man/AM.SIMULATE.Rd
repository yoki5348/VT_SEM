\name{AM.SIMULATE}
\alias{AM.SIMULATE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AM.SIMULATE
%%  ~~function to do ... ~~
}
\description{
Simulate trio dataset for the analysis of VT-SEM
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CV.INFO}{
  Information of causal variants, must be 2 columned data frame.
  Column 1: Minor allele frequencies (MAF) of causal variants, Column 2: Effect sizes of causal variants
  }
%%     ~~Describe \code{x} here~~
  \item{H2.TO}{
  Initial heritability
  }
  \item{NUM.GENERATIONS}{
  Number of generations
  }
  \item{POP.SIZE}{
  Population sizes
  }
  \item{MATE.COR}{
  Assortative mating correlation : If there is no AM, MATE.COR=0
  }
  \item{AVOID.INB}{
  Bullean which indicates whether avoiding inbreeding or not (default if FALSE)
  }
  \item{SAVE.EACH.GEN}{
  Bullean which indicates whether saving whole dataset from initial generation to current generation or only current generation.
  }
  \item{SAVE.COVS}{
  Bulean which indicates saving covariances of random variables or not.
  }
  \item{SEED}{
  Set seed value.
  }
  \item{VF.T0,PROP.H2.LATENT}{
  Initial variance induced by vertical transmission
  }
  \item{Unequal_AM}{
  Bulean which indicates using disequilibrium AM assumption or not.
  }


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{SUMMARY.RES}{
   Summary the simulation results from first generation to last.
  }
  \item{XO}{
  Generated genotype dataset of observed polygenic genetic scores (Additive coding)
  }
  \item{XL}{
  Generated genotype dataset of latent polygenic genetic scores (Additive coding)
  }
  \item{PHEN}{
  Generated phenotype dataset
  }
  \item{HISTORY}{
  Generated dataset which includes whole generation
  }
  \item{COVARIANCES}{
  Covariance matrices of phenotype dataset for each generation
  }
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Tahmasbi, Rasool, and Matthew C. Keller. "GeneEvolve: a fast and memory efficient forward-time simulator of realistic whole-genome sequence and SNP data." Bioinformatics 33.2 (2017): 294-296.
}
\author{
Yongkang Kim \email{yongkangkim87@gmail.com}, Jared Balbona, and Matthew C Keller
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

AM.DATA <- AM.SIMULATE(CV.INFO=CV.INFO, H2.T0=herit0, NUM.GENERATIONS=num.gen, POP.SIZE=pop.size*3, MATE.COR=am, AVOID.INB=avoid.inb, SAVE.EACH.GEN=save.history, SAVE.COVS=save.covariances, SEED=seed, VF.T0=vf0,PROP.H2.LATENT=prop.h2.latent,Unequal_AM=Unequal_AM) #Dataset generation
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "Assortative Mating (AM)" }% use one of  RShowDoc("KEYWORDS")
