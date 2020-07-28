\name{perform_SEM_model2_eq}
\alias{perform_SEM_model2_eq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run VT-SEM model 2eq
%%  ~~function to do ... ~~
}
\description{
Run VT-SEM model 2eq which assumes there is sequilibrium AM and latent polygenic genetic score(LGS)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  Data frame which includes offsprings' phenotype value (Yo), parental phenotype values (Ym,Yp), haplotypic transmitted PGS from parents (Tm,Tp), and haplotypic non-transmitted PGS of parents (NTm,NTp)
  }
%%     ~~Describe \code{x} here~~
  \item{dat,NTmlabel,Tmlabel,NTplabel,Tplabel,Ymlabel,Yplabel,Yolabel}{
  Column names of NTm,Tm,NTp,Tp,Ym,Yp, and Yo.
  }
  \item{max.cores}{
  Assign the number of cores to run OpenMx.
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{AFE.Fit}{
   OpenMx objective which includes fitted results
  }
  \item{result_summary}{
  Summaries of fitted values of model 0
  }
}
\references{
Kong, Augustine, et al. "The nature of nurture: Effects of parental genotypes." Science 359.6374 (2018): 424-428.
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
data(equil_exam)
model2eq=perform_SEM_model2eq(equil_exam,"NTMO","TMO","NTPO","TPO","YM","YP","Y")
summary(model2eq$AFE.Fit)
model2eq$result_summary
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "VT-SEM model 2eq" }% use one of  RShowDoc("KEYWORDS")
