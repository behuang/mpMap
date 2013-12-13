\name{fit}
\alias{fit}
\alias{fit.mpqtl}
\title{Fit a full model including all QTL and effects from base model}
\usage{
  fit(object, ...)
}
\arguments{
  \item{object}{Object of class \code{mpqtl}}

  \item{\dots}{Additional arguments to be used in
  \code{asreml}}
}
\value{
  An asreml model and summary table of QTL effects,
  p-values and Wald statistics from fitting the full model.
}
\description{
  Given the output from a scan of chromosomes with
  significant genetic variation, fits a full mixed model
  containing all effects in base model and all QTL effects.
}
\examples{
sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
fit(mpq.dat)
}
\seealso{
  \code{\link[mpMap]{mpIM}},
  \code{\link[mpMap]{summary.mpqtl}}
}

