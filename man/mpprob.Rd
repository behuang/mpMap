\name{mpprob}
\alias{calcmpprob}
\alias{mpprob}
\alias{print.mpprob}
\title{Compute founder probabilities for multi-parent crosses}
\usage{
  mpprob(object, chr, step = 0, impmap, mrkpos = TRUE,
    mapfx = c("haldane", "kosambi"), ibd = FALSE,
    threshold = 0.7,
    program = c("mpMap", "qtl", "happy", "mpMIM"),
    tempfiledirectory = "", generations = 5, est = TRUE,
    geprob = 1e-04)
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{chr}{Subset of chromosomes}

  \item{step}{Step size (in cM) to create grid of positions
  at which to compute probabilities. At default value of 0,
  probabilities are calculated at marker positions only}

  \item{impmap}{Full set of positions at which to compute
  probabilities. Overrides all other parameters for
  determining positions}

  \item{threshold}{Threshold for calling founder
  probabilities}

  \item{mapfx}{Map function used to convert map to
  recombination fractions}

  \item{ibd}{Flag to indicate whether to compute
  probabilities using IBD genotypes}

  \item{mrkpos}{Flag to indicate whether to compute
  probabilities at both marker positions and step size or
  just step size. Is overridden for step size of 0. Not
  used for progam = "mpMIM"}

  \item{program}{R package to use to compute probabilities}

  \item{tempfiledirectory}{Directory in which to output
  temporary files. Default is current working directory}

  \item{generations}{Number of generations to assume in
  HAPPY. see \code{happy}}

  \item{est}{Flag for whether to impute founder alleles}

  \item{geprob}{Probability of genotyping error}
}
\value{
  The input mpcross object is returned with two additional
  components: \item{prob}{A list with founder probabilities
  for each chromosome. Format is a matrix with n.founders *
  n.markers columns and n.lines rows. Each group of
  n.founders columns will add up to 1. Founder
  probabilities are in the order of founders in the input
  founder matrix. } \item{estfnd}{A list with estimated
  founders for each chromosome. Format is a matrix with
  n.markers columns and n.lines rows. Missing values
  indicate where no founder probability exceeded the input
  threshold. Numeric values for founders indicate the row
  in the input founder matrix corresponding to the
  estimated founder.}
}
\description{
  Using haplotype probabilities, computes the probability
  that each location on a genome was inherited from each
  founder. Locations are run either at markers only, at the
  midpoints of all intervals or at step sizes of x cM.
  Probabilities can be computed using internally, or with
  R/happy.hbrem or R/qtl.
}
\details{
  If \code{program=="mpMap"} then probabilities are
  computed using flanking markers at positions across the
  genome and represent 3-point haplotype probabilities. If
  \code{program=="happy"} then probabilities are computed
  using default values in R/happy.hbrem, which calculates
  ancestral haplotypes without using pedigree information.
  This only allows for probabilities to be computed at
  midpoints of intervals. If \code{program=="qtl"} then
  probabilities are computed from multipoint founder
  probabilities in R/qtl.

  If \code{step<0} for R/mpMap or R/qtl or \code{step==0}
  for R/happy.hbrem, then probabilities are computed at the
  midpoints of marker intervals. However, if \code{step==0}
  for R/qtl or R/mpMap, probabilities are computed only at
  marker locations.
}
\examples{
sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 50, .4, 0, 0, 0), nrow=1, ncol=6, byrow=TRUE), seed=1)
mpp.dat <- mpprob(sim.dat, program="qtl")
plot(mpp.dat)
summary(mpp.dat)
}
\seealso{
  \code{\link[mpMap]{plot.mpprob}},
  \code{\link[mpMap]{summary.mpprob}}
}

