\name{mpIM}
\alias{mpIM}
\title{(Composite) Interval Mapping for QTL detection in multi-parent crosses}
\usage{
  mpIM(baseModel, object, pheno, idname = "id",
    threshold = 0.001, chr, step = 0,
    responsename = "predmn", ncov = 1000, window = 10,
    mrkpos = TRUE, ...)
}
\arguments{
  \item{baseModel}{Base phenotypic model for analysis}

  \item{object}{Object of class \code{mpcross}}

  \item{pheno}{Phenotypic object}

  \item{idname}{The idname in phenotypic data for which to
  output predicted means. Should match rownames of the
  object$finals}

  \item{threshold}{Significance threshold for QTL p-values}

  \item{chr}{Subset of chromosomes for which to compute QTL
  profile}

  \item{step}{Step size at which to compute the QTL
  profile. See \code{\link[mpMap]{mpprob}} for further
  description of default values}

  \item{responsename}{Optional input of response name to
  look for in object$pheno}

  \item{ncov}{Number of marker covariates to search for -
  default is to search for as many as possible using
  stepAIC (forward/backward selection)}

  \item{window}{Window of cM on each side of markers where
  we exclude covariates in CIM}

  \item{mrkpos}{Flag for whether to consider both marker
  positions and step positions or just steps. Is overridden
  if step=0}

  \item{...}{Additional arguments}
}
\value{
  The original input object with additional component
  QTLresults containing the following elements:
  \item{pheno}{Input phenotype data} \item{pvalue}{Each
  component contains estimated p-values at each position on
  a given chromosome} \item{wald}{Each component contains
  Wald statistics at each position on a given chromosome}
  \item{fndrfx}{Each component contains founder effects
  estimated at each position on a given chromosome}
  \item{qtl}{Each component contains the position and
  effects of a detected QTL} \item{call}{Input arguments to
  function} and with attributes describing the number of
  QTL detected, and the threshold used for detection. Will
  only return one QTL per chromosome; to find more QTL see
  \code{\link[mpMap]{findqtl2}}
}
\description{
  Interval mapping in multi-parent crosses with options for
  single-stage mixed model approach; multi-stage approach
  using predicted means; multi-stage approach including
  cofactors (CIM)
}
\details{
  Depending on the options selected, different models will
  be fit for QTL detection. If the baseModel input does not
  include a term matching the idname input, it will be
  assumed that a single-stage QTL mapping approach is
  desired. In this case, no covariates will be added (ncov
  will be set to 0); all models will be fitted in asreml;
  and all phenotypic covariates and design factors
  specified in the baseModel will be fitted along with
  genetic covariates in mixed model interval mapping.

  If the baseModel input does include a term matching the
  idname, then it will be assumed that a two-stage QTL
  mapping approach is desired. In this case, the baseModel
  will be fit using asreml and predicted means will be
  output to be used as a response in linear model interval
  mapping. If \code{ncov>0} additional marker cofactors
  will be fit; otherwise simple interval mapping will be
  run. All phenotypic covariates and design factors
  specified in the baseModel will be fit in the first
  stage.

  Note that no weights are used in the second stage of
  analysis which may result in a loss of efficiency
  compared to a one-stage approach.

  If no baseModel is input, it will be assumed that
  predicted means have been included in \code{object} as a
  phenotypic variable named predmn. In this case
  \code{pheno} is not required and asreml does not need to
  be used. (composite) Interval mapping will proceed as in
  the two-stage case depending on the value of \code{ncov}.
}
\examples{
sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
sim.ped <- sim.mpped(4, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
mpp.dat <- mpprob(sim.dat, program="qtl")
## Two-stage simple interval mapping
mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
}
\seealso{
  \code{\link[mpMap]{plot.mpqtl}},
  \code{\link[mpMap]{summary.mpqtl}}
}

