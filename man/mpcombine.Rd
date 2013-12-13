\name{mpcombine}
\alias{mpcombine}
\title{Recode genotypes in informative bins}
\usage{
  mpcombine(object, thresh)
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{thresh}{The number of unique founder haplotypes
  which must be present in a bin in order to recode
  genotypes within that bin.}
}
\value{
  An object of class \code{mpcross}
}
\description{
  Recodes genotypes in bins where there are sufficiently
  diverse founder haplotypes.Markers within a bin are
  assumed to all have the same genotype. Genotypes within a
  bin are recoded based on which founder hapltotype the
  line matches. If \code{thresh} is set to the number of
  founders, this function will idetnify bins with unique
  founder haplotypes, and recode the object as if this
  cluster was a fully informative marker.If indiviudal
  lines do not match any of the bin haplotypes, genotypes
  in the bin are set to missing for that line. If genotypes
  are recoded in a bin, the bin is collapsed, and only one
  marker in that bin is retained in the final object.
}

