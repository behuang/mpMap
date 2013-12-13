\name{subset.binmpcross}
\alias{subset.binmpcross}
\title{Subset binned mpcross object}
\usage{
  \method{subset}{binmpcross} (x, groups = NULL,
    chr = NULL, markers = NULL, lines = NULL, ...)
}
\arguments{
  \item{x}{Object of class \code{binmpcross}}

  \item{chr}{Selected chromosomes TO KEEP. Default is all}

  \item{groups}{selected groups to KEEP. Default is all}

  \item{markers}{Selected markers TO KEEP. Default is all}

  \item{lines}{Selected lines TO KEEP. Default is all}

  \item{...}{Additional arguments}
}
\value{
  The original object with chromosomes/lines/markers
  removed which are not listed in the arguments.
}
\description{
  Reduces a binned mpcross object down to a specified set
  of chromosomes, markers and/or lines
}
\note{
  This is essentially the same function as subset.mpcross
  with the addition of dealing with the bins component.
}
\seealso{
  \code{\link[mpMap]{mpcollapse}},
  \code{\link[mpMap]{mpexpand}}
}

