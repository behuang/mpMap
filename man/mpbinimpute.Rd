\name{mpbinimpute}
\alias{mpbinimpute}
\title{Impute Missing Genotypes}
\usage{
  mpbinimpute(object, output.type = "summary")
}
\arguments{
  \item{object}{Object of class \code{mpcross}}

  \item{output.type}{What the function should return. If
  set to "summary", the function will return a 3 column
  matrix giving the location of bins where genotypes can be
  imputed, and the number of observations that can be
  imputed within the bin. The first column gives the
  chromosome, the secon gives the bin number and the third
  gives the count of imputable genotypes. Zero counts are
  omitted from the matrix. If set to "object" the function
  will return a copy of the original mpcross object with
  the missing genotypes imputed.}
}
\value{
  If output.type="summary" a 3 column matrix. If
  output.type="object" an mpcross object.
}
\description{
  Impute missing genotypes in a binned mpcross object
}

