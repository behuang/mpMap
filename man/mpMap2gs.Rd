\name{mpMap2gs}
\alias{mpMap2gs}
\title{Export objects to Richard Mott's genome_scan}
\usage{
mpMap2gs(mpcross, cm2bp = 250000, removeMono = FALSE)
}
\arguments{
  \item{mpcross}{Object to export}

  \item{cm2bp}{Factor for converting a genetic map to
  physical map (how many bp per centiMorgan)}

  \item{removeMono}{Remove monomorphic markers - not
  required}
}
\value{
Two directories of files ready for input to reconstruction;
VARIANT.TABLES (containing founder info) and LINES
(containing progeny info). Note that the file
reconstruction_lib.c may need modifying if the number of
chromosomes simulated is nonstandard.
}
\description{
Given an mpcross object, export it in files suitable for
input to Richard Mott's reconstruction and genome_scan
(http://mus.well.ox.ac.uk/19genomes/magic.html) algorithms
}

