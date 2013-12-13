\name{plotGTK}
\alias{plotGTK}
\title{Visually display mapping information}
\usage{
  plotGTK(mpcross, type = "theta", groups = NULL,
    output.env = new.env())
}
\arguments{
  \item{mpcross}{Object of class mpcross to be examined
  interactively}

  \item{type}{The type of data to be displayed graphically
  - "Theta" displays recombination fractions, "lod"
  displays likelihood ratios and "combined" uses values
  that combine the two.}

  \item{groups}{The linkage groups to display graphically.
  Defaults to all.}

  \item{output.env}{An environment which is used to output
  the changes / commands made by the user. These are output
  in a variable named \code{commands} inside this
  environment. If this is omitted it is impossible to
  record the changes made interactively.}
}
\value{
  Returns \code{NULL}.
}
\description{
  Produce a heatmap of recombination fraction values, and
  allow desired changes to be specified graphically
}
\details{
  \emph{Features} This experimental function is part of a
  new workflow for grouping and ordering markers. The
  plotGTK function allows a user to: View a heatmap to
  assess the current grouping into linkage groups and
  ordering within those groups Interactively join together
  two linkage groups\cr Interactively delete specific
  markers\cr Interactively delete specified linkage groups
  Specify linkage sub-groups, so that when the object is
  reordered the sub-groups are constrained to be placed
  consecutively (although as a consequence they are
  themselves unordered)\cr Specify that a certain
  consecutive set of markers should be reordered (ignoring
  any other unselected markers)\cr

  \emph{Usage} These features are split across four modes -
  "Standard", "single-marker", "multiple-markers" and
  "reorder". Switch modes by right-clicking to bring up a
  menu.

  Mode "standard" allows the user to delete and join
  linkage groups (shortcuts <del> and <ctrl + j>
  respectively) Mode "single-marker" allows the user to
  delete single markers. Select a marker using the mouse or
  left-right arrows on the keyboard. Delete a marker using
  the <del> shortcut key. Mode "multiple-markers" allows
  the user to create linkage sub-groups. There can be up to
  ten subgroups, corresponding to the number keys 0..9. A
  range of markers is selected by highlighting a marker,
  and selecting a second marker while holding down the
  <shift> key. A range is added to group number n by
  pressing <ctrl + n>, and removed from group number n by
  pressing <ctrl + n>. Hitting <n> by itself will highlight
  all markers currently contained in group number n. Mode
  "reorder-markers" uses the same keyboard commands as mode
  "multiple-markers"

  The output of the plotGTK function is the variable
  \code{commands} in the input environment. When the
  heatmap is closed by the user, a copy of all changes and
  commands made is saved into this environment, which can
  then be `played back' by the runCommands function, which
  actually performs these operations on an \code{mpcross}
  object.
}
\examples{
map <- sim.map(len=rep(100, 2), n.mar=11, eq.spacing=TRUE, include.x=FALSE)
sim.ped <- sim.mpped(8, 1, 500, 6, 1)
sim.dat <- sim.mpcross(map=map, pedigree=sim.ped, seed=1)
#delete map entry
sim.dat$map <- NULL
sim.dat <- mpestrf(sim.dat)
grouped <- mpgroup(sim.dat, groups=3, clusterBy = "combined", method="average")
ordered <- mporder(grouped, type="2", mapfx="haldane", criterion="AR_events", use.identity=FALSE)
output.env <- new.env()
\dontrun{
#This part is interactive
plotGTK(ordered, output.env=output.env)
altered <- runCommands(ordered, get("commands", env=output.env))}
}
\seealso{
  \code{\link[mpMap]{runCommands}}
}

