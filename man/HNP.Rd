\name{HNP}
\alias{HNP}
\alias{My}
\alias{Sestan}
\docType{data}
\title{
RRHO comparison data sets.
}
\description{
RRHO comparison data sets. See references for details.
}
\usage{data(lists)}

\format{
  Three data frames: \code{HNP}, \code{My} and \code{Sestan}.
  Each is a data.frame with gene identifiers and sorting values so that they can be used as inputes to \code{\link{RRHOComparison}}.
}
    
\references{
Stein JL*, de la Torre-Ubieta L*, Tian Y, Parikshak NN, Hernandez IA, Marchetto MC, Baker DK, Lu D, Lowe JK, Wexler EM, Muotri AR, Gage FH, Kosik KS, Geschwind DH. "A quantitative framework to evaluate modeling of cortical development by neural stem cells." Manuscript in press at Neuron.  
(*) Authors contributed equally to this work.
}
\examples{
data(lists)
str(HNP) ; str(Sestan); str(My)
}
\keyword{datasets}
