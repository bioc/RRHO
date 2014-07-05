\name{RRHO}
\alias{RRHO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rank-Rank Hypergeometric Overlap Test}

\description{The function tests for significant overlap between two sorted lists using the method in the reference.}

\usage{
  RRHO(
  list1, list2, 
  stepsize = defaultStepSize(list1, list2), 
  labels, 
  alternative,
  plots = FALSE, 
  outputdir = NULL, 
  BY = FALSE,
  log10.ind=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list1}{
   data.frame. 
   First column is the element (possibly gene) identifier, and the second is its value on which to sort.  For differential gene expression, values are often -log10(P-value) * sign(effect).
}
  \item{list2}{
   data.frame. Same as \code{list1}.
}
  \item{stepsize}{
	  Controls the resolution of the test: how many items between any two overlap tests. 
}
  \item{labels}{
	  Character vector with two elements: the labels of the two lists.
}
\item{alternative}{
  Either "enrichment" for a one sided test, or "two.sided" for a two sided test. See Details section.
}
  \item{plots}{
	  Logical. Should output plots be returned?
}
  \item{outputdir}{
	  Path name where plots ae returned.
}

  \item{BY}{
	  Logical. Should Benjamini-Yekutieli FDR corrected pvalues be computed?
}
\item{log10.ind}{
    Logical. Should pvalues be reported and plotted in -log10 scale and not -log scale?
}
}


\details{
	Following the method in the reference, the function computes the number of overlapping elements in the first \eqn{i*stepsize} and \eqn{j*stepsize} elements of each list, and return the observed significance of this overlap using a hypergeometric test (see \code{\link{fisher.test}}).
	The output is returned as a list of matrices including: the overlap in the first \eqn{i*stepsize,j*stepsize} elements and the significance of this overlap. 
  
  If \code{plots=TRUE} then plots of these matrices are stored in .jpg format. 
  In the case of \code{alternative='two.sided'} the pvalue plots are signed, just like in [1], thus distinguishing between over and under enrichment.}

\value{
\item{hypermat}{Matrix of \eqn{-log(pvals)} of the test for the first \eqn{i,j} elements of the lists.}
\item{hypermat.counts}{Counts of the number of agreements in the first \eqn{i,j} elements of the lists.}
\item{hypermat.by}{An optional output of the B-Y corrected p-values of \code{hypermat} }
\item{hypermat.signs}{Matrix of the type of deviation from the null. Negative for underenrichment and positive for overenrichment.}
}

\section{Notes}{
  By default, pvalues are reported in (minus) the natural log scale and not in (minus) log 10 scale. 
  This behaviour is governed by \code{log10.ind}.
  
  The p-values of the two-sided hypothesis test differ from those in reference [1]. This is because the two-sided p-values suggested in [1], are based on taking either the upper or lower tail of the distribution without appropriately using both tails.  This method does not correctly control the type I error rate.  In the implementation here, for a two-sided test we sum the probabilities from both tails of the hypergeometric distribution.
  See the package vignette for a small simulation.
}



\references{
[1] Plaisier, Seema B., Richard Taschereau, Justin A. Wong, and Thomas G. Graeber. "Rank-rank Hypergeometric Overlap: Identification of Statistically Significant Overlap Between Gene-expression Signatures." Nucleic Acids Research 38, no. 17(September 1, 2010)

[2] Benjamini, Y., and D. Yekutieli. 2001. "The Control of the False Discovery Rate in Multiple Testing Under Dependency." ANNALS OF STATISTICS 29 (4): 1165-1188.

[3] Stein JL(*), de la Torre-Ubieta L(*), Tian Y, Parikshak NN, Hernandez IA, Marchetto MC, Baker DK, Lu D, Lowe JK, Wexler EM, Muotri AR, Gage FH, Kosik KS, Geschwind DH. "A quantitative framework to evaluate modeling of cortical development by neural stem cells." Manuscript in press at Neuron.  
(*) Authors contributed equally to this work.
}
\author{
	Jonathan Rosenblatt and Jason Stein
}


\seealso{
  \code{\link{pvalRRHO}}; \code{\link{RRHOComparison}}
  }

\examples{
	list.length <- 100
	list.names <- paste('Gene',1:list.length, sep='')
	gene.list1<- data.frame(list.names, sample(100))
	gene.list2<- data.frame(list.names, sample(100))
  # Enrichment alternative
	RRHO.example <-  RRHO(gene.list1, gene.list2, alternative='enrichment')
	image(RRHO.example$hypermat)
  
  # Two sided alternative
  RRHO.example <-  RRHO(gene.list1, gene.list2, alternative='two.sided')
	image(RRHO.example$hypermat)

}
\keyword{htest}
