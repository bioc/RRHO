\name{RRHOComparison}
\alias{RRHOComparison}

\title{
Compares two RRHO maps to a third
}
\description{
Comparing two RRHO maps where one of the lists is shared between the two maps as in \{RRHO map 1: list1 vs list3\} vs \{RRHO map 2: list2 vs list3\}.
}
\usage{
RRHOComparison(list1, list2, list3, 
  stepsize, plots = FALSE, 
  labels, outputdir = NULL,
  hypermat.range)
}
%- maybe also 'usage' for other objects documented here. 
\arguments{
  \item{list1}{
A data.frame from experiment 1 with two columns, column 1 is the `Gene Identifier', column 2 is the signed ranking value (e.g. signed -log10 of p-value, or fold change).
}
  \item{list2}{
Same as \code{list1}.
}
  \item{list3}{
Same as \code{list1}.
}
  \item{stepsize}{
Integer indicating how many genes to increase by in each algorithm iteration.
}
  \item{labels}{
Character vector carrying the labels for the outputted plots.
}
  \item{plots}{
Logical. Should comparisons be plotted?
}
  \item{outputdir}{
Plot destination directory.
}
  \item{hypermat.range}{
Range of (-log) pvalues to be plotted for the BY adjusted pvalues of the difference in differences.
  }
}

\details{
The difference in \{overlap between \code{list1} and \code{list3}\} compared to the \{overlap between \code{list2} and \code{list3}\}.  This is useful for determining if there is a statistically significant difference between two RRHO maps.  In other words, this is useful for determining if the overlap between list1 and list3 is statistically different between the overlap between \code{list2} and \code{list3}.

RRHO Difference maps are produced by calculating for each pixel the normal approximation of difference in log odds ratio and standard error of overlap between the two RRHO maps. 
This Z score is then converted to a P-value and corrected for multiple comparisons across pixels [3].

The function will return a RRHO of the significance of overlap between \code{list1} and \code{list3} and \code{list2} and \code{list3}. A third RRHO gives the significance of the difference between these two overlap maps.

Note that all pvalues are outputted in -log scale.
}

\value{
A oject including:
  \item{hypermat1}{Pvalues of comparing \code{list1} to \code{list3} (-log scale).}
  \item{hypermat2}{Pvalues of comparing \code{list2} to \code{list3} (-log scale).}
  \item{Pdiff}{The pvalue of the test for a difference in difference between lists 1-3 and 2-3 (-log scale).}
  \item{Pdiff.by}{Pvalues, corrected for the search over all of the list using Benjamini-Yekutieli (-log scale).}
}
\references{
[1] Plaisier, Seema B., Richard Taschereau, Justin A. Wong, and Thomas G. Graeber. "Rank-rank Hypergeometric Overlap: Identification of Statistically Significant Overlap between Gene-Expression Signatures." Nucleic Acids Research 38, no. 17 (September 1, 2010): e169-e169.

[2] Benjamini, Y., and D. Yekutieli. "The Control of the False Discovery Rate in Multiple Testing under Dependency." ANNALS OF STATISTICS 29, no. 4 (2001): 1165-88.

[3] Stein JL*, de la Torre-Ubieta L*, Tian Y, Parikshak NN, Hernandez IA, Marchetto MC, Baker DK, Lu D, Lowe JK, Wexler EM, Muotri AR, Gage FH, Kosik KS, Geschwind DH. "A quantitative framework to evaluate modeling of cortical development by neural stem cells." Manuscript in press at Neuron.  
(*) Authors contributed equally to this work.
}
\author{
Jason Stein and Jonathan Rosenblatt
}
\seealso{
\code{\link{RRHO}}
}
\examples{
size<- 500
list1<- data.frame(GeneIdentifier=paste('gen',1:size, sep=''), 
RankingVal=-log(runif(size)))
list2<- data.frame(GeneIdentifier=paste('gen',1:size, sep=''), 
RankingVal=-log(runif(size)))
list3<- data.frame(GeneIdentifier=paste('gen',1:size, sep=''), 
RankingVal=-log(runif(size)))
(temp.dir<- tempdir())
RRHOComparison(list1,list2,list3,
                stepsize=10,
                labels=c("list1","list2","list3"),
                plots=TRUE,
                outputdir=temp.dir, 
                hypermat.range= c(-46,46)*log(10));
}