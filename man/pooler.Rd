\name{pooler}
\alias{pooler}
\title{
Pools low ion counts to make data complaint with chi square statistics for detecting the Parent-Fragment pairs.
}
\description{
\code{pooler} is used to pool the low ion counts in the data provided. The count k from the Poisson distributed random data which has values less than 5, are pooled together to make it valid for the chi square test. Any contamination would tend to inflate the resulting chi square statistics which in turns reject the validity of goodness of fit test and thus deduces the coeluting nature of metabolites.
}
\usage{
pooler(counts,truth=NULL)
}
\arguments{
  \item{counts}{the matrix value consisting of ion counts for a given Poisson distributed data.}
  \item{truth}{the value is a temporary estimate of the binomial probability.}
}
\value{\code{finalCounts}a vector containing pooled low ion counts}
\references{
1. Construction of Confidence Regions for Isotopic Abundance Patterns in LC/MS Data Sets for Rigorous Determination of Molecular Formulas
Andreas Ipsen, Elizabeth J. Want, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (17), 7319-7328.

2. A Statistically Rigorous Test for the Identification of Parent-Fragment Pairs in LC-MS Datasets
Andreas Ipsen, Elizabeth J. Want, John C. Lindon, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (5), 1766-1778
}
\author{
Singh, Sukhdeep \email{sukhdeepsingh.bio@gmail.com}

Ipsen, Andreas \email{andreas.ipsen06@imperial.ac.uk}
}
\note{
No notes yet.
}

\seealso{
	\code{\link{coelutionTest}}
	\code{\link{isoTest}}
	\code{\link{peakReader}}
}
\keyword{pooler}
\keyword{coelution}
