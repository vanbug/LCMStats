\name{qqchisq}
\alias{qqchisq}
\title{Plots the quantile-quantile plot for determining variation between the observed and theoretical distribution.}
\description{
\code{qqchisq} processes the chi square statistics pooled data from the chromatographic scans for ploting the quantile-quantile plot. It visualizes the variation between the observed and theoretical isotopic abundance patterns for the empirical data obtained from the \code{\link{isoTest}} function.
}
\usage{
qqchisq(x,df=1,new=TRUE,main="",sub="")
}
\arguments{
  \item{x}{list of pooled chi square statistics across the chromatographic scans as returned from \code{\link{isoTest}}.}
  \item{df}{degree of freedom is N-1, where N is the numbe of lowest mass isotopologues.}
  \item{new}{default is TRUE, creates a new plot if TRUE.}
  \item{main}{plot heading of the quantile-quantile plot.}
  \item{sub}{sub heading of the quantile-quantile plot.}
}
\value{
The output is a quantile-quantile plot between the emperical and theorertical abundance pattern with a p-value.
}
\references{
1. Construction of Confidence Regions for Isotopic Abundance Patterns in LC/MS Data Sets for Rigorous Determination of Molecular Formulas
Andreas Ipsen, Elizabeth J. Want, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (17), 7319-7328.

2. A Statistically Rigorous Test for the Identification of Parent-Fragment Pairs in LC-MS Datasets
Andreas Ipsen, Elizabeth J. Want, John C. Lindon, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (5), 1766-1778.
}
\author{
Singh, Sukhdeep \email{sukhdeepsingh.bio@gmail.com}

Ipsen, Andreas \email{andreas.ipsen06@imperial.ac.uk}
}
\note{
No notes yet.
}
\seealso{
	\code{\link{coelutionTest}},
	\code{\link{isoTest}},
	\code{\link{peakReader}},
	\code{\link{qchisq}}
}
\examples{
# fetching data
data(chenodeoxycholicAcid1, chenodeoxycholicAcid2, chenodeoxycholicAcid3)
data(chenodeoxycholicAcid_Dimer1, chenodeoxycholicAcid_Dimer2, chenodeoxycholicAcid_Dimer3)

# a 'list of list' of peakReader outputs corresponding to the two clusters of isotopologues
dataList=list(list(chenodeoxycholicAcid1, chenodeoxycholicAcid2, chenodeoxycholicAcid3), list(chenodeoxycholicAcid_Dimer1, chenodeoxycholicAcid_Dimer2, chenodeoxycholicAcid_Dimer3))

#approximate isotopic abundance patterns for chenodeoxycholoc acid (first vector) and its dimer (second vector) which can obtained via generators
isoPat=list(c(0.7647, 0.2031, 0.0322), c(0.58476609, 0.31062114, 0.09049629))

Q=isoTest(dataList, isoPat, thresh=300, robust=0.1, Np=915)

#moderate p-value
Q$p.value

#chi square quantile-quantile plot constructed via \code{\link{qqchisq}} and with robustness for visualizing the variation. Degree of freedom is 2 as three lowest mass isotopologues were used. The plot signifies that observed statistics follow the theoretical isotopic abundance pattern to a large extent
qqchisq(c(Q$robust[[1]],Q$robust[[2]]),2,main="Quantile-Quantile Plot")
}
\keyword{chi square test}
\keyword{chi statistics}
\keyword{pearson correlation}
\keyword{isotope}
