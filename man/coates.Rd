\name{coates}
\alias{coates}
\title{Coates algorithm for dead time correction induced by high ion arrivals in Mass Spectrometry instrument}
\description{
LC-TOF Mass Spectroscopy sytems employ a time-to-digital converter (TDC) to measure the time of flight and the count of incoming ions which is governed by the poisson distribution. Every single ion arrival triggers a period of dead time during which the TDC is incapable of registering further ions. Thus, when ion arrival rate is high , the data will produce a significant deviation from the poisson distribution which conflicts the goodness of fit statistics employed by the \code{\link{coelutionTest}} for detecting the coeluting ions.
Only a fraction of ions are counted when the rate of ion arrivals are high and those values too are biased towards the low ion counts which is referred to as detector saturation. Coates algorithm is applied for correcting the detector saturation so as to include more data which is in complaince with the poisson distribution.
}
\usage{
coates(X,Np)
}
\arguments{
  \item{X}{a vector containing ion counts.}
  \item{Np}{number of pulses to make chromatographic scan and is dependent on the type of detector system used.}
}
\value{
\item{X}{a centroided ion count data vector with applied correction algorithm.}
}
\references{
1. Analytical correction for dead time effects in the measurement of time-interval distributions
Coates, P.Rev. Sci. Instrum. 1991, 63(3), 2084-2088.

2. Construction of Confidence Regions for Isotopic Abundance Patterns in LC/MS Data Sets for Rigorous Determination of Molecular Formulas
Andreas Ipsen, Elizabeth J. Want, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (17), 7319-7328.

3. A Statistically Rigorous Test for the Identification of Parent-Fragment Pairs in LC-MS Datasets
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
\code{\link{coelutionTest}}
\code{\link{isoTest}}
\code{\link{peakReader}}
}
\examples{
# fetching data
data(chenodeoxycholicAcid1, chenodeoxycholicAcid2, chenodeoxycholicAcid3)
data(chenodeoxycholicAcid_Dimer1, chenodeoxycholicAcid_Dimer2, chenodeoxycholicAcid_Dimer3)

# a 'list of list' of peakReader outputs corresponding to the two clusters of isotopologues
dataList=list(list(chenodeoxycholicAcid1, chenodeoxycholicAcid2, chenodeoxycholicAcid3), list(chenodeoxycholicAcid_Dimer1, chenodeoxycholicAcid_Dimer2, chenodeoxycholicAcid_Dimer3))

# coates correction applied, Np value is provided dependent on the detector system used
coates(chenodeoxycholicAcid1$ic,Np=915)
coates(dataList[[2]][[1]]$ic,Np=915)
}
\keyword{coates}
