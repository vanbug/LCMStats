\name{LCMStats-package}
\alias{LCMStats-package}
\alias{LCMStats}
\docType{package}
\title{A package for detailed inspection and analysis of LCMS data}
\description{
A package for inspection and analysis of the LCMS data including the test for the identification of Parent-Fragment pairs, checking coelution between them and isotopic ratio test for checking variation between the theorertical and emperical formulae patterns. It includes a peak reading function which returns the low and missing ion counts excluded out by the \code{\link{xcmsRaw}} as a part of its preprocessing algorithm and a heat plot function which generates the image plot of profile data with colour matrix depending on the ion intensities. It also includes functions for coates algorithm for correcting the detector saturation bias by fixing the dead time induced at the time of high ion arrivals, for pooling the low ion counts to make it valid test data and a function for plotting quantile-quantile plots between the observed and theoretical chi square goodness of fit distributions.
}
\details{
\tabular{ll}{
Package: \tab LCMStats\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-07-15\cr
License: \tab GPL(>=3)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Sukhdeep Singh \email{sukhdeepsingh.bio@gmail.com}

Maintainer: Ipsen, Andreas \email{andreas.ipsen06@imperial.ac.uk}
}
\references{
1. Construction of Confidence Regions for Isotopic Abundance Patterns in LC/MS Data Sets for Rigorous Determination of Molecular Formulas
Andreas Ipsen, Elizabeth J. Want, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (17), 7319-7328.

2. A Statistically Rigorous Test for the Identification of Parent-Fragment Pairs in LC-MS Datasets
Andreas Ipsen, Elizabeth J. Want, John C. Lindon, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (5), 1766-1778.

3. Analytical correction for dead time effects in the measurement of time-interval distributions
Coates, P.Rev. Sci. Instrum. 1991, 63(3), 2084-2088.
}
\keyword{package}
\keyword{LCMS}
\keyword{LCMStats}
\seealso{
\code{\link{xcmsRaw-class}}
}
