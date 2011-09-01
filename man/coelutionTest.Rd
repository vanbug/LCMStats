\name{coelutionTest}
\alias{coelutionTest}
\title{
Test for identification of Parent-Fragment Pairs in LC-MS Datasets.
}
\description{
A statistically rigorous test for the identification of Parent-Fragment Pairs in the LC-MS datasets. The data of interest is provided to the \code{coelutionTest} as a pair of coeltuing ions X and Y which are essentially the lists with ion intensities, m/z range and scan ranges produced by the \link{peakReader} function. Any contamination would tend to inflate the resulting pearsons chi square statistics thus rejecting the validity of goodness of fit test and forms the basis of this test. Selections from heat plot can also be provided as an input to determine their coelution.
}
\usage{
coelutionTest(X,Y,thresh=300,robust=0.05,truth=NA,Np=NA,plotit=TRUE,pooler=1)
}
\arguments{
  \item{X}{list of ion intensities, m/z range and scan ranges produced by the \code{peakReader} for A region of interest.}
  \item{Y}{list of ion intensities, m/z range and scan ranges produced by the \code{peakReader} for B region of interest.}
  \item{thresh}{default value is 300 (excludes high ion counts greater than 300 effected by detector saturation).}
  \item{robust}{default value is 0.05 (excludes the heavy tails which might be contaminated , if the value is increased).}
  \item{truth}{default is NA and is used for pooling low ion counts by \code{pooler}.}
  \item{Np}{a parameter dependent on the detector used and is sum of pulses to make a chromatographic scan.}
  \item{plotit}{default is true and plots the scatter plots with p-values for the goodness of fit statistic of each data point, indicated by the colour code.}
  \item{pooler}{default is 1(TRUE) which pools low ion counts to  make them fit for the chi square statistics. pooler should be switched off(pooler=2) while working with centroided data as it already preserves high ion counts.
}
}
\value{
\item{p.value}{probability value which estimates the null hypothesis validity}
\item{robust.chis}{chi square statistic values for the robustness applied}
}
\details{
The test is of primary use at low and intermediate ion counts above which the detector saturation causes substantial bias to the recorded ion counts and the distributions observed do not follow Poisson distribution. Coates correction is applied to overcome the dead time induced by detector saturtaion and a variable 'thresh' is used to over come this bias as it filters out the datasets with values less than 300 (set as default) so as to be used for test. The test statistics by the \code{coelutionTest} is consistent with the null distribution. The method is made robust by excluding top 10\% of the statistics which can be a contamination from the mild mixture of adjacent isotopologues arising when the peak centroiding is done. This contamination induces a significant difference between the observed isotopic abundance pattern from the theorertical patterns and can be checked by plotting the quantile-quantile plots as done by isotopic ratio test called as \code{\link{isoTest}}. Default value of the robust parameter is kept 0.05 and for isotopologues, a value of 0.1 can used (excluding the extreme 10\% values). \code{\link{pooler}} function is used for pooling the low counts so as to make them in complaince with the chi square statistics.
}
\references{
1. Construction of Confidence Regions for Isotopic Abundance Patterns in LC/MS Data Sets for Rigorous Determination of Molecular Formulas
Andreas Ipsen, Elizabeth J. Want, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (17), 7319-7328

2. A Statistically Rigorous Test for the Identification of Parent-Fragment Pairs in LC-MS Datasets
Andreas Ipsen, Elizabeth J. Want, John C. Lindon, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (5), 1766-1778

3. Analytical correction for dead time effects in the measurement of time-interval distributions
Coates, P.Rev. Sci. Instrum. 1991, 63(3), 2084-2088.
}
\author{
Singh, Sukhdeep \email{sukhdeepsingh.bio@gmail.com}

Ipsen, Andreas \email{andreas.ipsen06@imperial.ac.uk}
}
\note{
No notes yet
}
\seealso{
  \code{\link{peakReader}},
  \code{\link{heat}},
  \code{\link{isoTest}},
  \code{\link{pooler}},
  \code{\link{coates}}
}
\examples{
data(uridine1)
#keep selection=TRUE , to manually select the regions (zooming in)
heat1=heat(uridine1,selection=FALSE)

data(uridine2)
heat2=heat(uridine2)

#switchs off the heat device to plot a new goodness of fit plot
dev.off()

# sending two selections to the coelutionTest for testing of Parent-Fragment pairs or the coelution determination. It ouputs a p-value of 0.29 which is not significant and thus accepts the null hypothesis of exact coelution
coelutionTest(heat1,heat2)

# fetching monomer and dimer data
data(uridine1,uridine2,uridine_Dimer1,uridine_Dimer2)

# test applied for checking the coelution between monomer1 and monomer2 of uridine with robust=0.05 (default for excluding small contamination buffer)
coelutionTest(uridine1,uridine2)

# test applied for checking the coelution between dimer1 and dimer2 of uridine which outputs a p-value of 0.855 which is highly non-significant and thus the dimer can be said to have exact coelution
coelutionTest(uridine_Dimer1,uridine_Dimer2)
}
\keyword{coelution}
\keyword{robustness}
\keyword{pooler}
