\name{isoTest}
\alias{isoTest}
\title{
Plots quantile-quantile plots for determining variation between observed and theorertical isotopic abundance patterns for given metabolites.
}
\description{
Isotopic ratio test indicates whether the theorertical isotopic abundance pattern of a given chemical formula is consistent with the observed data using Pearsons chi square goodness-of-fit statistic and thereby reduces the number of candidate formulas for unknown compounds. This test determines the degree by which the goodness of fit statistics derived from the scans adhere to the theoretical distributions for the given isotopologues or other region of interests provided as dataList. Any accounted substantial departure from the 45 degree line on the resulting quantile-quantile plots would indicate the deviations from the predicted distribution concluding a contamination or non parent-fragment pairs. 
}
\usage{
isoTest(dataList,isoPat,thresh=300,robust=0.1,Np=NA)
}
\arguments{
  \item{dataList}{list of fragments as returned by the \code{\link{peakReader}} function.}
  \item{isoPat}{isotopic abundance patterns for the fragments.}
  \item{thresh}{a parameter set at 300 (default), enables to include only low ion counts preventing the detector saturation bias which effects the high ion counts}
  \item{robust}{a parameter set at 0.1 (default), enables to exclude the top 10\% of goodness of fit test statistic preventing the contamination or undesirable results.}
  \item{Np}{a parameter dependent on the detector used and is sum of pulses to make a chromatographic scan.}
}
\value{
\item{p.value}{probability value infering the goodness of fit statistic for the exact coelution of the given list of metabolites.}
\item{raw.chis}{pearson's chi square test statistics values from the data without the robustness application}
\item{robust.chis}{pearson's chi square test statistics values from the data with the robustness application }
}
\details{
Observed isotopic pattern is provided by the user which can be easily obtained by providing the m/z values , polarities, retention times or the formulae depending on the type of pattern generator used.
Sometimes, the heavy tails of adjacent isotopologues are included while centroiding the data and thus distorts the isotopic ratio patterns. A robust procedure to deal with that is applied, whose parameter value is kept 0.1 which exluded the largest 10\% of goodness of fit statistic obtained. Coates correction is applied along with thresh paramter, whose value is kept at 300 (default) which enables to include only low and moderate ion counts for processing.
The isotopic abundance pattern data observed at distinct chromatographic scans and at distinct fragments, adducts, and polymers that are derived
from the same underlying metabolite are pooled  to place further constraints on the identity. This pooling of data allows the analyst to rule out potentially much large number of candidate formulas while guarding against false positives.
}
\references{
1. Sheffield Chemputer (Isotope Patterns Calculator) : \url{http://winter.group.shef.ac.uk/chemputer/isotopes.html}

2. Isotopic Pattern Generator : \url{http://fiehnlab.ucdavis.edu/projects/Seven_Golden_Rules/Isotopic_Pattern_Generator/}

3. Isotopic Pattern Calculator : \url{http://yanjunhua.tripod.com/pattern1.htm}

4. Metabolite Search : \url{http://137.131.20.83/metabo_search_alt2.php}

5. Construction of Confidence Regions for Isotopic Abundance Patterns in LC/MS Data Sets for Rigorous Determination of Molecular Formulas
Andreas Ipsen, Elizabeth J. Want, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (17), 7319-7328

6. A Statistically Rigorous Test for the Identification of Parent-Fragment Pairs in LC-MS Datasets
Andreas Ipsen, Elizabeth J. Want, John C. Lindon, Timothy M. D. Ebbels
Analytical Chemistry 2010 82 (5), 1766-1778

7. Analytical correction for dead time effects in the measurement of time-interval distributions
Coates, P.Rev. Sci. Instrum. 1991, 63(3), 2084-2088.
}
\author{
Singh, Sukhdeep \email{sukhdeepsingh.bio@gmail.com}

Ipsen, Andreas \email{andreas.ipsen06@imperial.ac.uk}
}
\note{
No Notes yet
}
\seealso{
\code{\link{coelutionTest}},
\code{\link{peakReader}},
\code{\link{pooler}},
\code{\link{pchisq}}
}
\examples{
# loading uridine monomer and dimer fragment data which is essentially the peakReader's output
data(uridine1,uridine2,uridine_Dimer1,uridine_Dimer2)

# making a list of the fragments
dataList=list(list(uridine1,uridine2),list(uridine_Dimer1,uridine_Dimer2))

# isotopic pattern list for both structure fragments of the uridine
isoPat=list(c(0.902,0.098),c(0.822,0.178))

# Np=915 for the spectrometer which was used to run this sample. The test was run using a robustness of 0.1 as some distortion is expected due to heavy tails of mass peaks
Q=isoTest(dataList,isoPat,thresh=300,robust=0.1,Np=915)

# the pvalue observed is 0.745 which is not significant. Thus the null hypothesis can be accepted confirming that they are pure fragments of uridine
Q$p.value

# chi square quantile-quantile plot constructed with a robustness of 0.1 (default) which trims the top 10% of the statistics using \code{\link{qqchisq}} (custom quantile-quantile plot generator for the abundance patterns following chis square distribution). The chi square distribution adheres to one degree of freedom.
qqchisq(c(Q$robust[[1]],Q$robust[[2]]),1,sub=paste('Pvalue=',round(Q$p.value,4)),main="Quantile-Quantile Plot")
}
\keyword{Coelution}
\keyword{Isotope test}
\keyword{LCMS}
