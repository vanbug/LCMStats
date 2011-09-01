\name{heat}
\alias{heat}
\title{
Generates the image plot of the processed data.}
\description{Heat generates the image plot of the continuum and centroid data. It is useful to show the parent fragment pairs of metabolities as a cluster of ions , one m/z distance apart called as isotopologues. Image selection is devised which allows the user to select a cluster of interest by a single line diagonal selection and then right click switches off the selection. After the selection is done, the selected plot area and the original heat plot opens visualizing the selected area. It essentially enables the user to zoom and record a particular region of interest in the heat plot which is very useful tool to select the significant ion cluster while excluding the long heavy tail which might cause distortions to the chi square plots and a significant devaition from true p-value is also observed. For a particular peak data, the scan and m/z range need not to be specified as it is directy fetched from peak data outputted by the \code{\link{peakReader}} function.}
\usage{
heat(Q,scanRange=NULL,massRange=NULL,scanRange1=c(1,nrow(Q[[1]])),xlim=NULL,ylim=NULL,selection=FALSE)
}
\arguments{
  \item{Q}{
  	\code{Q} is the list of ioncounts, m/z range, scan range and the retention times which can be used to generate heat image plot. The output of \code{\link{peakReader}} function is fed as parameter Q to the heat plot which has the filled ion intensities for very low ion values(~0) known as original ion counts along with the total m/z range and scan range or retention time as supplied.
}
  \item{scanRange}{scan number range.}
  \item{massRange}{mass/charge range.}
  \item{scanRange1}{ is a vector with ion count's total column number for heat code's internal processing}
  \item{xlim}{the x limits (x1,x2) of the plot.}
  \item{ylim}{the y limits (y1,y2) of the plot.}
  \item{selection}{a boolean(default is FALSE), if TRUE, enables the manual single line diagonal selection of the region of interest from the heat plot.}
}
\value{
Heat plot returns a list of values as:
\item{marks}{selection coordinates for the selection as finalized by the user.}
\item{ic}{ion counts of the selection.}
\item{mz}{m/z values of the selection.}
\item{scans}{scans of the selection.}
If no selection is made, the original input data is returned as the output.
}
\author{
Singh, Sukhdeep \email{sukhdeepsingh.bio@gmail.com},

Ipsen, Andreas \email{andreas.ipsen06@imperial.ac.uk}
}
\note{
No notes yet
}
\seealso{
\code{\link{peakReader}}
\code{\link{image}}
\code{\link{coelutionTest}}
}
\examples{
data(chenodeoxycholicAcid1)
#keep selection=TRUE to manually select the regions (zooming in)
heat1=heat(chenodeoxycholicAcid1,selection=FALSE)

data(chenodeoxycholicAcid2)
heat2=heat(chenodeoxycholicAcid2)

#switchs off the heat device to plot a new goodness of fit plot
dev.off()
# sending two selections to the coelutionTest for testing of Parent-Fragment pairs or the coelution determination.
coelutionTest(heat1,heat2)
}
\keyword{image}
