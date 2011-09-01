\name{peakReader}
\alias{peakReader}
\title{Retrieves the missing low ion counts from the preprocessed data for a given scan range, m/z and LC-MS dataset}
\description{
It scans the given netCDF/mzxml file with provided scan and m/z range for the retrieval of very low ion counts. The standard preprocessing steps that are applied to the LC-MS data by xcms make the data easier to analyse on a superficial level but the actual statistical distribution of preprocessed data is extremely complex and may not be suitable for rigorous statistical analysis. One step is done in direction by peakReader which calculates and returns the ions which are excluded by xcms as a part of its preprocessing algorithm because of low ion intensities so as to preserve more meaningful data. It can also plot the image graph for the data provided and the output can be sent to the \code{\link{heat}} function which plots the ion intensity image graph for displaying the particle ranges as heavy or light tails in case of the coeluting metabolities.
}
\usage{
peakReader(xraw,scan=NULL,mz,rt=NULL,plotit=FALSE,warnings=FALSE)
}
\arguments{
  \item{xraw}{
	  \code{xraw} is a \code{\link{xcmsRaw}} object which stores file information and other environment variables.
}
  \item{scan}{scan number range.}
  \item{mz}{m/z range.}
  \item{rt}{retention or acquisition time/scan time in seconds.}
  \item{plotit}{a boolean, which if is TRUE plots image graph between original ion counts returned by peakReader's processing, m/z range and the scan numbers}
	\item{warnings}{a boolean, which if is TRUE, outputs the wanings (scans variation error and empty scans for given mass) verbosely}
}
\value{
an object holding the reprocessed LCMS data as a list of
\item{ic}{original ions counts which was natively present in the LCMS data and is excluded by xcms preprocessing}
\item{mass}{m/z values of the processed compatible scans}
\item{scans}{total scans processed}
\item{mz}{m/z values of total processed scans}
\item{rt}{retention time range of total processed scans}
\item{warning1}{outputs the scan number displacement for each scan from the expected values}
\item{warning2}{outputs the scans where user inputted m/z values are not present}
}
\section{Methods}{
read$xraw=\link{xcmsRaw} file object

\code{peak=peakReader(read$xraw,scan,mz)}
}
\references{
xcms processing package: \url{http://www.bioconductor.org/packages/release/bioc/html/xcms.html}
}
\author{
Singh, Sukhdeep \email{sukhdeepsingh.bio@gmail.com}

Ipsen, Andreas \email{andreas.ipsen06@imperial.ac.uk}
}

\seealso{
\code{\link{readData}},
\code{\link{readFiles}},
\code{\link{coelutionTest}},
\code{\link{heat}}
}
