\name{readData}
\alias{readData}
\title{Reads single LC-MS (netCDF/mzxml/msn) data files}
\description{
Reads single netCDF, mzxml and msn liquid chromatography-mass spectrometry files using \code{\link{xcmsRaw}} file reader.
It stores the mass-charge(mz) ratios , corresponding signal intensities for each m/z value, scanindex (integer vector with starting positions of each scan in the mz and intensity variables) and scantimes (numeric vector with acquisition time (in seconds) for each scan). The \code{\link{xcmsRaw}} object is in scope which stores and segments the whole data of netCDF/mzxml files in various environment variables for later manipulation.
This function also displays the time taken by the reader to read the file in seconds called as runtime.
}
\usage{readData(x)}
\arguments{
  \item{x}{path of the netCDF or mzxml file}
}
\value{
a list object holding the raw segmented netCDF / mzxml data as
\item{xraw}{\code{\link{xcmsRaw}} object which stores file information and other environment variables}
\item{mz}{mass/charge values of the current file}
\item{ins}{signal intensities for the corresponding m/z values}
\item{scanidx}{scanindex (starting and ending positions of each scan)}
\item{scantime}{acquisition times for each scan}
}
\section{Methods}{
\code{read=readData('PATH_TO_CDF')}
}
\references{
  NetCDF file format:
  \url{http://www.unidata.ucar.edu/software/netcdf/}
  
  \url{http://www.astm.org/Standards/E2077.htm}
  
  \url{http://www.astm.org/Standards/E2078.htm}

  mzXML file format:
  \url{http://sashimi.sourceforge.net/software_glossolalia.html}

  PSI-MS working group who developed mzData and mzML file formats:
  \url{http://www.psidev.info/index.php?q=node/80}
  
xcms processing package: \url{http://www.bioconductor.org/packages/release/bioc/html/xcms.html}
}
\author{
Singh, Sukhdeep \email{sukhdeepsingh.bio@gmail.com}
}
\seealso{
\code{\link{xcmsRaw-class}},
\code{\link{xcmsRaw}},
\code{\link{readFiles}}
}
\keyword{read}
\keyword{file}
