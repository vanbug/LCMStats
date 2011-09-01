\name{readFiles}
\alias{readFiles}
\title{Reads multiple LCMS files}
\description{
Reads and stores multiple netCDF / mzxml files as a list. It is capable of reading all netCDF and mzxml files present in the specified directory.
It stores the mass-charge(mz) ratios , corresponding signal intensities for each m/z value, scanindex (integer vector with starting positions of each scan in the mz and intensity variables) and scantimes (numeric vector with acquisition time (in seconds) for each scan) of the inputted files. The \code{\link{xcmsRaw}} object is in scope which stores and segments the whole data of netCDF/mzxml files in various environment variables for later manipulation. It also displays the runtime after the process completion.
}
\usage{
readFiles()
}
\value{
a list with multiple objects holding the segmented raw netCDF / mzxml data as
\item{xraw}{\code{\link{xcmsRaw}} object which stores file information and other environment variables}
\item{mz}{mass/charge values of the files}
\item{ins}{signal intensities for the corresponding m/z values of the files}
\item{scanidx}{scanindex (starting and ending positions for each scan) of the files}
\item{scantime}{acquisition times for each scan of the files}
}
\section{Methods}{
\code{read=readFiles()}
}
\references{
	NetCDF file format: 
		\url{http://www.unidata.ucar.edu/software/netcdf/}
		
		\url{http://www.astm.org/Standards/E2077.htm}
		
		\url{http://www.astm.org/Standards/E2078.htm}
		
	mzXML file format: 
		\url{http://sashimi.sourceforge.net/software_glossolalia.html }
		
	xcms processing package: 
		\url{http://www.bioconductor.org/packages/release/bioc/html/xcms.html}
}
\author{
Singh, Sukhdeep \email{sukhdeepsingh.bio@gmail.com}
}
\section{Warning}{
Use with caution for reading directory, as owing to the complex nature and size of the data files , the system may respond very slowly untill the process is finished. A large amount of random access memory will be consumed for the process to complete depending on the data files and the system architecture.
}
\seealso{
\code{\link{xcmsRaw-class}},
\code{\link{xcmsRaw}},
\code{\link{readData}}
}
\keyword{read}
\keyword{file}
\keyword{LCMStats}
\keyword{lcms}
