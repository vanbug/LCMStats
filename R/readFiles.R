readFiles <-
function() {

	# fetches all files of current directory
	allFiles<-list.files(getwd())
	
	# pastes files ending with CDF (small and large caps)
	cdfData<-grep('CDF$',allFiles,value=TRUE,ignore.case=TRUE)
	mzxmlData<-grep('XML$',allFiles,value=TRUE,ignore.case=TRUE)
	
	parentDir=getwd() #stores the current directory

	# reads LC-MS files for processing
		wd<-readline('Please enter the directory name ? ')
		# string parsing, removing extra quotes
		wd=gsub('\'','',wd)
		setwd(wd)
		allFiles1<-grep('CDF$',list.files(getwd()),value=TRUE,ignore.case=TRUE)
		allFiles2<-grep('XML$',list.files(getwd()),value=TRUE,ignore.case=TRUE)
		print ('netCDF files present : ')
		print (allFiles1)
		print ('mzXML files present : ')
		print (allFiles2)
		cdfData=allFiles1
		mzxmlData=allFiles2
	
	if (length(cdfData)==0) {allData=mzxmlData}
	if (length(mzxmlData)==0) {allData=cdfData}
	if (length(cdfData)==0&&length(mzxmlData)==0) {stop("no data files present")}
	if (length(cdfData)!=0&&length(mzxmlData)!=0) {allData<-paste(cdfData,mzxmlData)}
	# variable decalarations
	xraw<-list()
	mz<-list()
	ins<-list()
	scanidx<-list()
	scantime<-list()

	# user confirmation 
	confirmA<-readline("If you wish to read all data files in current dir, (y/n) : ")
	if (confirmA=='n') {
	confirmB<-readline("If you wish to read all CDF files in current dir, (y/n) : ")
	confirmC<-readline("If you wish to read all mzxml files in current dir, (y/n) : ")
	if (confirmB=='y') {
		runtime=system.time(for (i in 1:length(cdfData)) {
			xraw[[i]]=xcmsRaw(cdfData[i])	
			mz[[i]]<-xraw[[i]]@env$mz
			ins[[i]]<-xraw[[i]]@env$intensity
			scanidx[[i]]<-xraw[[i]]@scanindex
			scantime[[i]]<-xraw[[i]]@scantime	
		})
		filesread=length(cdfData)
		}
		if (confirmC=='y') {
		runtime=system.time(for (i in 1:length(mzxmlData)) {
			xraw[[i]]=xcmsRaw(mzxmlData[i])	
			mz[[i]]<-xraw[[i]]@env$mz
			ins[[i]]<-xraw[[i]]@env$intensity
			scanidx[[i]]<-xraw[[i]]@scanindex
			scantime[[i]]<-xraw[[i]]@scantime	
		})
		filesread=length(mzxmlData)
		}
	}
	
	# loop for reading the files in, used system.time to show the process time
	if (confirmA=='y') {
		runtime=system.time(for (i in 1:length(allData)) {
			xraw[[i]]=xcmsRaw(allData[i])	
			mz[[i]]<-xraw[[i]]@env$mz
			ins[[i]]<-xraw[[i]]@env$intensity
			scanidx[[i]]<-xraw[[i]]@scanindex
			scantime[[i]]<-xraw[[i]]@scantime	
		})
		filesread=length(allData)
		}
	print ('Variables acquired are xraw,  mz , ins , scantime')
		if (confirmA=='n'&&confirmB=='n'&&confirmC=='n') {
	confirmD<-readline("Enter the first charaters of the file to read in. Enter S to read all files starting with S : ")
	if (confirmA=='n') {userFile=grep(paste('^',confirmD,sep=""),allData,value=TRUE,ignore.case=TRUE) }
	if (length(userFile)==0) {	if (confirmB=='n') {userFile=grep(paste('^',confirmD,sep=""),cdfData,value=TRUE,ignore.case=TRUE) }}
	if (length(userFile)==0) {	if (confirmC=='n') {userFile=grep(paste('^',confirmD,sep=""),mzxmlData,value=TRUE,ignore.case=TRUE) }}
		if (length(userFile)!=0) {
			runtime=system.time(for (i in 1:length(userFile)) {
				xraw[[i]]=xcmsRaw(userFile[i])
				mz[[i]]<-xraw[[i]]@env$mz
				ins[[i]]<-xraw[[i]]@env$intensity
				scanidx[[i]]<-xraw[[i]]@scanindex
				scantime[[i]]<-xraw[[i]]@scantime	
			})
			print ('Variables acquired are xraw,  mz , ins , scantime')
			filesread=length(userFile)
		} else {print ("No file is present with that starting characters"); filesread=0}
		}
	
	# formatting the runtime and filesread for display
	runtime1=format(runtime[3],digits=3)
	if (filesread==1) {(fileSen='file')} else (fileSen='files')
	print (paste('It took',runtime1,'secs','to read',filesread,fileSen ))
	runtime1=paste('Time taken by this process',runtime1,'secs')
	filesread=paste('Total number of netCDF files read:',filesread)
	
	setwd(parentDir)	#retrieves the current dir
	
	# sending results back which can be accessed later
	return (list(xraw=xraw,mz=mz,scanidx=scanidx,scantime=scantime,filesread=filesread,runtime=runtime1))
}

