readData <-
function(x) {
	# variable decalarations
	xraw<-list()
	mz<-list()
	ins<-list()
	scanidx<-list()
	scantime<-list()
	if (length(grep('^~',x,value=TRUE))!=0) {print ("Please enter complete filename, excluding shortcuts.")}
runtime=system.time((xraw=xcmsRaw(x))) ; runtime=format(runtime[3],digits=2)
mz<-xraw@env$mz ; ins<-xraw@env$intensity;scanidx<-xraw@scanindex;scantime<-xraw@scantime
		filesread=length(x)
		print (paste('It took ',runtime,'secs to read ',filesread,' file'))
	print ('Variables acquired are xraw,mz ,ins ,scanidx,scantime')
	return (list(xraw=xraw,mz=mz,ins=ins,scanidx=scanidx,scantime=scantime))
}

