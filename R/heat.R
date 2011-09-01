heat <-
function(Q,scanRange=NULL,massRange=NULL,scanRange1=c(1,nrow(Q[[1]])),xlim=NULL,ylim=NULL,selection=FALSE){
	
	# fetching scanRange and massRange from the peakReader output as Q.
	if (is.null(scanRange)==TRUE) {	scanRange<-c() ; 	scanRange[1]=Q[[3]][[1]] ; 	scanRange[2]=Q[[3]][[length(Q[[3]])]]} else {scanRange=scanRange}
	if (is.null(massRange)==TRUE) {massRange=Q[[4]]} else {massRange=massRange}
	Q[[1]]=round(Q[[1]])
	# stores the maximum Q[[1]] value which is Inferred ion count
	exoMax=max(Q[[1]][Q[[1]]!=Inf])

	# creates a range of exoSeq with 0.5 difference
	exoSeq=c(-0.5:(exoMax+0.5))

	# variable decalarations
	colourPower=2; selPoints<-list()

	# colscheme stores a series of colours in html colour code format using rgb which takes 3 numbers as input and a colour code is output
	colscheme =c("black",rgb(cbind((c(1: exoMax)/exoMax)^ colourPower,(4*c(1: exoMax)/exoMax*(1-c(1: exoMax)/exoMax))^1, (1-c(1: exoMax)/exoMax)^ colourPower+(c(1: exoMax)/exoMax)^(50*colourPower))))

	# function for heat plot generation and for locating points
	image.identify <- function(x,y,z, mark=TRUE, digits=3,switcher,zoom,colBar,selection){
	
if (colBar==TRUE) {
	# layout function used which devides the device space into the row and column width specified , plot division
	layout(matrix(c(1,2),nrow=1),c(0.20,1))

	# par function is used to send graphical parameters to the device such as cex : used for magnifying font
	par(cex=1)

	## mar controls the number of lines of margin to be specified on the four sides of the plot
	par(mar=c(3.5,3.5,2,1.55)+0.1) #a=width of bar, b=keeps it left,d=keeps it right , their ratio should be balanced (b+d(3.5+1.5)=b+d(3+2))

	plot(c(0,1),c(0,max(Q[[1]])),t="n",lab=c(1,5,7),xlab="",yaxt="n",ylab="",xaxt="n",yaxs="i",xaxs="i") #makes the bar boundry

	axis(2,las=1)
	par(srt=180)
	text(2.45,max(Q[[1]])/2,"ion count",srt=270,xpd=TRUE)
	if(exoMax<=10){
	for(i in 1:(length(colscheme))){
	
		loopSeq=seq(exoSeq[i],exoSeq[i+1],length.out=50)
		for(j in 1: 50){
			lines(c(0,1),rep(loopSeq[j],2), col=colscheme[i],lwd=1)
			}
		}	
		}
	if((exoMax<50)&(exoMax>10)){
	for(i in 1:(length(colscheme))){
	
		loopSeq=seq(exoSeq[i],exoSeq[i+1],length.out=25)
		for(j in 1: 25){
			lines(c(0,1),rep(loopSeq[j],2), col=colscheme[i],lwd=1)
			}
		#lines(c(-0,1),rep(exoSeq[i+1],2),col=colscheme[i],lwd=1)
		}	
		}
	if(exoMax>=50){
	for(i in 1:(length(colscheme))){
	
		loopSeq=seq(exoSeq[i],exoSeq[i+1],length.out=10)
		for(j in 1: 10){
			lines(c(0,1),rep(loopSeq[j],2), col=colscheme[i],lwd=1)
			}
		#lines(c(-0,1),rep(exoSeq[i+1],2),col=colscheme[i],lwd=1)
		}	
		}
	}

		# switch statements for creating image depending on the xlim and ylim occupancy 		
		if (switcher=='1') {image(z=Q[[1]],x=c(scanRange[1]:scanRange[2]),y= Q[[2]],ylab="m/z",xlab="scan number",col=colscheme,xlim=xlim,ylim=ylim)}
		if (switcher=='2') {image(z=Q[[1]],x=c(scanRange[1]:scanRange[2]),y= Q[[2]],ylab="m/z",xlab="scan number",col=colscheme,ylim=ylim)}
		if (switcher=='3') {image(z=Q[[1]],x=c(scanRange[1]:scanRange[2]),y= Q[[2]],ylab="m/z",xlab="scan number",col=colscheme,xlim=xlim)}
		if (switcher=='4') {image(z=Q[[1]],x=c(scanRange[1]:scanRange[2]),y= Q[[2]],ylab="m/z",main="Heat Plot",xlab=paste("scan number\nScan Range = ", scanRange[1],'-',scanRange[2]),col=colscheme)}

		if (selection==TRUE) {
		print ("Click the graph to pick Coordinates or press right click to exit or finish selection mode")
		res <- data.frame() 
		if (zoom==TRUE) {
		
		# locator function for locating clicked points on graph
		xyL <- locator(2,type="l",col="white") 
			
			# error handling : if clicked outside the plot
		while(!is.null(xyL)) {if (xyL$y[1]<massRange[1]||xyL$y[1]>massRange[2]||xyL$y[2]<massRange[1]||xyL$y[2]>massRange[2]||xyL$x[1]<scanRange[1]||xyL$x[1]>scanRange[2]||xyL$x[1]<scanRange[1]||xyL$x[1]>scanRange[2]) {print ("Warning: You clicked outside the heat plot")}
		
			xbin <- as.numeric(cut(xyL$x,x)) 
			ybin <- as.numeric(cut(xyL$y,y)) 
			
			# signal checker, 0 if click is inside the plot , 1 if its outside and skips the rest code
			GO=0 ; if(is.na(xbin[1])||is.na(xbin[2])||is.na(ybin[1])||is.na(ybin[2])) {GO=1}
		
			if (GO!=1) {
				if(mark) { 
					points(xyL$x,xyL$y,pch=10,cex=.50,col="grey")																# points properties
					#text(xyL$x,xyL$y,format(z[xbin,ybin],digits=digits),adj=-.2,cex=0.8,col="white")	# text properties controlling point text
				} 
				# output display properties
				#cat("Clicked Cordinate : ","[",xyL$x,",",xyL$y,"] ",", Ion Count=",z[xbin,ybin],"\n",sep='') 
				cat("[",xbin,",",ybin,"] = ",z[xbin,ybin],"\n",sep='') 
				res <- rbind(res,data.frame(i=xbin,j=ybin,x=xyL$x,y=xyL$y,z=z[xbin,ybin]))
				xyL <- locator(2,type='l',col="white")
			}
			if (GO==1){xyL <- locator(2,type='l',col="white")}
		}

		# image zooming : one selection is zoomed
		if (length(res$x)==2){
			# variable decalarations
			scn1<-c() ; scn2<-c() ; ic1<-c() ; ic2<-c() ; massbin1<-c() ; massbin2<-c() ; 
			boxicRangetemp<-list() ; boxicRange<-c() ; 	boxScanRange <-c(); boxMassRange <-c()
			ic1=res$i[1] ; ic2=res$i[2]
			scn1=res$x[1] ; scn2=res$x[2]
			massbin1=res$j[1] ; massbin2=res$j[2]
			if (res$j[1]>res$j[2]) {massbin1=res$j[2];massbin2=res$j[1]} #enables the upper back and forward selection
			boxMassRange=Q[[2]][massbin1:massbin2]
			if (round(res$x[1])>round(res$x[2])) {boxScanRange=c(round(scn2):round(scn1)) } #backslash selection
			if (round(res$x[1])<round(res$x[2])) {boxScanRange=c(round(scn1):round(scn2))}	#forwardslash selection
			for (i in massbin1:massbin2) {
				if (ic2>ic1) {
					if (length(ic1:ic2)>length(boxScanRange)) {ic2=ic2-(length(ic1:ic2)-length(boxScanRange))}
				boxicRangetemp[[i-massbin1+1]]=Q[[1]][,i][ic1:ic2]}
				if (ic1>ic2) {
					if (length(ic2:ic1)>length(boxScanRange)) {ic1=ic1-(length(ic2:ic1)-length(boxScanRange))}
				boxicRangetemp[[i-massbin1+1]]=Q[[1]][,i][ic2:ic1]}
			}
			for (i in 1:((massbin2+1)-massbin1)) {boxicRange=cbind(boxicRange,boxicRangetemp[[i]])}
		}
		
		if(length(res$x)>2) {stop("Please make a single selection")}
			#error handling : if selection is not made
			if (length(res$x)==0) {stop('You made an incomplete selection! ')}

		return (list(res1=res,boxScanRange=boxScanRange,boxMassRange=boxMassRange,selectedicRange=boxicRange))
			}
		}
	} #image.identify function finishes
	
	# calls the image.identify
	imageit=function(selection) {
		if((!is.null(xlim))&(!is.null(ylim))){image1=image.identify(x=c(scanRange[1]:scanRange[2]),y= Q[[2]],z=Q[[1]],switcher=1,zoom=TRUE,colBar=TRUE,selection=selection)}
		if(is.null(xlim)&(!is.null(ylim))){image1=image.identify(x=c(scanRange[1]:scanRange[2]),y= Q[[2]],z=Q[[1]],switcher=2,zoom=TRUE,colBar=TRUE,selection=selection)}
		if((!is.null(xlim))&(is.null(ylim))){image1=image.identify(x=c(scanRange[1]:scanRange[2]),y= Q[[2]],z=Q[[1]],switcher=3,zoom=TRUE,,colBar=TRUE,selection=selection)}
		if(is.null(xlim)&is.null(ylim)){image1=image.identify(x=c(scanRange[1]:scanRange[2]),y= Q[[2]],z=Q[[1]],switcher=4,zoom=TRUE,,colBar=TRUE,selection=selection)}
	
	if(selection==FALSE){return(list(ic=Q[[1]],mz=Q[[2]],scans=c(scanRange[1]:scanRange[2])))}
	if(selection==TRUE){
		# new device opened for plotting	
		dev.new()
		par(mfrow=c(1,3))
	
		#plots the zoomed image
		image(image1$boxScanRange,image1$boxMassRange,image1$selectedicRange,ylab="m/z",main="Heat Plot",xlab=paste("scan number\nScan Range = ", scanRange[1],'-',scanRange[2]),col=colscheme)
		image.identify(x=c(scanRange[1]:scanRange[2]),y= Q[[2]],z=Q[[1]],switcher=4,zoom=FALSE,colBar=FALSE,selection=TRUE)
		abline(h=image1$boxMassRange[1],v=image1$boxScanRange[1],col="white",lty=3)
		abline(h=tail(image1$boxMassRange,1),v=tail(image1$boxScanRange,1),col="white",lty=3)
	
		# vertical colour/IC bar code
		# par function is applied to make space for the colour bar
		par(mar=c(5,0.4,4,16)+0.1)

		# plot function then creates the bar vertically on the side heat graph
		plot(c(0,1),c(0,max(Q[[1]])),t="n",lab=c(1,5,7),xlab="",yaxt="n",ylab="",xaxt="n",yaxs="i",xaxs="i")

		# displayes the colour bar value range
		axis(4,las=1)
		par(srt=180)

		# text is used to display 'ioncount' text in the mid of vertical colour bar
		text(2.45,max(Q[[1]])/2,"ion count",srt=270,xpd=TRUE)

		# loops for drawing the vertical ion count line, range depending upon the conditions
		## 1
		if(exoMax<=10){
			for(i in 1:(length(colscheme))){
				loopSeq=seq(exoSeq[i],exoSeq[i+1],length.out=50)
				for(j in 1: 50){
					# lines join the x,y coordinates by drawing a line between them
					lines(c(0,1),rep(loopSeq[j],2), col=colscheme[i],lwd=1)
				}
			}	
		}
		## 2
		if((exoMax<50)&(exoMax>10)){
			for(i in 1:(length(colscheme))){
				loopSeq=seq(exoSeq[i],exoSeq[i+1],length.out=25)
				for(j in 1: 25){
					lines(c(0,1),rep(loopSeq[j],2), col=colscheme[i],lwd=1)
				}
			}	
		}
		## 3
		if(exoMax>=50){
			for(i in 1:(length(colscheme))){
				loopSeq=seq(exoSeq[i],exoSeq[i+1],length.out=10)
				for(j in 1: 10){
					lines(c(0,1),rep(loopSeq[j],2), col=colscheme[i],lwd=1)
				}
			}	
		}
	
		layout(matrix(c(1),nrow=1),1)
	
		# error handling : is user happy with selection!!
		ansImage=readline('Are you happy with your selection (y/n)?')

		# creates plot again or returns the clicked co-ordinates and zoomed ics
		if (ansImage=='n') {
			if (length(dev.list())==1) {dev.off(); imageit(selection)}
			if (length(dev.list())==2) {dev.off(); dev.off(); imageit(selection)}
			if (length(dev.list())==0) {imageit(selection)}
			if (length(dev.list())>2){imageit(selection)}
			return(list(marks=image1$res1,ic=image1$selectedicRange,mz=image1$boxMassRange,scans=image1$boxScanRange))
			} else {return(list(marks=image1$res1,ic=image1$selectedicRange,mz=image1$boxMassRange,scans=image1$boxScanRange))}
		}
	}
	# creates the initial heat plot 
	imageit(selection)
}

