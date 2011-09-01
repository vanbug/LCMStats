peakReader <-
function(xraw,scan=NULL, mz,rt=NULL,plotit=FALSE,warnings=FALSE) {
scanRange=scan; massRange=mz
# variable decalarations
fullMass<-c(); fullIc<-c(); valdex=c(); voidTracker=c(); voidIndex=c(); myIC=list(); myMass=list(); myScans_index<-list() ; myScans<-list() ; myRT=c(); inferredMass=c(); inferredIC=c() ; minLength<-c();myMass<-list() ; warn1<-c() ; warn2<-c() ; scan_list<-c() ; scans_in<-c()

	# fetching data
	mz<-xraw@env$mz
	ins<-xraw@env$intensity
	scanidx<-xraw@scanindex
	scantime<-xraw@scantime
	totalScanlength<-length(scanidx)

# calculates scanRange corresponding to retention/scan/ acquisition time sent
	# if scanRange is not provided
	if (length(scanRange)==0) {for(i in 1:length(scantime)) {if (scantime[i]<=rt[1]&&rt[1]<=scantime[i+1]) {scanRange[1]=i}
	if (scantime[i]<=rt[2]&&rt[2]<=scantime[i+1]) {scanRange[2]=i}}}

	# if scanRange is provided
	if (length(scanRange)!=0) {rt=c(scantime[scanRange[1]]:scantime[scanRange[2]])}


	# warning decalarations
	if((massRange[2]-massRange[1])>1){print("Warning: very large mass range")}
	if(scanRange[2]>totalScanlength){print("Input scanRange is out of bounds for this data")}

	# main loop driving the function structure
	for(i in scanRange[1]:scanRange[2]) {
		# stores m/z ratio of the corresponding scanRange from xcmsraw
		fullMass=mz[c(scanidx[i]:scanidx[i+1])]
		fullMass=fullMass[-1]

		# stores the ion intensity/ion counts of the corresponding scanRange from xcmsraw
		fullIc=ins[c(scanidx[i]:scanidx[i+1])]
		fullIc=fullIc[-1]

		# boolean filtering for true masses in range
		valdex=(fullMass>=massRange[1])&(fullMass<=massRange[2])

		# code division in 2 cases : 
		# if filtered IC is non zero
		if((sum(fullIc[valdex])>0)){
			# storing masses in range
			mass= fullMass[valdex]

			# storing ion count in range
			ic= fullIc[valdex]
			
			# storing the scan numbers
			#scans=scanidx[valdex]
			# storing the scan indices
			#scans_index=which(valdex)
			
			if (length(mass)<=1) (stop("The input data is not compatible with the code, please widen your mz/scan range selection or enter some other data"))
			# variability check in adjacent m/z values
			massDiffs=diff(mass)
	
			# checking nearly similar values to identify the cluster of ticks
			trueTicks=(round(massDiffs/min(massDiffs))==1)

			# checking for discrete values
			nonTicks=c(1:length(trueTicks))[trueTicks!=1]

			# specifying range
			nonTicks=c(0,nonTicks,length(mass))

			# range of true clusters
			tickRange=c(min(massDiffs[trueTicks]),max(massDiffs[trueTicks]))
	
			# median can be used for the robustness
			iniTick=mean(massDiffs[trueTicks])

			# inferredTicks storing discrete mass variation
			inferredTicks=massDiffs[!trueTicks]/iniTick

				# calulates the variation in observed scan number from expected
				if(max(massDiffs/iniTick-round(massDiffs/iniTick))>0.1){
				if (warnings==FALSE) {warn1=rbind(warn1,(paste("warning: scan ",i," has an error of ", max(massDiffs/iniTick-round(massDiffs/iniTick)))))}
				else {warn1=rbind(warn1,print(paste("warning: scan ",i," has an error of ", max(massDiffs/iniTick-round(massDiffs/iniTick)))))}
				}
				# variable decalarations
				inferredMass=c(); inferredIC=c()

				# calculates inferredMass and inferredIC				
				if(sum(trueTicks)!=length(trueTicks)){
					for(v in 1:length(inferredTicks)){
						inferredMass=c(inferredMass , mass[(nonTicks[v]+1):(nonTicks[v+1])], seq(mass[nonTicks[v+1]]+iniTick,mass[nonTicks[v+1]+1]-iniTick, length.out=round(inferredTicks[v]-1)))
						inferredIC =c(inferredIC, ic[(nonTicks[v]+1):(nonTicks[v+1])], seq(0,0,length.out=round(inferredTicks[v]-1)))	
					}
				}
				
			inferredMass=c(inferredMass, mass[c((nonTicks[length(nonTicks)-1]+1):(nonTicks[length(nonTicks)]))])
			inferredIC =c(inferredIC, ic[c((nonTicks[length(nonTicks)-1]+1):(nonTicks[length(nonTicks)]))])
			tick= iniTick
			nStartTicks=floor((min(inferredMass)-massRange[1])/tick)
			nEndTicks=floor((massRange[2]-max(inferredMass))/tick)

			# totalMass
			totalMass=c(seq(min(inferredMass)-nStartTicks*tick, min(inferredMass)-1*tick, length.out=nStartTicks),	inferredMass,
			seq(max(inferredMass)+1*tick, max(inferredMass)+ nEndTicks*tick, length.out= nEndTicks))
				
			# shiftedMass
			shiftedMass=c(totalMass, totalMass[length(totalMass)]+tick)-tick/2

			# shiftedIC
			shiftedIC=c(rep(0, nStartTicks),inferredIC,rep(0, nEndTicks))

			# final mass and IC calculations
			myMass[[i-scanRange[1]+1]]= shiftedMass
			myIC[[i-scanRange[1]+1]]= shiftedIC
			#myScans[[i-scanRange[1]+1]]=scans
			#myScans_index[[i-scanRange[1]+1]]=scans_index

			# boolean is false for IC's whose sum !=0
			voidTracker=c(voidTracker, FALSE)
		} # ending first condition
			
		# second checking, true for scans with sum(Ic)=0,
		if((sum(fullIc[valdex])==0)){

			voidTracker=c(voidTracker, TRUE)
		}
	}# Main for loop ends


	# stores the scan length
	voidIndex=c(1:length(c(scanRange[1]:scanRange[2])))

	# returns myScans as a list
	#for (i in 1:length(myMass)) {scan_list=c(scan_list,myScans[[i]])}

# Removing null elements from myMass
# Counts the number of non-null values in myMass
#myMass=myMass[-(which(sapply(myMass,is.null),arr.ind=TRUE))]
#myIC=myIC[-(which(sapply(myIC,is.null),arr.ind=TRUE))]

# stores scans with no masses in users range
noMassScan=which(sapply(myMass,is.null),arr.ind=TRUE)

# filtering the processed scans , this can be a single line using which and is.null but it doesn't seems to be working somehow
a<-c(); d<-c() # real scan indices
for (i in 1:length(myMass)) {a[i]=is.null(myMass[[i]])}
for (i in 1:length(a)) {if (a[i]==FALSE) {d=c(d,i)}}
scany=c(scanRange[1]:scanRange[2])
for (i in 1:length(d)) {scans_in=c(scans_in,scany[d[i]])}


# prints error for empty mass's scans
for (i in 1:length(noMassScan)) {
	if (length(noMassScan)!=0) {
	if (warnings==FALSE) {warn2=rbind(warn2,(paste("Warning : Scan",noMassScan[i]+scanRange[1]-1," has no mass/charge value in the specified massRange")))}
	else {warn2=rbind(warn2,print(paste("Warning : Scan",noMassScan[i]+scanRange[1]-1," has no mass/charge value in the specified massRange")))}
	}
}

# this loop is ony true for the sum(Ic)=0, (true values count to sum)
if(sum(voidTracker)>0){
	for(i in 1:sum(voidTracker)){
		myIC[[(voidIndex[voidTracker])[i]]]=rep(0, length(myIC[[(voidIndex[!voidTracker])[1]]]))
		myMass[[(voidIndex[voidTracker])[i]]]=myMass[[(voidIndex[!voidTracker])[1]]]
	}
}

# removes first and last element and calculates the average of the rest
for(i in 1:length(myMass)){
	myMass[[i]]=(myMass[[i]][-1]+myMass[[i]][-length(myMass[[i]])])/2
}

# initial values are assigned to the minLength and template variables
minLength=length(myMass[[1]]); template=myMass[[1]];
for(i in 1:length(myMass)){
	if(length(myMass[[i]])<minLength){
		minLength=length(myMass[[i]])
		template=myMass[[i]]
	}
}

for(i in 1:length(myMass)){
	if(length(myMass[[i]])>minLength){
		starterInd=which.min(abs(myMass[[i]]-template[1]))
		myMass[[i]]=myMass[[i]][c(starterInd:(starterInd+minLength-1))]
		myIC[[i]]=myIC[[i]][c(starterInd:(starterInd+minLength-1))]
	}	
}
	
temp=length(myMass[[1]])
OIC=(myIC[[1]])
# this loop prints crud for masses, in which the length of first myMass list is not equal to the length of consequtive list
for(i in 2:length(myMass)) {
	if(temp!=length(myMass[[i]])){
		print("crud!"); print(i)
	}
	temp=(length(myMass[[i]]))
	# OIC row binds the OIC and myIC for each scan
	OIC=rbind(OIC,(myIC[[i]]))
}	
# this loop will run if plotit is true and will plot a graph
if(plotit) {
	colourPower=1
	cols=rgb(cbind((c(0:1000)/1000)^ colourPower,(4*c(0:1000)/1000*(1-c(0:1000)/1000))^1, (1-c(0:1000)/1000)^ colourPower+(c(0:1000)/1000)^(50*colourPower),rep(0,length(c(0:1000)))))
	image(z=OIC,x=c(scanRange[1]:scanRange[2]),y= template,ylab="m/z",xlab="rt",col=cols)
}
# peakReader1 function returns a list of OIC and template which can be sent to the heat function
return(list(ic=OIC,mass=template,scans=c(scanRange[1]:scanRange[2]),mz=massRange,rt=rt,warning1=warn1,warning2=warn2))
}

