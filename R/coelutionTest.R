coelutionTest <-
function(X, Y, thresh=300, robust=0.05, truth=NA, Np=NA, plotit=TRUE,pooler=1){

	#sanity check
	if(!is.list(X)|!is.list(Y)){
		stop("X and Y must be lists")
		}
	if(!is.numeric(thresh)){
		stop("thresh is not a numeric")
		}
	if(!is.numeric(robust)){
		stop("robust is not a numeric")
		}
	if(!is.na(truth)&!is.numeric(truth)){
		stop("truth is not a numeric")
		}
	if(!is.na(Np)&!is.numeric(Np)){
		stop("Np is not a numeric")
		}
	if(!is.logical(plotit)){
		stop("plotit is not a boolean")
		}
	
	clusterCounts<-list()#list of centroided data

	#centroid continuum data and apply Coates' if required
	if(is.matrix(X$ic)){
		if(is.numeric(Np)){
			clusterCounts[[1]]<-coates(X$ic, Np)
			}
		if(!is.numeric(Np)){
			clusterCounts[[1]]<-apply(X$ic,1,sum)
			}	
		}
	if(is.matrix(Y$ic)){
		if(is.numeric(Np)){
			clusterCounts[[2]]<-coates(Y$ic, Np)
			}
		if(!is.numeric(Np)){
			clusterCounts[[2]]<-apply(Y$ic,1,sum)
			}	
		}	
	
	if(is.vector(X$ic)){
		clusterCounts[[1]]<-X$ic
		}
	if(is.vector(Y$ic)){
		clusterCounts[[2]]<-Y$ic
		}
	
	#find the scan-range common to both data-sets
	minScan<-max(min(X$scans), min(Y$scans)); maxScan<-min(max(X$scans), max(Y$scans))

	#reduce data to common scan-range
	clusterCounts[[1]]<- clusterCounts[[1]][(X$scans>=minScan)&(X$scans<=maxScan)]
	clusterCounts[[2]]<- clusterCounts[[2]][(Y$scans>=minScan)&(Y$scans<=maxScan)]
	
	# removing null elements from clusterCounts as sometimes zooming from heat plot indices it
	if (length(clusterCounts[[1]][-(which(sapply(clusterCounts[[1]],is.na),arr.ind=TRUE))])!=0) {
		clusterCounts[[1]]=clusterCounts[[1]][-(which(sapply(clusterCounts[[1]],is.na),arr.ind=TRUE))]
		clusterCounts[[2]]=clusterCounts[[2]][-(which(sapply(clusterCounts[[1]],is.na),arr.ind=TRUE))]
	}
	if (length(clusterCounts[[2]][-(which(sapply(clusterCounts[[2]],is.na),arr.ind=TRUE))])!=0) {
		clusterCounts[[1]]=clusterCounts[[1]][-(which(sapply(clusterCounts[[2]],is.na),arr.ind=TRUE))]
		clusterCounts[[2]]=clusterCounts[[2]][-(which(sapply(clusterCounts[[2]],is.na),arr.ind=TRUE))]
	}
	clusterMat<-matrix(c(clusterCounts[[1]], clusterCounts[[2]]), nrow=2, byrow=TRUE)
		
	if(!is.na(truth[1])){
		if (pooler==1) {clusterMat<-pooler(clusterMat, truth) }#pool data
		if (pooler==2) {clusterMat=clusterMat}
		valIndex<-c(1:ncol(clusterMat))[(apply(clusterMat,2,sum)<=thresh)]#find non-saturated data
		nCounts<-apply(clusterMat,2,sum)#sum of counts in each scan
		est<-truth[1]/sum(truth)#define binomial probability
		}
	if(is.na(truth[1])){
		if (pooler==1) {clusterMat<-pooler(clusterMat)}#pool data
		if (pooler==2) {clusterMat=clusterMat}
		valIndex<-c(1:ncol(clusterMat))[(apply(clusterMat,2,sum)<=thresh)]#find non-saturated data
		nCounts<-apply(clusterMat,2,sum)#sum of counts in each scan
		est<-sum(clusterMat[1, valIndex])/sum(clusterMat[, valIndex])#define binomial probability
		}

	#calculate chi2 statistics and apply robustness step
	chi<-((clusterMat[1,valIndex]-est* nCounts[valIndex]))^2/(nCounts[valIndex]* est*(1-est))
	chiRbstIndex<-order(chi, decreasing=FALSE)[c(1:round(length(valIndex)*(1-robust)))]
	chiRbst<-chi[chiRbstIndex]
	
	#calculate p-values for individual data points and for combined data-set
	pVals<-1-pchisq(chiRbst, df=1)
	pVal<-1-pchisq(sum(chiRbst), df=length(chiRbst)-1)
	
	if(plotit){
		
		#define confidence interval by inverting GOF test
		z<-seq(0,2*max(clusterMat[1,]),length.out=1000)
		chi99top<-sqrt((1-est)*qchisq(0.99,1)*((4*z)+(1-est)*qchisq(0.99,1))/(2* est))-((est-1)*(2*z+qchisq(0.99,1)))/(2* est)
		chi99bot<--sqrt((1-est)*qchisq(0.99,1)*((4*z)+(1-est)*qchisq(0.99,1))/(2* est))-((est-1)*(2*z+qchisq(0.99,1)))/(2* est)

		#define colour scheme
		colSigRGB<-matrix(rep(0,3*length(clusterMat[1,])),ncol=3)
		colSigRGB[valIndex[chiRbstIndex],]<-cbind(pchisq((chi[chiRbstIndex]),1)^1, 4*pchisq((chi[chiRbstIndex]),1)*(1-pchisq((chi[chiRbstIndex]),1)), (1-pchisq((chi[chiRbstIndex]),1))^1+pchisq((chi[chiRbstIndex]),1)^100)
		pValueChiCols <-rgb(colSigRGB)
		
		plot(clusterMat[1,], clusterMat[2,],col= pValueChiCols, main=paste("GOF p-value = ",round(pVal,4)), xlim=c(0,max(clusterMat[1,])), ylim=c(0,max(clusterMat[2,])), xlab="ion count", ylab="ion count")
		lines(z,chi99top,lty=2)
		lines(z,chi99bot,lty=2)
		abline(a=0,b=(est ^-1-1))
			
		}

	return(list(p.value=pVal, robust.chis= chiRbst))

	}

