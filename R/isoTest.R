isoTest <-
function(dataList, isoPat, thresh=300, robust=0.1, Np=NA){
	
	#sanity check
	if(!is.na(Np)&!is.numeric(Np)){
		stop("Np must be either NA or numeric")
		}		
	if(!is.numeric(thresh)){
		stop("thresh is not a numeric")
		}
	if(!is.list(isoPat)){
		stop("isopat is not a list")
		}
	if(!is.list(dataList)){
		stop("dataList is not a list")
		}
	if(!is.numeric(robust)){
		stop("robust is not a numeric")
		}

	nCluster<-length(isoPat)#the number of clusters of isotopologues
	df=c()#degree of freedom for the various clusters
	chis=list()#chi2 statistics for individual scans
	chisRbst=c()#chi2 statistics after the specified fration of largest values has been removed
	chiRbst=c()#sum of chi2 statistics for each cluster of isotopologues
	sSize=c()#number of effective scans for each cluster
	
	#calculate chi2 statistics for each cluster sequentially
	for(i in 1:nCluster){
		df[i]<-length(isoPat[[i]])-1
		clusterCounts<-list()#list of the centroided data from this cluster
		
		isoPat[[i]]<-isoPat[[i]]/sum(isoPat[[i]])#normalise isotope pattern

		minScan=0; maxScan=Inf#used to keep track of scan-ranges in order to find the largest scan range common to all data-sets in cluster
					
		for(j in 1:(df[i]+1)){
			#centroid continuum data and apply Coates' if required
			if(is.matrix(dataList[[i]][[j]][[1]])){
				if(is.numeric(Np)){
					clusterCounts[[j]]<-coates(dataList[[i]][[j]][[1]], Np)
					}
				if(is.na(Np)){
					clusterCounts[[j]]<-apply(dataList[[i]][[j]][[1]],1,sum)
				}	
			}
			if(is.vector(dataList[[i]][[j]][[1]])){
					clusterCounts[[j]]=dataList[[i]][[j]][[1]]
				}

	
			#find the largest scan-range common to all data-sets
			if(min(dataList[[i]][[j]]$scans)>minScan){
				minScan =min(dataList[[i]][[j]]$scans)
				}
			if(max(dataList[[i]][[j]]$scans)<maxScan){
				maxScan =max(dataList[[i]][[j]]$scans)
				}		
		}
		
		#reduce data to common scan-range
		for(j in 1:(df[i]+1)){
			clusterCounts[[j]]<-clusterCounts[[j]][(dataList[[i]][[j]]$scans>=minScan)&(dataList[[i]][[j]]$scans<=maxScan)]
		}	
		
		#turn data into matrix
		clusterMat=c()
		for(j in 1:(df[i]+1)){
			clusterMat=rbind(clusterMat, clusterCounts[[j]])
			}
		
		# pool together scans with low counts in order for chi2 test to be valid
		clusterMat=pooler(clusterMat, isoPat[[i]])
	
		#remove high ion counts
		clusterMat= clusterMat[,apply(clusterMat,2,sum)<=thresh]
	
		#sum of counts in each scan
		nCounts=apply(clusterMat,2,sum)
		
		#calculate and sort chi2 statistics
		chis[[i]]=sort(apply((clusterMat-outer(isoPat[[i]], nCounts))^2/(outer(isoPat[[i]], nCounts)),2,sum))
		
		#robustness step
		chisRbst[[i]]=chis[[i]][c(1:round(length(chis[[i]])*(1-robust)))]
		
		#pool statistics from across chromatographic scans
		chiRbst[i]= sum(chisRbst[[i]])
		
		#define effective number of chromatographic scans in order to calculate degree of freedom for final chi2 statistic
		sSize[i]=length(chisRbst[[i]])
	}
	
	#calculate final p-value
	pVal=1-pchisq(sum(chiRbst), df=sum(sSize*df))
	
	return(list(p.value=pVal, robust.chis= chisRbst, raw.chis= chis))
}

