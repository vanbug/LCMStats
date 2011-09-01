pooler <-
function(counts, truth=NULL){
	
	if(!is.matrix(counts)){
		stop("counts is not a matrix")
		}

	
	if(!is.null(truth)){
	if(nrow(counts)!=length(truth)){
		stop("dimensions of counts do not match those of truth")
		}
	}
		
	#calculate a temporary estimate of binomial probability
	if(is.null(truth)){
		truth<-apply(counts[,apply(counts,2,sum)<100],1,sum)/sum(apply(counts[,apply(counts,2,sum)<100],1,sum))
		}
		
	#define counts that are to be pooled
	lowCounts=counts[,apply(counts,2,sum)*min(truth)<5]

	# in case of zooming
	if (length(lowCounts)==0||(is.na(truth)==TRUE)) {finalCounts=counts} 
	if (length(lowCounts)!=0&&length(ncol(lowCounts))!=0&&(is.na(truth)==FALSE)){sums=rep(0, nrow(counts))#keep track of how many data have been pooled
	newCounts=matrix(nrow=nrow(counts), ncol=0)#matrix to store the pooled data
	for(i in 1:ncol(lowCounts)){
		sums=sums+ lowCounts[,i]
		if(sum(sums)*min(truth)>=5){newCounts=cbind(newCounts, sums); sums=rep(0, nrow(counts))}
		}
	index=apply(counts,2,sum)*min(truth)>=5#index of counts that were high enough to begin with
	finalCounts=cbind(counts[,index], newCounts)
	}
		
	return(finalCounts)
}

