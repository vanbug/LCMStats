coates <-
function(X, Np){
	if(!is.matrix(X)){
		stop("X must be a matrix")
		}
		p=matrix(ncol=ncol(X),nrow=nrow(X))
		p[,1]=X[,1]/Np
		p[,2]=X[,2]/(Np*(1-p[,1]))
		for(i in 3:ncol(X)){
			p[,i]=X[,i]/(Np*apply((1-p[,1:(i-1)]),1,prod))
			}
		p[is.nan(p)]=0.5
		p[p>1]=1
		p[p<0]=0
		X=apply(-Np*log(1-p), 1, sum)
	return(X)
	}

