qqchisq <-
function(x, df=1, new=TRUE, main="",sub=""){

if(sum(is.na(x))!=0){
	print("warning: NAs present")
	}

x=x[!is.na(x)]
n=length(x)

if(new){
	plot(sort(x),qchisq(c(1:n)/(n+1),df),ylim=c(min(x,qchisq(c(1:n)/(n+1),df)),max(x,qchisq(c(1:n)/(n+1),df))), xlim=c(min(x,qchisq(c(1:n)/(n+1),df)),max(x,qchisq(c(1:n)/(n+1),df))), xlab="Empirical Quantiles", ylab="Theoretical Quantiles", main=main,sub=sub)}
if(!new){
	points(sort(x),qchisq(c(1:n)/(n+1),df), col="red")
	}

abline(coef=c(0,1), lty=2, lwd=2, col="red")
}

