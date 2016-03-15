# Project: Copula Package
# 
# Author: lucp8394
###############################################################################


########################################################
## 1-stage: Clayton copula model with Weibull margins ##
########################################################


# p : lambda, rho, theta, beta's 
# NOTE THAT ORDER WAS CHANGED!!!!

loglik.1stage_weibCL <- function(p,data,covariates,status,time,cluster,ClusterData,ClusterDataList){
	
	varphi         <- function(theta,t){(1+theta*t)^(-1/theta)}
	varphi.inverse <- function(theta,t){(t^(-theta)-1)/theta}
	varphi.prime   <- function(theta,t){-(1+theta*t)^(-1/theta-1)}
	
	lambda <- exp(p[1])
	rho    <- exp(p[2])
	theta  <- exp(p[3])
		
	betas   <- p[4:length(p)]
	names(betas) <- covariates
	
	cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){betas %*% t(matrix(row,nrow=1,ncol=length(row)))})
	
	s <- exp(-lambda*exp(cov.lincomb)*data[,time]^rho) #u_ij
	f <- lambda*rho*data[,time]^(rho-1)*exp(cov.lincomb)*s #-du_ij/dy_ij
	
	
	sumG <- 1:length(ClusterDataList)
	sumH <- 1:length(ClusterDataList)
	
	for(i in 1:length(ClusterDataList)){
		ClusterDataList[[i]]$S <- s[data[,cluster]==i]
		ClusterDataList[[i]]$F <- f[data[,cluster]==i]
		
		ClusterDataList[[i]]$G <- ClusterDataList[[i]][,status]*log(-ClusterDataList[[i]]$F/varphi.prime(theta,varphi.inverse(theta,ClusterDataList[[i]]$S)))
		
		
#		C[[i]]*log(-F[[i]]/varphi.prime(theta,varphi.inverse(theta,S[[i]])))
		
		ClusterDataList[[i]]$H <- varphi.inverse(theta,ClusterDataList[[i]]$S)
		
		sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1])
		sumH[i] <- sum(ClusterDataList[[i]]$H )
	}
	
	loglik <- sapply(ClusterData$ClusterEvents,function(x) ifelse(x==0,0,sum(log(1/theta+seq(0,x-1)))))-
			(1/theta)*log(theta)+ sumG + (-ClusterData$ClusterEvents-1/theta)*log(sumH+1/theta);
	
	-sum(loglik)
}



loglik.1stage_pweCL <- function(p,cutpoints,num_pieces,datapiece){
	
	lambdas <- exp(p[1:num_pieces])
	theta  <- exp(p[num_pieces+1])
	beta   <- p[num_pieces+2:length(p)] 
	
 	# I AM HERE, WHY DATAPIECE NEVER USED?
	# Remember that theta and beta are switched  + beta became betas
	
	haz <- approx(cutpoints[1:num_pieces],lambdas,xout=time,method='constant',rule=2)$y
	DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
	Inter <- c(0,cumsum(lambdas*DiffInter))
	cumhaz <- approx(cutpoints,Inter,xout=time,method='linear',rule=2)$y
	
	s <- exp(-cumhaz*exp(beta*Z)); #u_ij 
	f <- haz*exp(beta*Z)*s; #-du_ij/dy_ij
	
	S <- vector("list",length=K)
	F <- vector("list",length=K)
	
	for(i in 1:K){
		S[[i]] <- s[cluster==i]
		F[[i]] <- f[cluster==i]}
	
	G <- vector("list",K)
	H <- vector("list",K)
	sumG <- 1:K
	sumH <- 1:K 
	
	for (i in 1:K){
		G[[i]] <- log(-F[[i]]/varphi.prime(theta,varphi.inverse(theta,S[[i]])))
		H[[i]] <- varphi.inverse(theta,S[[i]])
		sumG[i] <- sum(G[[i]][C[[i]]==1]);
		sumH[i] <- sum(H[[i]])}
	
	loglik <- sapply(d,function(x) ifelse(x==0,0,sum(log(1/theta+seq(0,x-1)))))-
			(1/theta)*log(theta)+ sumG + (-d-1/theta)*log(sumH+1/theta);
	-sum(loglik)
}
