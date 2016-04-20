# Project: Copula Package
# 
# Author: Gebruiker
###############################################################################


#Robust variance via grouped jackknife estimator

jack_2stage_pweCL <- function(data,covariates,status,time,clusters,ClusterData,ClusterDataList,num_pieces,init.values,cutpoints){
	
#	beta2jack_pwe   <- 1:length(DataClustList)
	betas2jack_pwe   <-  matrix(NA,nrow=length(ClusterDataList),ncol=length(covariates))
	lambdas2jack_pwe <- matrix(NA,nrow=length(ClusterDataList),ncol=num_pieces)
	theta2jack_pweCL <- 1:length(ClusterDataList)
	
	
	# some logllh
	loglik2.pwe_stage1jack <- function(p,num_pieces,status,time.k,Z.k,cutpoints){
		
		lambdas.k <- exp(p[1:num_pieces])
		betas.k <- p[(num_pieces+1):length(p)]
		
		
		haz.k <- approx(cutpoints[1:num_pieces],lambdas.k,xout=time.k,method='constant',rule=2)$y
		DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
		Inter.k <- c(0,cumsum(lambdas.k*DiffInter))
		cumhaz.k <- approx(cutpoints,Inter.k,xout=time.k,method='linear',rule=2)$y
		
		cov.lincomb <- apply(Z.k,MARGIN=1,FUN=function(row){betas.k %*% t(matrix(row,nrow=1,ncol=length(row)))})
		
		s.k <- exp(-cumhaz.k*exp(cov.lincomb)); #u_ij 
		f.k <- haz.k*exp(cov.lincomb)*s.k; #-du_ij/dy_ij
		
		loglik <- data.k[,status]*log(f.k)+(1-data.k[,status])*log(s.k)
		return(-sum(loglik))
	}
	
	
	loglik.2stagejack_pweCL <- function(p,ClusterData,ClusterDataList,status){
		
		theta  <- exp(p)
		
		
		sumG <- 1:length(ClusterDataList)
		sumH <- 1:length(ClusterDataList)
		
		for (i in 1:length(ClusterDataList)){
			ClusterDataList[[i]]$G <- ClusterDataList[[i]][,status]*log(-1/varphi.prime(theta,varphi.inverse(theta,ClusterDataList[[i]]$Sk)))
			ClusterDataList[[i]]$H <- varphi.inverse(theta,ClusterDataList[[i]]$Sk)
			sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1]);
			sumH[i] <- sum(ClusterDataList[[i]]$H)
		}
		
		loglik <- sapply(ClusterData$ClusterEvents,function(x) ifelse(x==0,0,sum(log(1/theta+seq(0,x-1)))))-
				(1/theta)*log(theta)+ sumG + (-ClusterData$ClusterEvents-1/theta)*log(sumH+1/theta)
		return(-sum(loglik[-k]))
		
	}
	
	# applying it
	for(k in 1:length(ClusterDataList)){
		
		data.k <- data[data[,clusters]!=k,]
		Z.k <- data.k[,covariates,drop=FALSE]
		time.k <- data.k[,time]
		
		
		res2.pwe_stage1jack <- optim(init.values[-(num_pieces+1)],loglik2.pwe_stage1jack,
				num_pieces=num_pieces,status=status,time.k=time.k,Z.k=Z.k,cutpoints=cutpoints,
				method="BFGS")
		
		lambdas2jack_pwe[k,] <- exp(res2.pwe_stage1jack$par[1:num_pieces])
		betas2jack_pwe[k,] <- res2.pwe_stage1jack$par[(num_pieces+1):length(res2.pwe_stage1jack$par)]
		
		
		#inpluggen in Clayton loglikelihood
		haz.kk <- approx(cutpoints[1:num_pieces],lambdas2jack_pwe[k,],xout=time.k,method='constant',rule=2)$y
		DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
		Inter.kk <- c(0,cumsum(lambdas2jack_pwe[k,]*DiffInter))
		cumhaz.kk <- approx(cutpoints,Inter.kk,xout=time.k,method='linear',rule=2)$y
		
		cov.lincomb <- apply(Z.k,MARGIN=1,FUN=function(row){betas2jack_pwe[k,] %*% t(matrix(row,nrow=1,ncol=length(row)))})
		
		s.kk <- exp(-cumhaz.kk*exp(cov.lincomb))
		
		for(i in 1:length(ClusterDataList)){
			if(i!=k){
				ClusterDataList[[i]]$Sk <- s.kk[data.k[,clusters]==i]
			}
			else {
				ClusterDataList[[i]]$Sk <- rep(0.123,length(ClusterDataList[[i]][,time]))
			}
		}
		

		res2jack_pweCL <- optim(init.values[num_pieces+1],loglik.2stagejack_pweCL,
				ClusterData=ClusterData,ClusterDataList=ClusterDataList,status=status,
				method="BFGS")
		theta2jack_pweCL[k] <- exp(res2jack_pweCL$par) 
		
		
	}
	return(list(lambdas=lambdas2jack_pwe,theta=theta2jack_pweCL,betas=betas2jack_pwe))

}



#######################################################################################
#######################################################################################

