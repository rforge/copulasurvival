# Project: Copula Package
# 
# Author: lucp8394
###############################################################################


## EXAMPLE DATA ##
#' Insemination Data
#'
#' Description of the data.
#' Can have multiple lines. (Note that other fields like references, authors,... can be added)
#'
#' @format A dataframe with 10513 rows and 5 columns.
#' @name insem
NULL


## IMPORTS ##

# roxygen code to import 
#library(survival)
#library(matrixcalc)

#insem <- read.table(file="insemination.txt",header=T,sep="")
#head(insem)
#
#data <- insem
#time <- "Time"
#status <- "Status"
#cluster <- "Herd"
#covariates <- "Heifer"

#init.values <- c(log(0.05),log(0.5),log(0.5),-0.05)   # lambda, rho, theta, beta's
# result <- CopulaModel_1stage(data=insem,time="Time",status="Status",cluster="Herd",covariates="Heifer",
# 		init.values=init.values,
#		marginal="Weibull",copula="Clayton")
#result

#init.values <- c(rep(log(0.02),20),log(0.5),-0.05)
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",cluster="Herd",covariates="Heifer",
#		init.values=init.values,
#		marginal="PiecewiseExp",copula="Clayton")
#result  ## OKAY THAT IT TAKES LONG? SAME COMPUTATION TIME?

#init.values <- c(log(0.05),log(0.5),log(0.5/(1-0.5)),-0.05)
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",cluster="Herd",covariates="Heifer",
# 		init.values=init.values,
#		marginal="Weibull",copula="GH")
#result  
#
#init.values <- c(rep(log(0.02),10),log(0.5/(1-0.5)),-0.05)
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",cluster="Herd",covariates="Heifer",
# 		init.values=init.values,n.piecewise=10,
#		marginal="PiecewiseExp",copula="GH")
#result  # Takes pretty long!

#' Title of this function (do later, function still temporary)
#' 
#' Description of this function. (Note: you can also add other fields like references, authors,... See which ones are necessary)
#' 
#' @export
#' @param data Input dataframe containing all variables.
#' @param time Which variable name is the time covariate?
#' @param status Which variable name is the status covariate? (more explanation?)
#' @param cluster Which variable name is the cluster covariate?
#' @param covariates Vector of 1 or more other covariates to be included in the model. (NOTE: MULTIPLE COVARIATES NEED TO BE TESTED!!!)
#' @param init.values Initial values for parameters: 
#' \itemize{
#' \item Clayton & Weibull: \deqn{\lambda, \rho, \theta and \beta's}
#' \item Clayton & Piecewise Exponential: \deqn{\lambda's, \theta and \beta's} Note that the number of lambda's is equal to the number of pieces of the piecewise exponential.
#' }
#' @param marginal Which marginal to use? Can either be \code{"Weibull"} or \code{"PiecewiseExp"} for Piecewise Exponential.
#' @param copula Which copula to user? Can either be \code{"Clayton"} or \code{"GH"} for Gumbel-Hougaard.
#' @param n.piecewise For \code{marginal="PiecewiseExp"}, how many pieces should the Piecewise Exponential have? (Default = 20)
#' @return Not yet final, so far: list object with: parameters (est & se) + kendal + cov matrix. Add more information?(NOTE: EXAMPLES CAN LATER BE RUN)
#' @examples
#' \dontrun{
#' data("insem",package="Sunclarco")
#' result <- CopulaModel_1stage(data=insem,time="Time",status="Status",cluster="Herd",covariates="Heifer",
#' 		init.values=c(log(0.05),log(0.5),log(0.5),-0.05),
#' 		marginal="Weibull",copula="Clayton")
#' result
#' }
CopulaModel_1stage <- function(data,time,status,cluster,covariates,init.values,marginal="Weibull",copula="Clayton",n.piecewise=20){
	### Check Dataframe
	
	# TO DO: check if data is data.frame + time,status,cluster, covariates are characters -> Include in S4 method
	
	for(i in c(time,status,cluster,covariates)){
		if(!(i %in% colnames(data))){stop(paste0("\"",i,"\" is not a variable of the data frame."),call.=FALSE)}
	}
	
	### Check parameters  (check if combination possible, length of init.values (lambda's and betas!) )
	if(!(marginal %in% c("Weibull","PiecewiseExp"))){stop(paste0("Parameter 'marginal' can not be ",marginal,". It should be either \"Weibull\" or \"PiecewiseExp\""),call.=FALSE)}
	if(!(copula %in% c("Clayton","GH"))){stop(paste0("Parameter 'copula' can not be ",copula,". It should be either \"Clayton\" or \"GH\""),call.=FALSE)}
	
	if(marginal=="Weibull"){
		correct_length <- 3+length(covariates)
	}
	else if(marginal=="PiecewiseExp"){
		correct_length <- n.piecewise+1+length(covariates)
	}
	if(length(init.values)!=correct_length){stop(paste0("Parameter 'init.values' has an incorrect length. With given parameters it should be of length ",correct_length),call.=FALSE)}
	
	
	### CLUSTER VARIABLE NEEDS TO BE CHECKED! IT SHOULD BE A NUMBER FROM 1 TO N. IF NEEDED IT HAS TO BE TRANSFORMED!
	cluster_labels <- sort(unique(data[,cluster]))
	if(!identical(cluster_labels,c(1:length(cluster_labels)))){
		
		newlabels <- c(1:length(cluster_labels))
		
		temp <- sapply(data[,cluster],FUN=function(x){
					index <- which(x == cluster_labels)
					return(newlabels[index])
				})
		data[,cluster] <- temp
	}
	
	### Sort Data: Cluster (needed for rle)
	data <- data[order(data[,cluster]),]
				
	
	### Get Information
#	Time <- data[,time]
#	Status <- data[,status]
	Cluster <- data[,cluster]
#	Covariates <- data[,covariates,drop=FALSE]  # keep structure
	
	nclusters <- length(levels(as.factor(Cluster)))
	
	ClusterData <- data.frame(Cluster=c(1:nclusters),ClusterSubjects=rle(Cluster)$lengths,ClusterEvents=aggregate(data[,status],by=list(Cluster),FUN=sum)[,2])

	ClusterDataList <- vector("list",length=nclusters)
	
	for(i.list in 1:nclusters){
		ClusterDataList[[i.list]] <- data[Cluster==i.list,c(time,status,covariates)]
	}
	
	if(marginal=="PiecewiseExp"){
		
		#define cutpoints for piecewise exponential
		
		# Sort data on time
		data_time <- data[order(data[,time]),]
		#define cutpoints such that all intervals contain same number of events
		cutpoints <- quantile(data_time[,time][data_time[,status]==1], probs = seq(0, 1, length.out=(n.piecewise+1)),names=FALSE)
		cutpoints[1] <- 0
		cutpoints[length(cutpoints)] <- max(cutpoints)+1000  
		num_pieces <- length(cutpoints)-1 #number of intervals between cutpoints
		
		
		# DELETE FOLLOWING LINE IF NEVER USED !!
		datapiece <- split(data_time, cut(data_time[,time], cutpoints, include.lowest=TRUE))
	}
	
	
	
	
	if(marginal=="Weibull" & copula=="Clayton"){ ## THIS IS ONLY TEMPORARY, NEED TO CHECK IF BETTER OVERLAP WITH OTHER METHODS
		
		res1_weibCL <- optim(init.values,
				loglik.1stage_weibCL,
				data=data,time=time,status=status,cluster=cluster,covariates=covariates,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				hessian=TRUE,control=list(maxit=3000)) 
		
		# Get Estimates and put in dataframe
		lambda1_weibCL <- exp(res1_weibCL$par[1]) 
		rho1_weibCL    <- exp(res1_weibCL$par[2]) 
		theta1_weibCL  <- exp(res1_weibCL$par[3]) 
		
		betas_weibCL   <- res1_weibCL$par[4:length(res1_weibCL$par)] 
		
		
		# Preparation for standard errors
		VarCov1_weibCL <- diag(c(lambda1_weibCL,rho1_weibCL,theta1_weibCL,rep(1,length(betas_weibCL))))%*%solve(res1_weibCL$hessian)%*%diag(c(lambda1_weibCL,rho1_weibCL,theta1_weibCL,rep(1,length(betas_weibCL))))
		
		stderrtheta1_weibCL <- sqrt(VarCov1_weibCL[3,3])
		
		#standard errors for lambda, rho,theta, betas
		stderr_all <- sqrt(diag(VarCov1_weibCL))
		
	
		#kendall's tau
		tau1_weibCL <- theta1_weibCL/(2+theta1_weibCL) #0.09600836
		stderrtau1_weibCL <- (2/(theta1_weibCL+2)^2)*stderrtheta1_weibCL #0.006113899
		
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambda1_weibCL,rho1_weibCL,theta1_weibCL,betas_weibCL),StandardErrors=c(stderr_all))
		rownames(out1) <- c("lambda","rho","theta",paste0("beta_",covariates))
		
		out2 <- c(tau1_weibCL,stderrtau1_weibCL)
		names(out2) <- c("Estimates","StandardErrors")
		
		out3 <- VarCov1_weibCL
		rownames(out3) <- colnames(out3) <- rownames(out1)
		

	}
	
	if(marginal=="PiecewiseExp" & copula=="Clayton"){
		
		# Why different method here? Should this be decidable for the user?
		res1_pweCL <- optim(init.values,
				loglik.1stage_pweCL,cutpoints=cutpoints,num_pieces=num_pieces,data=data,status=status,time=time,covariates=covariates,cluster=cluster,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				hessian=TRUE,method="BFGS")
		
		# Get Estimates and put in dataframe
		lambdas1_pweCL <- exp(res1_pweCL$par[1:num_pieces]) 
		theta1_pweCL <- exp(res1_pweCL$par[num_pieces+1]) 
		betas_pweCL <- res1_pweCL$par[(num_pieces+2):length(res1_pweCL$par)] 
		
		# Preparation for standard errors
		VarCov1_pweCL <- diag(c(lambdas1_pweCL,theta1_pweCL,rep(1,length(betas_pweCL))))%*%solve(res1_pweCL$hessian)%*%diag(c(lambdas1_pweCL,theta1_pweCL,rep(1,length(betas_pweCL))))
		stderrtheta1_pweCL <- sqrt(VarCov1_pweCL[num_pieces+1,num_pieces+1])

		#standard errors for lambdas, theta, betas
		stderr_all <- sqrt(diag(VarCov1_pweCL))
	
		#kendall's tau
		tau1_pweCL <- theta1_pweCL/(2+theta1_pweCL)
		stderrtau1_pweCL <- (2/(theta1_pweCL+2)^2)*stderrtheta1_pweCL
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambdas1_pweCL,theta1_pweCL,betas_pweCL),StandardErrors=c(stderr_all))
		rownames(out1) <- c(paste0("lambda",c(1:length(lambdas1_pweCL))),"theta",paste0("beta_",covariates))
		
		out2 <- c(tau1_pweCL,stderrtau1_pweCL)
		names(out2) <- c("Estimates","StandardErrors")
				
		out3 <- VarCov1_pweCL
		rownames(out3) <- colnames(out3) <- rownames(out1)
		

	}
	
	
	if(marginal=="Weibull" & copula=="GH"){
		
#		res1_weibGH <- optim(init.values,
#				loglik.1stage_weibGH,data=data,time=time,status=status,cluster=cluster,covariates=covariates,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
#				hessian=TRUE,control=list(maxit=3000)) 
		
		res1_weibGH <- optim(init.values,
				loglik.1stage_GH,data=data,cutpoints=NULL,num_pieces=NULL,time=time,status=status,cluster=cluster,covariates=covariates,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				marginal=marginal,
				hessian=TRUE,control=list(maxit=3000)) 
		
		# Get Estimates and put in dataframe
		lambda1_weibGH <- exp(res1_weibGH$par[1]) 
		rho1_weibGH    <- exp(res1_weibGH$par[2]) 
		theta1_weibGH  <- exp(res1_weibGH$par[3])/(1+exp(res1_weibGH$par[3]))
		betas1_weibGH   <- res1_weibGH$par[4:length(res1_weibGH$par)] 
		
		# Preparation for standard errors
		VarCov1_weibGH <- diag(c(lambda1_weibGH,rho1_weibGH,exp(res1_weibGH$par[3])/(1+exp(res1_weibGH$par[3]))^2,rep(1,length(betas1_weibGH))))%*%solve(res1_weibGH$hessian)%*%diag(c(lambda1_weibGH,rho1_weibGH,exp(res1_weibGH$par[3])/(1+exp(res1_weibGH$par[3]))^2,rep(1,length(betas1_weibGH))))
		stderrtheta1_weibGH <- sqrt(VarCov1_weibGH[3,3]) 
		
		#standard errors for lambda, rho, theta, betas
		stderr_all <- sqrt(diag(VarCov1_weibGH))
		
		#kendall's tau
		tau1_weibGH <- 1-theta1_weibGH
		stderrtau1_weibGH <- stderrtheta1_weibGH
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambda1_weibGH,rho1_weibGH,theta1_weibGH,betas1_weibGH),StandardErrors=c(stderr_all))
		rownames(out1) <- c("lambda","rho","theta",paste0("beta_",covariates))
	
		out2 <- c(tau1_weibGH,stderrtau1_weibGH)
		names(out2) <- c("Estimates","StandardErrors")
		
		out3 <- VarCov1_weibGH
		rownames(out3) <- colnames(out3) <- rownames(out1)
		
	}
	
	if(marginal=="PiecewiseExp" & copula=="GH"){
#		res1_pweGH <- optim(init.values,
#				loglik.1stage_pweGH,,cutpoints=cutpoints,num_pieces=num_pieces,data=data,status=status,time=time,covariates=covariates,cluster=cluster,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
#				hessian=TRUE,method="BFGS")

		res1_pweGH <- optim(init.values,
				loglik.1stage_GH,,cutpoints=cutpoints,num_pieces=num_pieces,data=data,status=status,time=time,covariates=covariates,cluster=cluster,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				marginal=marginal,
				hessian=TRUE,method="BFGS")
		
		# Get Estimates and put in dataframe
		lambdas1_pweGH <- exp(res1_pweGH$par[1:num_pieces])
		theta1_pweGH <- exp(res1_pweGH$par[num_pieces+1])/(1+exp(res1_pweGH$par[num_pieces+1]))  
		betas1_pweGH <- res1_pweGH$par[(num_pieces+2):length(res1_pweGH$par)] 

		# Preparation for standard errors
				
		VarCov1_pweGH <- diag(c(lambdas1_pweGH,exp(res1_pweGH$par[num_pieces+1])/(1+exp(res1_pweGH$par[num_pieces+1]))^2,rep(1,length(betas1_pweGH))))%*%solve(res1_pweGH$hessian)%*%diag(c(lambdas1_pweGH,exp(res1_pweGH$par[num_pieces+1])/(1+exp(res1_pweGH$par[num_pieces+1]))^2,rep(1,length(betas1_pweGH))))
		stderrtheta1_pweGH <- sqrt(VarCov1_pweGH[num_pieces+1,num_pieces+1]) 
		
		#standard errors for lambdas, theta, betas
		stderr_all <- sqrt(diag(VarCov1_pweGH))
				
		#kendall's tau
		tau1_pweGH <- 1-theta1_pweGH
		stderrtau1_pweGH <- stderrtheta1_pweGH
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambdas1_pweGH,theta1_pweGH,betas1_pweGH),StandardErrors=c(stderr_all))
		rownames(out1) <- c(paste0("lambda",c(1:length(lambdas1_pweGH))),"theta",paste0("beta_",covariates))
		
		out2 <- c(tau1_pweGH,stderrtau1_pweGH)
		names(out2) <- c("Estimates","StandardErrors")
		
		out3 <- VarCov1_pweGH
		rownames(out3) <- colnames(out3) <- rownames(out1)
	}
	
		out <- list(Parameters=out1,Kendall_Tau=out2,ParametersCov=out3)
		class(out) <- "Sunclarco"
		return(out)

}


##################################################
##################################################



