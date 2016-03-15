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
# 		init.values=c(log(0.05),log(0.5),log(0.5),-0.05),
#		marginal="Weibull",copula="Clayton")


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
#' \item Clayton & Piecewise Exponential: \deqn{\lambda's, \theta and \beta's}
#' }
#' @param marginal Which marginal to use? (ADD OPTIONS)
#' @param copula Which copula to user? (ADD OPTIONS)
#' @return Not yet final, so far only dataframe with estimates and standard errors. (NOTE: EXAMPLES CAN LATER BE RUN)
#' @examples
#' \dontrun{
#' data("insem",package="Sunclarco")
#' result <- CopulaModel_1stage(data=insem,time="Time",status="Status",cluster="Herd",covariates="Heifer",
#' 		init.values=c(log(0.05),log(0.5),log(0.5),-0.05),
#' 		marginal="Weibull",copula="Clayton")
#' result
#' }
CopulaModel_1stage <- function(data,time,status,cluster,covariates,init.values,marginal="Weibull",copula="Clayton"){
	### Check Dataframe
	
	# TO DO
	
	### Check parameters  (check if combination possible, length of init.values)
	
	# TO DO
	
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
		cutpoints <- quantile(data_time[,time][data_time[,status]==1], probs = seq(0, 1, by = 0.05),names=FALSE)
		cutpoints[1] <- 0
		cutpoints[length(cutpoints)] <- max(cutpoints)+1000  
		num_pieces <- length(cutpoints)-1 #number of intervals between cutpoints
		
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
		
		#standard errors for lambda, rho, beta, theta
		stderr_all <- sqrt(diag(VarCov1_weibCL))
		
	
		#kendall's tau
		tau1_weibCL <- theta1_weibCL/(2+theta1_weibCL) #0.09600836
		stderrtau1_weibCL <- (2/(theta1_weibCL+2)^2)*stderrtheta1_weibCL #0.006113899
		
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambda1_weibCL,rho1_weibCL,theta1_weibCL,betas_weibCL,tau1_weibCL),StandardErrors=c(stderr_all,stderrtau1_weibCL))
		rownames(out1) <- c("lambda","rho","theta",paste0("beta_",covariates),"Kendall_Tau")
		

		return(out1)
	}
	
	if(marginal=="PiecewiseExp" & copula=="Clayton"){
		
		
		
		
	}
	

}


##################################################
##################################################



