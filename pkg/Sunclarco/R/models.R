# Project: Copula Package
# 
# Author: lucp8394
###############################################################################



#Wat de optimalisatiemethode betreft: de BFGS is zeker nodig in de PWE context, dus daar zou ik die altijd gebruiken. 
#Ik ga nog uittesten of het zou lukken om overal BFGS te gebruiken.

## TO DO:
# - CHANGE INPUT to model="<copula>-<marginal>"
# - verbose voor jackknife? Add verbose in general (to know at which step you are) -> DO AT THE END
# - testing all methods for multiple betas
# - automatic initial values
# - clusterlabel code should be testen
# - fuse stage 1 and 2 together
# - correct parameters for optim
# - faster alternative than optim? -> bottleneck (maxLik, glm, nlm?)
# - covariance matrix for 2-stage piecewise clayton?
# - robust solution instead of jackknife 2stage piecewise clayton
# - ... (see questions)

# - add jackknife threshold when it can not be used!!


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

#' @importFrom survival survreg Surv cluster survfit coxph
NULL

# roxygen code to import 
#library(survival)
#library(matrixcalc)

#insem <- read.table(file="insemination.txt",header=T,sep="")
#head(insem)
#
#data <- insem


#init.values <- c(log(0.05),log(0.5),log(0.5),-0.05)   # lambda, rho, theta, beta's
# result <- CopulaModel_1stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
# 		init.values=init.values,
#		marginal="Weibull",copula="Clayton")
#result


#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
#			marginal="Weibull",copula="Clayton")
#result


#init.values <- c(rep(log(0.02),20),log(0.5),-0.05)
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
#		init.values=init.values,
#		marginal="PiecewiseExp",copula="Clayton")
#result  ## OKAY THAT IT TAKES LONG? SAME COMPUTATION TIME?
#
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
#		marginal="PiecewiseExp",copula="Clayton")
#result
#
#init.values <- c(log(0.05),log(0.5),log(0.5/(1-0.5)),-0.05)
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
# 		init.values=init.values,
#		marginal="Weibull",copula="GH")
#result  
#
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
#		marginal="Weibull",copula="GH")
#result  

#
#init.values <- c(rep(log(0.02),10),log(0.5/(1-0.5)),-0.05)
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
# 		init.values=init.values,n.piecewise=10,
#		marginal="PiecewiseExp",copula="GH")
#result  # Takes pretty long!
#
#result <- CopulaModel_1stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
#		n.piecewise=10,
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
#' @param clusters Which variable name is the cluster covariate?
#' @param covariates Vector of 1 or more other covariates to be included in the model. (NOTE: MULTIPLE COVARIATES NEED TO BE TESTED!!!)
#' @param init.values Initial values for parameters (if \code{NULL}, automatic choice, see Details for more information): 
#' \itemize{
#' \item Clayton & Weibull: \deqn{\lambda, \rho, \theta and \beta's}
#' \item Clayton & Piecewise Exponential: \deqn{\lambda's, \theta and \beta's} Note that the number of lambda's is equal to the number of pieces of the piecewise exponential.
#' }
#' @param marginal Which marginal to use? Can either be \code{"Weibull"} or \code{"PiecewiseExp"} for Piecewise Exponential.
#' @param copula Which copula to user? Can either be \code{"Clayton"} or \code{"GH"} for Gumbel-Hougaard.
#' @param n.piecewise For \code{marginal="PiecewiseExp"}, how many pieces should the Piecewise Exponential have? (Default = 20)
#' @return Not yet final, so far: list object with: parameters (est & se) + kendal + cov matrix. Add more information? MAYBE ADD CHOSEN PARAMETERS (CHOSEN MODEL)! (NOTE: EXAMPLES CAN LATER BE RUN)
#' @details SOME MORE INFORMATION ABOUT AUTOMATIC INITIAL PARAMETERS
#' @examples
#' \dontrun{
#' data("insem",package="Sunclarco")
#' result <- CopulaModel_1stage(data=insem,time="Time",status="Status",cluster="Herd",covariates="Heifer",
#' 		init.values=c(log(0.05),log(0.5),log(0.5),-0.05),
#' 		marginal="Weibull",copula="Clayton")
#' result
#' }
CopulaModel_1stage <- function(data,time,status,clusters,covariates,init.values=NULL,marginal="Weibull",copula="Clayton",n.piecewise=20){
	
	
	## DATA CHECK AND INPUT SAME FOR COPULAMODEL 1stage 2stage
	# - > when merging in bigger function, move this part to there!!
	
	
	### Check Dataframe
	
	# TO DO: check if data is data.frame + time,status,clusters, covariates are characters -> Include in S4 method
	
	for(i in c(time,status,clusters,covariates)){
		if(!(i %in% colnames(data))){stop(paste0("\"",i,"\" is not a variable of the data frame."),call.=FALSE)}
	}
	
	## Automatic initial parameters + naming parameters correctly
	if(is.null(init.values)){
		init.values <- estimate_parameters(data=data,time=time,status=status,clusters=clusters,covariates=covariates,n.piecewise=n.piecewise,marginal=marginal,copula=copula,stage=1)
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
	clusters_labels <- sort(unique(data[,clusters]))
	if(!identical(clusters_labels,c(1:length(clusters_labels)))){
		
		newlabels <- c(1:length(clusters_labels))
		
		temp <- sapply(data[,clusters],FUN=function(x){
					index <- which(x == clusters_labels)
					return(newlabels[index])
				})
		data[,clusters] <- temp
	}
	
	### Sort Data: Cluster (needed for rle)
	data <- data[order(data[,clusters]),]
				
	
	### Get Information
#	Time <- data[,time]
#	Status <- data[,status]
	Cluster <- data[,clusters]
#	Covariates <- data[,covariates,drop=FALSE]  # keep structure
	
	nclusterss <- length(levels(as.factor(Cluster)))
	
	ClusterData <- data.frame(Cluster=c(1:nclusterss),ClusterSubjects=rle(Cluster)$lengths,ClusterEvents=aggregate(data[,status],by=list(Cluster),FUN=sum)[,2])

	ClusterDataList <- vector("list",length=nclusterss)
	
	for(i.list in 1:nclusterss){
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
				data=data,time=time,status=status,clusters=clusters,covariates=covariates,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
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
		
		out4 <- res1_weibCL$value
	}
	
	if(marginal=="PiecewiseExp" & copula=="Clayton"){
		
		# Why different method here? Should this be decidable for the user?
		res1_pweCL <- optim(init.values,
				loglik.1stage_pweCL,cutpoints=cutpoints,num_pieces=num_pieces,data=data,status=status,time=time,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
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
		
		out4 <- res1_pweCL$value
	}
	
	
	if(marginal=="Weibull" & copula=="GH"){
		
#		res1_weibGH <- optim(init.values,
#				loglik.1stage_weibGH,data=data,time=time,status=status,clusters=clusters,covariates=covariates,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
#				hessian=TRUE,control=list(maxit=3000)) 
		
		res1_weibGH <- optim(init.values,
				loglik.1stage_GH,data=data,cutpoints=NULL,num_pieces=NULL,time=time,status=status,clusters=clusters,covariates=covariates,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
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
		
		out4 <- res1_weibGH$value
		
	}
	
	if(marginal=="PiecewiseExp" & copula=="GH"){
#		res1_pweGH <- optim(init.values,
#				loglik.1stage_pweGH,,cutpoints=cutpoints,num_pieces=num_pieces,data=data,status=status,time=time,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
#				hessian=TRUE,method="BFGS")

		res1_pweGH <- optim(init.values,
				loglik.1stage_GH,,cutpoints=cutpoints,num_pieces=num_pieces,data=data,status=status,time=time,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
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
		
		out4 <- res1_pweGH$value
	}
		
		# Give init.values correct naming + save all parameters
		names(init.values) <- rownames(out1)
		parameters <- list(time=time,status=status,clusters=clusters,covariates=covariates,init.values=init.values,
			marginal=marginal,copula=copula,n.piecewise=n.piecewise)
	
		out <- list(Parameters=out1,Kendall_Tau=out2,ParametersCov=out3,logllh=out4,parameter.call=parameters)
		class(out) <- "Sunclarco"
		return(out)

}



#init.values <- log(0.2)
# result <- CopulaModel_2stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
# 		init.values=init.values,
#		marginal="Weibull",copula="Clayton")
#result
#result <- CopulaModel_2stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
#		marginal="Weibull",copula="Clayton")
#result

#init.values <- c(rep(log(0.05),10),log(0.2),-0.05) # lambda's,theta, beta's,
#result <- CopulaModel_2stage(data=insem,time="Time",status="Status",clusters="Herd",covariates="Heifer",
#		init.values=init.values,
#		marginal="PiecewiseExp",copula="Clayton",n.piecewise=10)
#result

#' Title of this function (do later, function still temporary)
#' 
#' Description of this function. (Note: you can also add other fields like references, authors,... See which ones are necessary)
#' 
#' @export
#' @param data Input dataframe containing all variables.
#' @param time Which variable name is the time covariate?
#' @param status Which variable name is the status covariate? (more explanation?)
#' @param clusters Which variable name is the cluster covariate?
#' @param covariates Vector of 1 or more other covariates to be included in the model. (NOTE: MULTIPLE COVARIATES NEED TO BE TESTED!!!)
#' @param init.values Initial values for parameters: 
#' \itemize{
#' \item Clayton & Weibull: \deqn{\theta}
#' \item Clayton & Piecewise Exponential: \deqn{\lambda 's, \theta, \beta 's}
#' }
#' @param marginal Which marginal to use? Can either be \code{"Weibull"} or \code{"PiecewiseExp"} for Piecewise Exponential.
#' @param copula Which copula to user? Can either be \code{"Clayton"} or \code{"GH"} for Gumbel-Hougaard.
#' @param n.piecewise For \code{marginal="PiecewiseExp"}, how many pieces should the Piecewise Exponential have? (Default = 20)
#' @return Not yet final, so far: list object with: parameters (est & se) + kendal + cov matrix. Add more information? MAYBE ADD CHOSEN PARAMETERS (CHOSEN MODEL)! (NOTE: EXAMPLES CAN LATER BE RUN)
#' @examples
#' \dontrun{
#' data("insem",package="Sunclarco")
#' result <- CopulaModel_2stage()
#' }
CopulaModel_2stage <- function(data,time,status,clusters,covariates,init.values=NULL,marginal="Weibull",copula="Clayton",n.piecewise=20){
	
	
	## DATA CHECK AND INPUT SAME FOR COPULAMODEL 1stage 2stage
	# - > when merging in bigger function, move this part to there!!
	

	### Check Dataframe
	
	# TO DO: check if data is data.frame + time,status,clusters, covariates are characters -> Include in S4 method
	
	for(i in c(time,status,clusters,covariates)){
		if(!(i %in% colnames(data))){stop(paste0("\"",i,"\" is not a variable of the data frame."),call.=FALSE)}
	}
	
	## Automatic initial parameters + naming parameters correctly
	if(is.null(init.values)){
		init.values <- estimate_parameters(data=data,time=time,status=status,clusters=clusters,covariates=covariates,n.piecewise=n.piecewise,marginal=marginal,copula=copula,stage=2)
	}
		
	
	### Check parameters  (check if combination possible, length of init.values (lambda's and betas!) )
	if(!(marginal %in% c("Weibull","PiecewiseExp","Cox"))){stop(paste0("Parameter 'marginal' can not be ",marginal,". It should be either \"Weibull\" or \"PiecewiseExp\""),call.=FALSE)}
	if(!(copula %in% c("Clayton","GH"))){stop(paste0("Parameter 'copula' can not be ",copula,". It should be either \"Clayton\" or \"GH\""),call.=FALSE)}
	
	if(marginal=="Weibull"){
		correct_length <- 1
	}
	else if(marginal=="PiecewiseExp"){
		correct_length <- n.piecewise+1+length(covariates)
	}
	else if(marginal=="Cox"){
		correct_length <- 1
	}
	if(length(init.values)!=correct_length){stop(paste0("Parameter 'init.values' has an incorrect length. With given parameters it should be of length ",correct_length),call.=FALSE)}
	
	
	### CLUSTER VARIABLE NEEDS TO BE CHECKED! IT SHOULD BE A NUMBER FROM 1 TO N. IF NEEDED IT HAS TO BE TRANSFORMED!
	clusters_labels <- sort(unique(data[,clusters]))
	if(!identical(clusters_labels,c(1:length(clusters_labels)))){
		
		newlabels <- c(1:length(clusters_labels))
		
		temp <- sapply(data[,clusters],FUN=function(x){
					index <- which(x == clusters_labels)
					return(newlabels[index])
				})
		data[,clusters] <- temp
	}
	
	### Sort Data: Cluster (needed for rle)
	data <- data[order(data[,clusters]),]
	
	
	### Get Information
	Cluster <- data[,clusters]
	
	nclusterss <- length(levels(as.factor(Cluster)))
	
	ClusterData <- data.frame(Cluster=c(1:nclusterss),ClusterSubjects=rle(Cluster)$lengths,ClusterEvents=aggregate(data[,status],by=list(Cluster),FUN=sum)[,2])
	
	ClusterDataList <- vector("list",length=nclusterss)
	
	for(i.list in 1:nclusterss){
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
	
	if(marginal=="Weibull" & copula=="Clayton"){
		
		temp_formula <- Surv(data[,time],data[,status]) ~ cluster(data[,clusters])
		for(i.covariate in covariates){
			temp_formula <- update(temp_formula, ~ . + data[,i.covariate])
		}
		
		Surv_result   <- survreg(temp_formula,dist="weibull",robust=TRUE)
		
		mu <- Surv_result$coeff[1]
		gammas <- Surv_result$coeff[2:length(Surv_result$coeff)]
		sigma <- Surv_result$scale
		
		lambda2_weib <- as.numeric(exp(-mu/sigma)) 
		rho2_weib    <- as.numeric(1/sigma) 
		betas2_weib   <- as.numeric(-gammas/sigma) 
		
		
		
#		Surv_result$var
		#Surv_result$var gives ROBUST variance-covariance matrix for coef[1],coef[2] and log(scale)
		#use delta method to transform log(scale) to scale
		varmu <-Surv_result$var[1,1]
		
		vargammas <- 	matrix(Surv_result$var[c(2:length(Surv_result$coeff)),c(2:length(Surv_result$coeff))],nrow=length(gammas),ncol=length(gammas))
		
		covmugammas <- sapply(c(2:length(Surv_result$coeff)),FUN=function(x){Surv_result$var[1,x]})
		
		covmusigma <- Surv_result$var[1,dim(Surv_result$var)[2]]*Surv_result$scale
	
		covgammassigma <- sapply(c(2:length(Surv_result$coeff)),FUN=function(x){Surv_result$var[x,dim(Surv_result$var)[2]]*Surv_result$scale})
		
		varsigma <- Surv_result$var[dim(Surv_result$var)[2],dim(Surv_result$var)[2]]*Surv_result$scale^2
		
		#Variance-covariance matrix for marginal parameters lambda, rho, beta
		#Klein & Moeschberger (1997) p.378-379

		varbetas <- matrix(NA,nrow=length(gammas),ncol=length(gammas))
		for (j in 1:length(gammas)){
				for (k in 1:length(gammas)){
					varbetas[j,k] <- vargammas[j,k]/sigma^2-gammas[j]*covgammassigma[j]/sigma^3-gammas[k]*covgammassigma[k]/sigma^3+gammas[j]*gammas[k]*varsigma/sigma^4
				}
		}


		varlam=exp(-2*mu/sigma)*(varmu/sigma^2-2*mu*covmusigma/sigma^3+mu^2*varsigma/sigma^4)
		varrho=varsigma/sigma^4
		covbetaslam=exp(-mu/sigma)*(covmugammas/sigma^2-mu*covgammassigma/sigma^3-gammas*covmusigma/sigma^3+gammas*mu*varsigma/sigma^4)
		covbetasrho=covgammassigma/sigma^3-gammas*varsigma/sigma^4
		covlamrho=exp(-mu/sigma)*(covmusigma/sigma^3-mu*varsigma/sigma^4)
		
	
		# Covariance Matrix: Lambda, Rho, Betas

		VarCov2_lambda_rho <- matrix(c(
						varlam,covlamrho,covbetaslam,
						covlamrho,varrho,covbetasrho
									),nrow=2,ncol=(2+length(varbetas)),byrow=TRUE)
		
		VarCov2_betas_1 <- matrix(c(covbetaslam,covbetasrho),nrow=length(varbetas),ncol=2,byrow=FALSE)
		VarCov2_betas_2 <- varbetas
		
		VarCov2 <- rbind(VarCov2_lambda_rho,cbind(VarCov2_betas_1,VarCov2_betas_2))
		rownames(VarCov2) <- colnames(VarCov2) <- c("lambda","rho",paste0("beta_",covariates))
		
		IVI <- VarCov2
		
		stderrbetas2_weib <- sapply(c(3:dim(VarCov2)[2]),FUN=function(x){sqrt(VarCov2[x,x])})
		
		
		#vartheta2 <- 1/I22+I21%*%IVI%*%I12/(I22^2) #I12 en I22 komen nog uit one-stage 
		#stderrtheta2 <- sqrt(vartheta2)
		
		
		cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){betas2_weib %*% t(matrix(row,nrow=1,ncol=length(row)))})
		s <- exp(-lambda2_weib*exp(cov.lincomb)*data[,time]^rho2_weib) #u_ij
		f <- lambda2_weib*rho2_weib*data[,time]^(rho2_weib-1)*exp(cov.lincomb)*s #-du_ij/dy_ij
		
		
		for (i in 1:length(ClusterDataList)){
			ClusterDataList[[i]]$S <- s[data[,clusters]==i]
			ClusterDataList[[i]]$F <- f[data[,clusters]==i]
		}
		
		res2_weibCL <- optim(init.values,
				loglik.2stage_CL,
				ClusterData=ClusterData,ClusterDataList=ClusterDataList,copula="Clayton",status=status,
				hessian=TRUE,control=list(maxit=3000))
		
				
		theta2_weibCL <- exp(res2_weibCL$par) 
		
		#get estimates for I12 en I22 (not using complete one-stage results)
		#estimate for I22
		vartheta2stage_weibCL <- theta2_weibCL*solve(res2_weibCL$hessian)*theta2_weibCL #estimate for 1/I22
		stderrtheta2stage_weibCL <- sqrt(vartheta2stage_weibCL) #estimate for sqrt(1/I22)
		
		#estimate for I12
		#use hessian from one-stage loglikelihood, evaluated at 2-stage solution
		

		res12_weibCL <- optim(c(log(lambda2_weib),log(rho2_weib),log(theta2_weibCL),betas2_weib),
						loglik.1stage_weibCL,data=data,covariates=covariates,status=status,time=time,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
						hessian=TRUE,control=list(maxit=1))
		
		
		VarCov12_weibCL <- diag(c(lambda2_weib,rho2_weib,theta2_weibCL,rep(1,length(betas2_weib))))%*%solve(res12_weibCL$hessian)%*%diag(c(lambda2_weib,rho2_weib,theta2_weibCL,rep(1,length(betas2_weib))))
		
		
		II_weibCL.inv <- VarCov12_weibCL
		II_weibCL <- solve(VarCov12_weibCL)
		II_weibCL11 <- II_weibCL[-3,-3]
		II_weibCL21 <- II_weibCL[3,-3]
		II_weibCL12 <- II_weibCL[-3,3]
		II_weibCL22 <- II_weibCL[3,3]
		
		vartheta2_weibCL <- vartheta2stage_weibCL+II_weibCL21%*%IVI%*%II_weibCL12*vartheta2stage_weibCL^2 
		stderrtheta2_weibCL <- sqrt(vartheta2_weibCL)  
		
				
		#kendall's tau
		tau2_weibCL <- theta2_weibCL/(2+theta2_weibCL) 
		stderrtau2_weibCL <- (2/(theta2_weibCL+2)^2)*stderrtheta2_weibCL # OR stderrtheta2_weibCLa FROM COMPARING?!
		
		# VarCov2 holds the SE of the others
		stderr_all <- c(sqrt(diag(VarCov2)),stderrtheta2_weibCL)
		

		# Correct estimates and standard errors? No transformation required?
		# Which estimates do you take? From the first or from res12_weibCL ?
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambda2_weib,rho2_weib,betas2_weib,theta2_weibCL),StandardErrors=c(stderr_all))
		rownames(out1) <- c(rownames(VarCov2),"theta")

		out2 <- c(tau2_weibCL,stderrtau2_weibCL)
		names(out2) <- c("Estimates","StandardErrors")

		# DOES THETA VARIANCE HAVE TO BRE REPLACED?
		out3 <- VarCov2
#		rownames(out3) <- colnames(out3) <- rownames(out1)
		# WHICH VARIANCE MATRIX SHOULD BE OUTPUTTED?

		out4 <- res2_weibCL$value
		
		
#		Give init.values correct naming 
		names(init.values) <- "theta"
	}
	
	if(marginal=="PiecewiseExp" & copula=="Clayton"){

		
		# ADD JACKKNIFE THRESHOLD CHECK
		# TO DO

		
		# input lambda's, beta's	
		res2.pwe_stage1 <- optim(init.values[-(num_pieces+1)],loglik.1stage_pweCL,
				cutpoints=cutpoints,num_pieces=num_pieces,data=data,time=time,status=status,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,stage2part=TRUE,
				hessian=TRUE,method="BFGS")
	
		# (estimates)
		lambdas2_pwe <- exp(res2.pwe_stage1$par[1:num_pieces])
		betas2_pwe <- res2.pwe_stage1$par[(num_pieces+1):(length(res2.pwe_stage1$par))]

	
		#plug in into Clayton loglikelihood
		haz <- approx(cutpoints[1:num_pieces],lambdas2_pwe,xout=data[,time],method='constant',rule=2)$y
		DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
		Inter <- c(0,cumsum(lambdas2_pwe*DiffInter))
		cumhaz <- approx(cutpoints,Inter,xout=data[,time],method='linear',rule=2)$y
	
		cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){betas2_pwe %*% t(matrix(row,nrow=1,ncol=length(row)))})
			
		s <- exp(-cumhaz*exp(cov.lincomb)) #u_ij 
		f <- haz*exp(cov.lincomb)*s #-du_ij/dy_ij
	

		for(i in 1:length(ClusterDataList)){
			ClusterDataList[[i]]$S <- s[data[,clusters]==i] #contains estimated survival probabilities
			ClusterDataList[[i]]$F <- f[data[,clusters]==i] #contains estimated distribution
		}
		
	
		#plug marginal estimates into loglikelihood
	
		# theta input
		res2_pweCL <- optim(init.values[num_pieces+1],
				loglik.2stage_CL,status=status,ClusterData=ClusterData,ClusterDataList=ClusterDataList,copula="Clayton",
				hessian=TRUE,method="BFGS")
		theta2_pweCL <- exp(res2_pweCL$par) 
		
				
		jack_out <- jack_2stage_pweCL(data=data,covariates=covariates,status=status,time=time,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,num_pieces=num_pieces,init.values=init.values,cutpoints=cutpoints)
	
		betas2jack_pwe <- jack_out$betas
		theta2jack_pweCL <- jack_out$theta
		lambdas2jack_pwe <- jack_out$lambdas
		
#		stderrbetas2_pwe <- sqrt((length(ClusterDataList)-(num_pieces+1))/length(ClusterDataList)*sum((betas2jack_pwe-betas2_pwe)^2))  
		stderrbetas2_pwe <- apply(betas2jack_pwe,MARGIN=2,FUN=function(x){sqrt((length(ClusterDataList)-(num_pieces+length(covariates)))/length(ClusterDataList)*sum((x-betas2_pwe)^2))})  
		## More than 1 beta, still num_pieces+1
		
		stderrtheta2_pweCL <- sqrt((length(ClusterDataList)-(num_pieces+length(covariates)+1))/length(ClusterDataList)*sum((theta2jack_pweCL-theta2_pweCL)^2)) 
		# still num_pieces+2 (see above)
		
		# How about SE for lambdas?
		stderrlambdas2_pwe <- apply(lambdas2jack_pwe,MARGIN=2,FUN=function(x){sqrt((length(ClusterDataList)-(num_pieces+length(covariates)))/length(ClusterDataList)*sum((x-lambdas2_pwe)^2))})

		
		stderr_all <- c(stderrlambdas2_pwe,stderrtheta2_pweCL,stderrbetas2_pwe)
		
		#kendall's tau
		tau2_pweCL <- theta2_pweCL/(2+theta2_pweCL) 
		stderrtau2_pweCL <- (2/(theta2_pweCL+2)^2)*stderrtheta2_pweCL 
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambdas2_pwe,theta2_pweCL,betas2_pwe),StandardErrors=c(stderr_all))
		rownames(out1) <- c(paste0("lambda",c(1:n.piecewise)),"theta",paste0("beta_",covariates))
		
		out2 <- c(tau2_pweCL,stderrtau2_pweCL)
		names(out2) <- c("Estimates","StandardErrors")
		
		# SHOULD COVARIANCE BE HERE?
		out3 <- NULL
##		rownames(out3) <- colnames(out3) <- rownames(out1)
#		

		# CORRECT LIKELIHOOD???
		out4 <- res2_pweCL$value
		
#		Give init.values correct naming 
		names(init.values) <- rownames(out1)

	}
	
	
	if(marginal=="Cox" & copula=="Clayton"){
		
		
		temp_formula <- Surv(data[,time],data[,status]) ~ cluster(data[,clusters])
		for(i.covariate in covariates){
			temp_formula <- update(temp_formula, ~ . + data[,i.covariate])
		}
		
		# note: as.factor(heifer)
		
		PH_COX <- coxph(temp_formula)
		betas2_semipar <- as.numeric(PH_COX$coeff) 
		stderrbetas2_semipar <- sqrt(PH_COX$var) 
		
		
		# I AM HERE 
		
		
		fit <- survfit(PH_COX) # not needed
		fit0 <- survfit(PH_COX,newdata=data.frame(heifer=0)) # put all covariates at 0
		fit1 <- survfit(PH_COX,newdata=data.frame(heifer=1)) # not needed
		SS <- fit$surv
		S0 <- fit0$surv
		S1 <- fit1$surv
		
		S <- vector("list",length=K)
		
		for(i in 1:K){
			for(j in 1:length(T[[i]])){
				S[[i]][j] <- (S0[fit0$time==T[[i]][j]])^exp(beta2_semipar*HEIF[[i]][j])
			}
		}
		
		#plug marginal estimates into loglikelihood 
		
		loglik.semiparCL <- function(p,exp=FALSE){
			
			if (exp==TRUE){
				theta  <- exp(p)}
			
			else{ #default
				theta  <- p}
			
			G <- vector("list",K)
			H <- vector("list",K)
			sumG <- 1:K
			sumH <- 1:K 
			
			for (i in 1:K){
				G[[i]] <- C[[i]]*log(-1/varphi.prime(theta,varphi.inverse(theta,S[[i]])))
				H[[i]] <- varphi.inverse(theta,S[[i]])
				sumG[i] <- sum(G[[i]][C[[i]]==1])
				sumH[i] <- sum(H[[i]])}
			
			loglik <- sapply(d,function(x) ifelse(x==0,0,sum(log(1/theta+seq(0,x-1)))))-
					(1/theta)*log(theta)+ sumG + (-d-1/theta)*log(sumH+1/theta);
			-sum(loglik)
		}
		
		res2_semiparCL <- optim(c(log(0.5)),loglik.semiparCL,exp=TRUE,control=list(maxit=3000))
		theta2_semiparCL <- exp(res2_semiparCL$par)
		
		#jackknife procedure to get standard error
		
		
	}
	
	# save all parameters
	parameters <- list(time=time,status=status,clusters=clusters,covariates=covariates,init.values=init.values,
			marginal=marginal,copula=copula,n.piecewise=n.piecewise)
	
	out <- list(Parameters=out1,Kendall_Tau=out2,ParametersCov=out3,logllh=out4,parameter.call=parameters)
	class(out) <- "Sunclarco"
	return(out)
}


##################################################
##################################################



