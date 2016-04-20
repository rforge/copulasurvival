# Project: Copula Package
# 
# Author: lucp8394
###############################################################################


estimate_parameters <- function(data,time,status,clusters,covariates,n.piecewise,marginal,copula,stage=1){

	theta <- ifelse(copula=="Clayton",0.5,0.55)
	
	if(stage==1){
		if(marginal=="Weibull"){
			#lambda, rho, beta's
						
			temp_formula <- Surv(data[,time],data[,status]) ~ cluster(data[,clusters])
			for(i.covariate in covariates){
				temp_formula <- update(temp_formula, ~ . + data[,i.covariate])
			}
			
			Surv_result   <- survreg(temp_formula,dist="weibull",robust=TRUE)
			
			mu <- Surv_result$coeff[1]
			gammas <- Surv_result$coeff[2:length(Surv_result$coeff)]
			sigma <- Surv_result$scale
			
			lambda <- as.numeric(exp(-mu/sigma)) 
			rho    <- as.numeric(1/sigma) 
			betas   <- as.numeric(-gammas/sigma) 
			
			out <- c(log(lambda),log(rho),log(theta),betas)
			cat("Initial Parameters:\n")
			cat("lambda =",out[1],"\n")
			cat("rho =",out[2],"\n")
			cat("theta =",out[3],"\n")
			cat("beta's =",paste(betas,collapse="; "),"\n")
		}
		
		if(marginal=="PiecewiseExp"){
		
#			lambda's, beta's
			temp_formula <- Surv(data[,time],data[,status]) ~ data[,covariates[1]]
			if(length(covariates)>1){
				for(i.covariate in covariates[-1]){
					temp_formula <- update(temp_formula, ~ . + data[,i.covariate])
				}
			}
		
			Surv_result   <- survreg(temp_formula,dist="weibull",scale=1)
					
			
			mu <- Surv_result$coeff[1]
			gammas <- Surv_result$coeff[2:length(Surv_result$coeff)]
			
			
			lambdas <- rep(as.numeric(exp(-mu)),n.piecewise)
			betas <- as.numeric(-gammas) 
			
		
			out <- c(log(lambdas),log(theta),betas)
			
			cat("Initial Parameters:\n")
			cat("lambda's =",paste(out[1],collapse="; "),"\n")
			cat("theta =",log(theta),"\n")
			cat("beta's =",paste(betas,collapse="; "),"\n")
		}
		
	}
	
	
	
	# Output
	if(stage==1){
#		print(out)
		return(out)
	}else if(stage==2){
		cat("Initial Parameters:\n")
		cat("theta = ")
		cat(theta)
		return(log(theta))
	}
	
	#also put a cat of the estimated parameters
}


#1-stage Weibull Clayton: Voor lambda, rho en beta's: gebruik survreg zoals in de eerste stap van de two-stage Weibull-Clayton (of Weibul-Gumbel: daar is de eerste stap ook hetzelfde aangezien dit het "Weibull-deel" betreft).

#		1-stage PWE Clayton: Voor lambda's en beta's: gebruik survreg als volgt:
#		Surv   <- survreg(Surv(time,stat)~heifer,dist="weibull",scale=1,data=insem)
#		mu <- Surv$coeff[1]
#		gamma <- Surv$coeff[2]
#		lambda_startvalue <- as.numeric(exp(-mu))
#		beta_startvalue <- as.numeric(-gamma) 
#		(de beta dien je nog uit te breiden naar een vector van beta's, zoals je bij de Weibull al hebt gedaan)

#1-stage Weibull-Gumbel: idem als 1-stage Weibull Clayton

#1-stage PWE Gumbel: idem als 1-stage PWE Clayton
#Voor de 2-stage procedures moet je enkel startwaarden meegeven voor de theta's, dus 0.5 en 0.55.
#		




