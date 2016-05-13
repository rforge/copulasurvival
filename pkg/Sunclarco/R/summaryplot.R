# Project: Copula Package
# 
# Author: Gebruiker
###############################################################################


## ADD SUMMARY, PRINT AND PLOT

#' @export
summary.Sunclarco <- function(object,...){
	cat("Execution Time:",object$parameter.call$elapsedtime,"mins\n")
	cat("Copula:",object$parameter.call$copula,"\n")
	cat("Marginal Survival Distribution:",object$parameter.call$marginal,"\n")
	cat("Loglikelihood:",object$logllh,"\n\n")
	print(object$Parameters)
	cat("\n")
	temp <- matrix(object$Kendall_Tau,nrow=1,ncol=2)
	rownames(temp) <- "Kendall's Tau"
	colnames(temp) <- c("Estimate","StandardError")
	print(temp)
}


#' @export
print.Sunclarco <- function(x,...){
	summary(x)
}