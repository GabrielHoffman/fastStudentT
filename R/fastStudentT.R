
#' Estimate parameters of Student t distribution
#' 
#' Estimation is substantially faster for large datasets than other methods 
#' 
#' @param y vector of observed data 
#' @param df degrees of freedom used to initialize the algorithm
#' @param maxiter maximum number of iterations
#' @param tol convergence tolerance for the log likelihood
#' @param verbose show messages
#'
#' @import Rcpp
#' @importFrom stats optim
#' @useDynLib fastStudentT
#'
#' @export
fastStudentT = function(y, df=50, maxiter=1000, tol = 1e-3, verbose=FALSE){

	nll = function(df){
		# fast evaluation of Student T log likelihood
		-1*dt_fast(y, fit$m, fit$s, df)
	}

	converged = FALSE
	logLikValue = Inf
	for(i in 1:maxiter){
		# Estimate mu and s given df using EM algorithm
		fit = studentTEM(y, df, 100)

		# Estimate df given mu and s using numerical optimization
		fit2 = optim(df, nll, lower=1, upper=1e5 ,method= "L-BFGS-B")
		ll = fit2$value
		df = fit2$par

		delta = logLikValue - ll
		if( verbose & i > 1) message("iter", i, ": ", delta)
		if( delta < tol ){
			converged = TRUE
			break
		}
		logLikValue = ll
	}

	if( ! converged ) warning("Optimization did not converge")

	data.frame( m = fit$m, s = fit$s, df = df, converged=converged)
}



