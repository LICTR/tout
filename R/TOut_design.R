#' Find optimal sample size and progression criteria
#' 
#' Given a null and alternative hypothesis, this function finds the 
#' lowest sample size such that a design with optimal progression criteria (as 
#' determined by the function `opt_pc`) satisfies upper constraints on three
#' operating characteristics.
#'
#' @param rho_0 null hypothesis.
#' @param rho_1 alternative hypothesis.
#' @param alpha_nom nominal upper constraint on alpha.
#' @param beta_nom nominal upper constraint on beta.
#' @param gamma_nom nominal upper constraint on gamma. Defaults to 1.
#' @param eta_0 probability of an incorrect decision under the null hypothesis
#' after an intermediate result. Defaults to 0.5.
#' @param eta_1 probability of an incorrect decision under the alternative hypothesis
#' after an intermediate result. Defaults to eta_0.
#' @param tau two element vector denoting lower and upper limits of the 
#' effect of adjustment.
#' @param max_n optional upper limit to use in search over sample sizes.
#' @param n optional sample size (optimised if left unspecified).
#' @param x optional vector of decision thresholds (optimised if left unspecified).
#' @param sigma standard deviation of outcome. If left unspecified, a binary outcome is assumed.
#' 
#' @return An object of class `tout`, which is a list containing the following components:
#' 
#' \item{`valid`}{boolean indicating if the nominal constraints are met.}
#' \item{`n`}{sample size.}
#' \item{`thesholds`}{numeric vector of the two decision thresholds.}
#' \item{`alpha`}{attained value of operating characteristic alpha.} 
#' \item{`beta`}{attained value of operating characteristic beta.} 
#' \item{`gamma`}{attained value of operating characteristic gamma.} 
#'
#' @examples
#' rho_0 <- 0.5
#' rho_1 <- 0.7
#' alpha_nom <- 0.05
#' beta_nom <- 0.2
#'
#' tout_design(rho_0, rho_1, alpha_nom, beta_nom)
#' 
#' # Allowing for adjustment effects:
#' 
#' tout_design(rho_0, rho_1, alpha_nom, beta_nom, tau = c(0.08, 0.12))
#' 
#' # Allowing for different error probabilities following a pause decision
#' 
#' tout_design(rho_0, rho_1, alpha_nom, beta_nom, eta_0 = 0.3)
#' 
#' # Designs for continuous outcomes:
#' 
#' tout_design(rho_0 = 0, rho_1 = 0.4, alpha_nom, beta_nom, sigma = 1)
#' 
#' @export
tout_design <-  function(rho_0, rho_1, alpha_nom, beta_nom, gamma_nom = 1, eta_0 = 0.5, eta_1 = eta_0, tau = c(0,0), max_n = NULL, n = NULL, x = NULL, sigma = NULL){

  validate_tout(new_tout(FALSE, n, NA, NA, NA, NA, NA, alpha_nom, beta_nom, rho_0, rho_1, tau, eta_0, eta_1, sigma))
  
  if(is.null(max_n)){
    # Get sample size for standard two outcome design taking a normal approx
    # and making conservative assumption on variance (maximised at rho = 0.5 for binary case )
    if(is.null(sigma)){ 
      n_two <- 0.25*(stats::qnorm(1 - alpha_nom) - stats::qnorm(beta_nom))^2/(rho_1 - rho_0)^2
    } else {
      n_two <- sigma^2*(stats::qnorm(1 - alpha_nom) - stats::qnorm(beta_nom))^2/(rho_1 - rho_0)^2
    }
    # Set max_n at an (arbitrarily) large multiple of this
    max_n <- floor(5*n_two)
  }
  
  if(!is.null(n)){
    valid <- TRUE
    final_design <- opt_pc(n, rho_0, rho_1, alpha_nom, beta_nom, 
                         tau = tau, eta_0 = eta_0, eta_1 = eta_1, x = x, sigma= sigma)
  } else {
    # An exhaustive search here due to the non-monotonic relationship with n
    n <- 0
    valid <- FALSE
    while(!valid & n <= max_n){
      n <- n + 1
      design <- opt_pc(n, rho_0, rho_1, alpha_nom, beta_nom, 
                            tau = tau, eta_0 = eta_0, eta_1 = eta_1, x = x, sigma= sigma)
      if(design$valid & (design$gamma <= gamma_nom)) {
        final_design <- design
        valid <- TRUE
      }
    }
  }
  
  if(valid){
    return(final_design)
  } else {
    warning("No valid design found. Consider increasing the maximum sample size (max_n).")
  }
}
