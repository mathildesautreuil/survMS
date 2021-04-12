#' Get Parameters for weibull distribution
#'
#' @param int_a : grid of first parameter of the distribution
#' @param med : median of real dataset
#' @param mu : mean of real dataset
#'
#' @return : list with parameter values of the distribution
#' @import stats
#' @export
#'
#' @examples
#' library(survMS)
get_param_weib = function(int_a = c(0.1,11), med, mu){

  med_opt = function(x, med){
    a <- x
    abs(med - med_fct_avec_mu(a, mu))
  }

  res = optimise(med_opt, interval = int_a, med = med)

  a = res$minimum
  lambda = ((1/mu)*gamma(1+1/a))^(a)
  return(list(a = a, lambda = lambda))
}

med_fct_avec_mu = function(a, mu=2325){
 (log(2))^(1/a)*(mu/gamma(1+1/a))
}


#' Getting parameters of log-normal distribution
#'
#' @param int_a : grid of first parameter of the distribution
#' @param med : median of real dataset
#' @param mu : mean of real dataset
#'
#' @return : list with parameter values of the distribution
#' @import stats
#' @export
#'
#' @examples
#' library(survMS)
get_param_ln2 = function(int_a = c(0.1,11), med, mu){

  mu_opt = function(x, med, mu){
    a <- x
    abs(mu - exp(log(med) + a^2/2))
  }

  res = optimise(mu_opt, interval = int_a, med = med, mu = mu)

  a = res$minimum
  lambda = log(mu) - a^2/2
  return(list(a = a, lambda = lambda))
}

#' Getting parameters of log-normal distribution
#'
#' @param var : variance of real dataset
#' @param mu : mean of real dataset
#'
#' @return list with parameter values of the distribution
#' @export
#'
#' @examples
#' library(survMS)
get_param_ln = function(var = 170000, mu = 2325){
  a2 = log(1 + var/mu^2)
  lambda = log(mu) - (1/2)*a2
  a = sqrt(a2)

  return(list(a = a, lambda = lambda))
}


#' Getting parameters of exponential distribution
#'
#' @param int_a : grid of first parameter of the distribution
#' @param med : median of real dataset
#' @param mu : mean of real dataset
#'
#' @return : list with parameter values of the distribution
#' @import stats
#' @export
#'
#' @examples
#' library(survMS)
get_param_exp = function(int_a = c(0.000001,110), med, mu){

  mu_med_opt = function(x, med, mu){
    lambda <- x
    abs(mu - 1/lambda) + abs(med - log(2)/lambda)
  }

  res = optimise(mu_med_opt, interval = int_a, med = med, mu = mu)
  lambda = res$minimum
  return(list(lambda = lambda))
}


# SurvTimesAHWeib = function(){
#
# }
#
# SurvTimesAHLN = function(){
#
# }

