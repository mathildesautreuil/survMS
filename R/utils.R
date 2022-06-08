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


#' Heatmap of Covariate matrix
#'
#' @param x output of modelSim function (must be of type modSim)
#' @param k number of column split
#' @param ind indices of columns to keep
#' @param ... supplementary parameters
#'
#' @return Heatmap x
#' @import circlize
#' @importFrom ComplexHeatmap Heatmap
#' 
#' @examples
#' library(survMS)
Heatmap.modSim<- function(x, k, ind = NULL, ...){
  
  col_fun = colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
  if(!is.null(ind)){
    ind = ind
  }else{
    ind = which(x$betaNorm != 0)
  }
  rows_info <- rep("High", nrow(x$Z))
  rows_info[which(x$TC > median(x$TC))] <- "Low"
  colnames(x$Z) <- paste0("X", 1:ncol(x$Z))
  Heatmap(as.matrix(x$Z)[,ind], name = "expression", 
          row_split = rows_info, 
          # column_split = c(rep("Sign", sum(x$betaNorm[ind] != 0)),
          #                  rep("No Sign", 
          #                      ncol(x$Z[,ind]) - sum(x$betaNorm[ind] != 0))),
          col = col_fun, #row_km = 2,
          show_column_names = TRUE, column_km = k)
  
}
