#' Simulation survival times from Cox/Weibull model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvTimesCoxWeib = function(Z, beta, Y, pp, hazParams){

  Ts=(-(1/hazParams[2])*exp((1/sqrt(pp))*(-Z %*% beta))*log(1-Y))^(1/hazParams[1]) #(1/(p*sum(beta)))* (1/(p))* (1/(p*sum(beta)))* (1/(sqrt(5*p)))* (1/sqrt(10000*p))*
  return(Ts)
}

#' Internal heatmap object
#'
#' @param x modSim obkect
#' @param k number of clusters for the heatmap
#' @param ind number of column covariate matrix
#' @param ... supplementary arguments
#'
#' @return Heatmap x
#' @export
#'
#' @examples
#' library(survMS)
Heatmap <- function(x, k, ind = NULL, ...) {
  
  UseMethod("Heatmap")
}


#' Draw Heatmap of modSim object
#'
#' @param x modSim object
#' @param k number of split for heatmap's columns
#' @param ind number of columns to keep
#' @param ... supplementary parameters
#'
#' @return draw x
#' @export
#' @import circlize
#' @importFrom ComplexHeatmap draw
#' 
#' @examples
#' library(survMS)
#' res_paramW = get_param_weib(med = 2.5, mu = 1.2)
#' listCoxSimCor_n500_p1000 <- modelSim(model = "cox", matDistr = "mvnorm", 
#'                                      matParam = c(0,0.6), n = 500, 
#'                                      p = 1000, pnonull = 20, betaDistr = 1, 
#'                                      hazDistr = "weibull", 
#'                                      hazParams = c(res_paramW$a, res_paramW$lambda), 
#'                                      seed = 1, d = 0)
#' print(listCoxSimCor_n500_p1000)
#' hist(listCoxSimCor_n500_p1000)
#' #draw(listCoxSimCor_n500_p1000, k = 3)
draw.modSim<- function(x, k, ind = NULL, ...){
  
  col_fun = colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
  if(!is.null(ind)){
    ind_col = ind
  }else{
    ind_col = which(x$betaNorm != 0)
  }
  rows_info <- rep("High", nrow(x$Z))
  rows_info[which(x$TC > median(x$TC))] <- "Low"
  colnames(x$Z) <- paste0("X", 1:ncol(x$Z))
  ht <- Heatmap(as.matrix(x$Z)[,ind_col], name = "expression", 
                row_split = rows_info, 
                # column_split = c(rep("Sign", sum(x$betaNorm[ind] != 0)),
                #                  rep("No Sign", 
                #                      ncol(x$Z[,ind]) - sum(x$betaNorm[ind] != 0))),
                col = col_fun, #row_km = 2,
                show_column_names = TRUE, column_km = k)
  draw(ht, heatmap_legend_side = "bottom")
  
}


#' Simulation survival times from Cox/Log-normal model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvTimesCoxLN = function(Z, beta, Y, pp, hazParams){

  Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
  # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
  return(Ts)
}

#' Simulation survival times from Cox/exponential model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvTimesCoxExp = function(Z, beta, Y, pp, hazParams){

  Ts <- (-log(1-Y) / (hazParams[1] * exp((1/(pp))*Z%*%beta)))
  # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
  # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
  return(Ts)
}

#' Simulation survival times from Cox/gompertz model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvTimesCoxGomp = function(Z, beta, Y, pp, hazParams){

  check <- ((-hazParams[1]*log(1-Y)) / (hazParams[2]*exp((1/(pp))*Z%*%beta))) + 1
  if (check < 0)
    return(Inf)
  Ts <- (1 / hazParams[1]) * log(check)
  # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
  # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
  return(Ts)
}


#' Simulation survival times from AH/Weibull model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvTimesAHWeib = function(Z, beta, Y, pp, hazParams){


  eta = exp((1/sqrt(pp))*(Z %*% beta))
  Ts = (1/eta)*(((-log(1-Y)*eta)/hazParams[2])^(1/hazParams[1]))
  # Ts=(-(1/hazParams[2])*exp((1/sqrt(pp))*(-Z %*% beta))*log(1-Y))^(1/hazParams[1]) #(1/(p*sum(beta)))* (1/(p))* (1/(p*sum(beta)))* (1/(sqrt(5*p)))* (1/sqrt(10000*p))*
  # stop("log-normale distribution must be used")
  return(Ts)

}

#' Simulation survival times from AH/Log-normal model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return Ts Observed times
#' @export
#'
#' @examples
#' library(survMS)
SurvTimesAHLN = function(Z, beta, Y, pp, hazParams){

  eta = exp((1/sqrt(pp))*(Z %*% beta))
  Ts<- (1/eta)*(exp(hazParams[1]*qnorm(1-exp((log(1-Y)*eta)),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2]) ) #+r*eta2
  return(Ts)

}
#' Simulation survival times from AFT/Weibull model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvTimesAFTWeib = function(Z, beta, Y, pp, hazParams){

  # phi2 = 0

  # b2 =  c(runif(pp, -1.5, 1.5), rep(0, ncol(Z)-pp))#runif(ncol(Z), -1.5, 1.5)
  # phi2 = -Z %*% b2
  # phi2 = 200*(1/sqrt(pp))*(-Z %*% b2)
  phi2 = 0
  eta = exp((1/sqrt(pp))*(Z %*% beta))
  Ts = (1/eta)*((-log(1-Y)/hazParams[2])^(1/hazParams[1])-phi2)
  # Ts=(-(1/hazParams[2])*exp((1/sqrt(pp))*(-Z %*% beta))*log(1-Y))^(1/hazParams[1]) #(1/(p*sum(beta)))* (1/(p))* (1/(p*sum(beta)))* (1/(sqrt(5*p)))* (1/sqrt(10000*p))*
  # stop("log-normale distribution must be used")
  return(Ts)

}

#' Simulation survival times from AFT/Log-normal model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvTimesAFTLN = function(Z, beta, Y, pp, hazParams){

  phi2 = 0
  eta = exp((1/sqrt(pp))*Z%*%beta)#
  Ts<- (1/eta)* (exp(hazParams[1]*qnorm(1-Y,mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])-phi2) #runif(n,0,1) ((1/(p))* ((1/(p))* ((1/(p))*Z%*%b0) (1/(sqrt(p))* (1/p)*
  return(Ts)

}

#' Simulation survival times from Shift AFT/Weibull model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#' @param beta2 vector of regression parameter or distribution of regression parameter
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords inetrnal
#'
#' @examples
#' library(survMS)
SurvTimesAFTshiftWeib = function(Z, beta, beta2, Y, pp, hazParams){


  b2 = beta2#c(runif(pp, -1.5, 1.5), rep(0, ncol(Z)-pp))
  # phi2 = -Z %*% b2#300* (1/sqrt(pp))*(-Z %*% beta)
  phi2 = 500* (1/(sqrt(pp))) * Z%*%b2 #200*(-Z %*% b2) #200*(1/sqrt(pp))*
  eta = exp((1/sqrt(pp))*(Z %*% beta))
  Ts = (1/eta)*((-log(1-Y)/hazParams[2])^(1/hazParams[1])+phi2)
  # Ts=(-(1/hazParams[2])*exp((1/sqrt(pp))*(-Z %*% beta))*log(1-Y))^(1/hazParams[1]) #(1/(p*sum(beta)))* (1/(p))* (1/(p*sum(beta)))* (1/(sqrt(5*p)))* (1/sqrt(10000*p))*
  # stop("log-normale distribution must be used")
  return(Ts)

}

#' Simulation survival times from Shift AFT/Log-normal model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param hazParams distribution parameters of baseline hazard risk
#' @param beta2 vector of regression parameter or distribution of regression parameter
#'
#' @return Ts Observed times
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvTimesAFTshiftLN = function(Z, beta, beta2, Y, pp, hazParams){

  b2 =  beta2#, c(runif(pp, -1.5, 1.5), rep(0, ncol(Z)-pp))#runif(ncol(Z), -1.5, 1.5)
  # phi2 = -Z %*% b2
  # phi2 = 200*(-Z %*% b2) #200*(1/sqrt(pp))*
  phi2 = 500* (1/(sqrt(pp))) * Z%*%b2
  # phi2 = 300* (1/sqrt(pp))*(-Z %*% beta)
  eta = exp((1/sqrt(pp))*Z%*%beta)#
  Ts<- (1/eta)* (exp(hazParams[1]*qnorm(1-Y,mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])+phi2) #runif(n,0,1) ((1/(p))* ((1/(p))* ((1/(p))*Z%*%b0) (1/(sqrt(p))* (1/p)*
  return(Ts)

}

# SurvTimesAHWeib = function(){
#
# }
#
# SurvTimesAHLN = function(){
#
# }

#' Survival curves of simulated data with Cox/Weibull model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return SurvFctCoxWeib returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctCoxWeib = function(Z, beta, pp, Ts, hazParams){

  tau = max(Ts)
  pas=100
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  h0_t = hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))
  h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(1/eta_i)
  H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
  F_t = 1 - exp(-H0_t*as.vector(eta_i))
  S_t = exp(-H0_t*as.vector(eta_i))

  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = h, grillet = grille_ti, tau))
}

#' Survival curves of simulated data with Cox/Log-normal model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return SurvFctCoxWeib returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctCoxLN = function(Z, beta, pp, Ts, hazParams){

  tau = max(Ts)
  pas=100
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  mat_grille_ti = matrix(grille_ti, nrow = nrow(Z), ncol = pas, byrow = T)
  h0_t = ((1/hazParams[2])*dlnorm(grille_ti, meanlog = hazParams[1], sdlog = hazParams[2]))/(1-plnorm(grille_ti, meanlog = hazParams[1], sdlog = hazParams[2]))
  h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(1/eta_i)
  H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
  F_t = 1 - exp(-H0_t*as.vector(1/eta_i))
  S_t = exp(-H0_t*as.vector(1/eta_i))

  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = h, grillet = grille_ti, tau))
}


#' Survival curves of simulated data with Cox/Gompertz model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return SurvFctCoxGomp returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctCoxGomp = function(Z, beta, pp, Ts, hazParams){

  #lambda = hazParams[1]
  #a = hazParams[2]
  tau = max(Ts)
  pas=100
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  h0_t = hazParams[1]*exp(hazParams[2]*grille_ti)
  h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
  H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
  F_t = 1 - exp(-H0_t*as.vector(1/eta_i))
  S_t = exp(-H0_t*as.vector(1/eta_i))
  # stop("to code")

  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = h, grillet = grille_ti, tau))
}

#' Survival curves of simulated data with Cox/Exponential model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return SurvFctCoxExp returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctCoxExp = function(Z, beta, pp, Ts, hazParams){

  tau = max(Ts)
  pas=100
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  h0_t = hazParams[1]
  h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(1/eta_i)
  H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
  F_t = 1 - exp(-H0_t*as.vector(1/eta_i))
  S_t = exp(-H0_t*as.vector(1/eta_i))
  # stop("to code")

  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = h, grillet = grille_ti, tau))
}

#' Survival curves of simulated data with AFT/Log-normal model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return SurvFctCoxWeib returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctAFTLN = function(Z, beta, pp, Ts, hazParams){
  tau = max(Ts)
  pas=100
  # a = hazParams[1]
  # lambda = hazParams[2]
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  # h0_t = ((1/hazParams[2])*dlnorm(grille_ti, meanlog = hazParams[1], sdlog = hazParams[2]))/(1-plnorm(grille_ti, meanlog = hazParams[1], sdlog = hazParams[2]))
  # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
  # H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
  # F_t = 1 - exp(-H0_t*as.vector(eta_i))
  # S_t = exp(-H0_t*as.vector(eta_i))
  #
  # H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
  # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
  # H0_t = matrix(-log(1 - pnorm((log(grille_ti*eta_i)-hazParams[2])/hazParams[1])),
  #               nrow = nrow(Z), ncol = length(grille_ti))
  mat_grille_ti = matrix(grille_ti, nrow = nrow(Z), ncol = pas, byrow = T)
  # (1/(hazParams[1]*((1/as.vector(eta_i))*mat_grille_ti)))*
  h0_t = (dlnorm(mat_grille_ti*as.vector(1/eta_i),
                                  meanlog = hazParams[2],
                                  sdlog = hazParams[1]))/(pnorm(-(log(mat_grille_ti*(1/as.vector(eta_i)))-hazParams[2])/hazParams[1]))

    # (1-plnorm(mat_grille_ti*as.vector(1/eta_i),
    #                                                                meanlog = hazParams[2],
    #                                                                sdlog = hazParams[1]))#hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))
  ht = h0_t*(1/as.vector(eta_i))
  ht[which(is.na(ht))]=0
  # ht[which(is.infinite(ht))]=0
  H0_t = matrix(-log(1 - plnorm(mat_grille_ti*as.vector(1/eta_i),meanlog = hazParams[2], sdlog = hazParams[1])),
                nrow = nrow(Z), ncol = length(grille_ti))
  H0_t[which(is.na(H0_t))]=0
  # H0_t[which(is.infinite(H0_t))]=0
  F_t = 1 - exp(-H0_t)
  S_t = exp(-H0_t)

  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = ht, grillet = grille_ti, tau))
}

#' Survival curves of simulated data with AFT/Weibull model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return SurvFctAFTWeib returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctAFTWeib = function(Z, beta, pp, Ts, hazParams){
  tau = max(Ts)
  pas=100
  # a = hazParams[1]
  # lambda = hazParams[2]
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  mat_grille_ti = matrix(grille_ti, nrow = nrow(Z), ncol = pas, byrow = T)
  h0_t = (1/as.vector(eta_i))*hazParams[1]*hazParams[2]*((mat_grille_ti*as.vector(1/eta_i))^(hazParams[1]-1))
  ht = h0_t*(1/as.vector(eta_i))
  ## to modify
  # h = NULL
  # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
  # H0_t = -log(1 - plnorm(mat_grille_ti*as.vector(eta_i),meanlog = hazParams[1], sdlog = hazParams[2]))
  H0_t = hazParams[2]*(mat_grille_ti*as.vector(1/eta_i))^hazParams[1]
  F_t = 1 - exp(-H0_t)
  S_t = exp(-H0_t)

  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = ht, grillet = grille_ti, tau))
}

#' Survival curves of simulated data with Shifted AFT/Weibull model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#' @param beta2 vector of regression parameter or distribution of regression parameter
#'
#' @return SurvFctCoxWeib returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctAFTshiftWeib = function(Z, beta, beta2, pp, Ts, hazParams){
  tau = max(Ts)
  pas=100
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  # phi2 = 300*eta_i
  b2 =  beta2 #c(runif(pp, -1.5, 1.5), rep(0, ncol(Z)-pp))#runif(ncol(Z), -1.5, 1.5)
  # phi2 = -Z %*% b2
  # phi2 = 200*(-Z %*% b2) #200*(1/sqrt(pp))*
  phi2 = 500* (1/(sqrt(pp))) * Z%*%b2
  # phi2 = 300* (1/sqrt(pp))*(-Z %*% beta)
  # a = hazParams[1]
  # lambda = hazParams[2]
  grille_ti=tau*(1/pas)*c(1:(pas))
  # h0_t = hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))
  # ## to modify
  # h = NULL
  # # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
  # H0_t = matrix(hazParams[2]*(grille_ti*eta_i)^hazParams[1] + phi2,
  #               nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
  mat_grille_ti = matrix(grille_ti, nrow = nrow(Z), ncol = pas, byrow = T)
  h0_t = (1/as.vector(eta_i))*hazParams[1]*hazParams[2]*((mat_grille_ti*as.vector(1/eta_i) - matrix(phi2, nrow(Z), ncol = pas))^(hazParams[1]-1))
  ## to modify
  # h = NULL
  # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i
  # H0_t = -log(1 - plnorm(mat_grille_ti*as.vector(eta_i),meanlog = hazParams[1], sdlog = hazParams[2]))
  ht = h0_t*(1/as.vector(eta_i))

  H0_t = hazParams[2]*(mat_grille_ti*as.vector(1/eta_i)-matrix(phi2, nrow(Z), ncol = pas))^hazParams[1]

  F_t = 1 - exp(-H0_t)
  S_t = exp(-H0_t)

  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = ht, grillet = grille_ti, tau))
}

#' Survival curves of simulated data with Shifted AFT/Log-normal model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#' @param beta2 vector of regression parameter or distribution of regression parameter
#'
#' @return SurvFctCoxWeib returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctAFTshiftLN = function(Z, beta, beta2, pp, Ts, hazParams){

  tau = max(Ts)
  pas=100
  # a = hazParams[1]
  # lambda = hazParams[2]
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  # phi2 = 300*eta_i
  b2 =  beta2 #c(runif(pp, -1.5, 1.5), rep(0, ncol(Z)-pp))#runif(ncol(Z), -1.5, 1.5)
  # phi2 = 200*(-Z %*% b2) #200*(1/sqrt(pp))*
  phi2 = 500* (1/(sqrt(pp))) * Z%*%b2
  # phi2 = 300* (1/sqrt(pp))*(-Z %*% beta)
  # h0_t = ((1/hazParams[2])*dlnorm(grille_ti, meanlog = hazParams[1], sdlog = hazParams[2]))/(1-plnorm(grille_ti, meanlog = hazParams[1], sdlog = hazParams[2]))
  # # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
  # h = NULL
  # H0_t = matrix(-log(1 - pnorm((log(grille_ti*eta_i)-hazParams[2])/hazParams[1])) + phi2,
  #               nrow = nrow(Z), ncol = length(grille_ti))
  mat_grille_ti = matrix(grille_ti, nrow = nrow(Z), ncol = pas, byrow = T)
  # (1/(hazParams[1]*((1/as.vector(eta_i))*mat_grille_ti- matrix(phi2, nrow(Z), ncol = pas))))*
  h0_t = ((dlnorm(mat_grille_ti*as.vector(1/eta_i) - matrix(phi2, nrow(Z), ncol = pas),
                                  meanlog = hazParams[2],
                                  sdlog = hazParams[1]))/pnorm(-(log(mat_grille_ti*(1/as.vector(eta_i))- matrix(phi2, nrow(Z), ncol = pas))-hazParams[2])/hazParams[1]))
    # (1-plnorm(mat_grille_ti*as.vector(1/eta_i) - matrix(phi2, nrow(Z), ncol = pas),
    #                                                                meanlog = hazParams[2],
    #                                                                sdlog = hazParams[1]))#hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))

  ht = h0_t*(1/as.vector(eta_i))
  ht[which(is.na(ht))]=0
  # ht[which(is.infinite(ht))]=0
  H0_t = matrix(-log(1 - plnorm(mat_grille_ti*as.vector(1/eta_i) - matrix(phi2, nrow(Z), ncol = pas), meanlog = hazParams[2], sdlog = hazParams[1])),
                nrow = nrow(Z), ncol = length(grille_ti))
  H0_t[which(is.na(H0_t))]=0
  # H0_t[which(is.infinite(H0_t))]=0
  F_t = 1 - exp(-H0_t)
  S_t = exp(-H0_t)

  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = ht, grillet = grille_ti, tau))
}


#' Survival curves of simulated data with AH/Weibull model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return SurvFctCoxWeib returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctAHWeib = function(Z, beta, pp, Ts, hazParams){
  tau = max(Ts)
  pas=100
  # a = hazParams[1]
  # lambda = hazParams[2]
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  mat_grille_ti = matrix(grille_ti, nrow = nrow(Z), ncol = pas, byrow = T)
  h0_t = hazParams[1]*hazParams[2]*((mat_grille_ti*as.vector(1/eta_i))^(hazParams[1]-1))
  ## to modify
  # h = NULL
  # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i
  # H0_t = -log(1 - plnorm(mat_grille_ti*as.vector(eta_i),meanlog = hazParams[1], sdlog = hazParams[2]))
  H0_t = hazParams[2]*(mat_grille_ti*as.vector(1/eta_i))^hazParams[1]
  F_t = 1 - exp(-H0_t*as.vector(eta_i))
  S_t = exp(-H0_t*as.vector(eta_i))
  # S_t <- F_t <- H0_t <- h <- grille_ti <- NULL
  # stop("to code")
  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = h0_t, grillet = grille_ti, tau))
}

#' Survival curves of simulated data with AH/Log-normal model
#'
#' @param Z Matrix of covariates
#' @param beta regression parameter
#' @param Y random uniform
#' @param pp number of pertinent covariates
#' @param Ts observed times
#' @param hazParams distribution parameters of baseline hazard risk
#'
#' @return SurvFctCoxWeib returns a list containing: \itemize{
#' \item{St}{ Matrix of survival functions (rows: individuals, columns: time grid)}
#' \item{Ft}{ Matrix of cumulative functions (rows: individuals, columns: time grid)}
#' \item{H0t}{ Matrix of cumulative hazard functions (rows: individuals, columns: time grid)}
#' \item{ht}{ Matrix of hazard risk functions (rows: individuals, columns: time grid)}
#' \item{grillet}{ Time grid}
#' \item{tau}{ maximum of observed times}
#' }
#' @export
#'
#' @keywords internal
#'
#' @examples
#' library(survMS)
SurvFctAHLN = function(Z, beta, pp, Ts, hazParams){
  tau = max(Ts)
  pas=100
  # a = hazParams[1]
  # lambda = hazParams[2]
  grille_ti=tau*(1/pas)*c(1:(pas))
  eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
  mat_grille_ti = matrix(grille_ti, nrow = nrow(Z), ncol = pas, byrow = T)
  # (((1/(alpha*grille_ti*eta_i[i])) * dnorm((log(grille_ti*eta_i[i])-lambda)/alpha)) / (pnorm(-(log(grille_ti*eta_i[i])-lambda)/alpha)))

  h0_t = ((1/(hazParams[1]))*as.vector(eta_i)*dlnorm(mat_grille_ti*as.vector(1/eta_i),
                                  meanlog = hazParams[2],
                                  sdlog = hazParams[1]))/(1-plnorm(mat_grille_ti*as.vector(1/eta_i),
                                                                   meanlog = hazParams[2],
                                                                   sdlog = hazParams[1]))#hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))
  H0_t = matrix(-log(1 - plnorm(mat_grille_ti*as.vector(1/eta_i), meanlog = hazParams[2], sdlog = hazParams[1])),
                nrow = nrow(Z), ncol = length(grille_ti))
  F_t = 1 - exp(-H0_t*as.vector(eta_i))
  S_t = exp(-H0_t*as.vector(eta_i))
  # S_t <- F_t <- H0_t <- h <- grille_ti <- NULL
  # stop("to code")
  return(list(St = S_t, Ft = F_t, H0t = H0_t, ht = h0_t, grillet = grille_ti, tau))
}
