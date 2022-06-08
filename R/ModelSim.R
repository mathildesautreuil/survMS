#' Data simulation from different survival models
#'
#' @param model Survival model: "cox", "AFT", "AFTshift" or "AH"
#' @param matDistr Distribution of matrix
#' @param n size of sample
#' @param p number of parameters
#' @param pourc pourcents
#' @param seed seed
#' @param matParam Parameters of matrix
#' @param hazParams Parameters of baseline hazard
#' @param pnonull number of partinent covariates
#' @param betaDistr Distribution of beta or vector of beta
#' @param Phi nonlinearity (not coded)
#' @param d censorship
#' @param hazDistr distribution of baseline hazard
#'
#' @details This function simulates survival data from different models: Cox model, AFT model and AH model.
#' 1. The Cox model is defined as:
#' \eqn{
#' \lambda(t|X) = \alpha_0(t) \exp(\beta^T X_{i.}),
#' }
#' with \eqn{\alpha_0(t)} is the baseline risk and \eqn{\beta} is the vector of coefficients.
#' Two distributions are considered for the baseline risk:
#' \itemize{
#'   \item Weibull: \eqn{\alpha_0(t) = \lambda a t^{(a-1)}};
#'   \item Log-normal: \eqn{\alpha_0(t) = (1/(\sigma\sqrt(2\pi t) \exp[-(\log t - \mu)^2 /2 \sigma^2]))/(1 - \Phi[(\log t - \mu)/\sigma}]);
#'   \item Exponential: \eqn{\alpha_0(t) = \lambda};
#'   \item Gompertz: \eqn{\alpha_0(t) = \lambda \exp(\alpha t)}.
#'   }
#'   To Simulate the covariates, two distributions are also proposed:
#'     \itemize{
#'     \item Uniform
#'     \item Normal
#'     }
#'     and the choice of parameters
#'     The Phi parameter enables to simulate survival data in a linear framework with no interaction, but its future implementation will take into account a non-linear framework with interactions. If the parameter Phi is NULL (to complete...).
#'
#' 2. The AFT model is defined from a linear regression of the interest covariate:
#' \eqn{
#'  Y_i = X_{i.} \beta + W_i,
#' }
#' with \eqn{X_{i.}} the covariates, \eqn{\beta} the vector of regression coefficients et \eqn{\epsilon_i} the error term
#' AFT model can also be defined from the baseline survival function \eqn{S_0(t)}, corresponding distribution tail \eqn{\exp(\epsilon_i)}. Survival function of AFT model is written as:
#' \eqn{
#'  S(t|{X_{i.}}) = S_0(t\exp{(\beta^T X_{i.})}),
#'  }
#'  and the expression of hazard risk is the form of:
#'    \eqn{
#'    \lambda(t|X_{i.}) = \exp(\beta^T X_{i.}) \alpha_0(t\exp(\beta^T X_{i.})). \label{eq:riskAFT}
#'    }
#' with \eqn{\alpha_0(t)} is the baseline risk and \eqn{\beta} is the vector of coefficients.
#' The advantage of AFT model is that the variables have a multiplicative effect on \eqn{t} rather than on the risk function, as is the case in Cox model.
#' Two distributions are considered for the baseline risk:
#' \itemize{
#'   \item Weibull: \eqn{\alpha_0(t) = \lambda a t^{(a-1)}};
#'   \item Log-normal: \eqn{\alpha_0(t) = (1/(\sigma\sqrt(2\pi t) \exp[-(\log t - \mu)^2 /2 \sigma^2]))/(1 - \Phi[(\log t - \mu)/\sigma}])}.
#'
#'   To Simulate the covariates, two distributions are also proposed:
#'     \itemize{
#'     \item Uniform
#'     \item Normal
#'     }
#'     and the choice of parameters
#'     The Phi parameter enables to simulate survival data in a linear framework with no interaction, but its future implementation will take into account a non-linear framework with interactions. If the parameter Phi is NULL (to complete...).
#' 3. The hazard risk of the AH model is defined for an individual \eqn{i} as:
#' \eqn{
#'   \lambda_{AH}(t|X_{i.}) = \alpha_0(t\exp(\beta^T X_{i.})),
#'   }
#' with \eqn{\alpha_0} the baseline risk and \eqn{\beta} the vector of regression parameters. In a model with only one binary variable considered that corresponds to the treatment, the hazard risk is written as follows:
#' \eqn{
#'     \lambda_1(t) = \alpha_0(\beta t).
#'     }
#' with \eqn{\alpha_0} the baseline risk and \eqn{\beta} the vector of regression parameters. In a model with only one binary variable considered that corresponds to the treatment, the hazard risk is written as follows:
#' \eqn{
#'     \lambda_1(t) = \alpha_0(\beta t).
#'     }
#' The regression vector \eqn{\beta} characterizes the influence of variables on the survival time of individuals, and \eqn{\exp(\beta^TX_{i.})} is a factor altering the time scale on hazard risk.
#' The positive or negative value of \eqn{\beta^T X_{i.}} will respectively imply an acceleration or deceleration of the risk.The AH model is defined from a linear regression of the interest covariate:
#' Two distributions are considered for the baseline risk:
#' \itemize{
#'   \item Weibull: \eqn{\alpha_0(t) = \lambda a t^{(a-1)}};
#'   \item Log-normal: \eqn{\alpha_0(t) = (1/(\sigma\sqrt(2\pi t) \exp[-(\log t - \mu)^2 /2 \sigma^2]))/(1 - \Phi[(\log t - \mu)/\sigma])}.
#'   }
#' To Simulate the covariates, two distributions are also proposed:
#' \itemize{
#'     \item Uniform
#'     \item Normal
#' }
#' and the choice of parameters
#' The Phi parameter enables to simulate survival data in a linear framework with no interaction, but its future implementation will take into account a non-linear framework with interactions. If the parameter Phi is NULL (to complete...).
#'
#'sim$model <- model
#' @return modelSim returns a list containing: \itemize{
#' \item{model}{ model (Cox, AFT, AFTshift, AH)}
#' \item{Z}{ Matrix of covariates}
#' \item{Y}{ random covariates}
#' \item{TC}{ Vector of survival times}
#' \item{delta}{ Vector of censorship indicator}
#' \item{betanorm}{ Vector of normalized regression parameter}
#' \item{crate}{ Censorship rate}
#' \item{crate_delta}{ Censorship rate}
#' \item{vecY}{ Vector of number of individuals at risk at time \eqn{t_i}}
#' \item{hazParams}{ Vector of parameter distribution of the baseline hazard function}
#' \item{hazDistr}{ Distribution of the baseline hazard function}
#' \item{St}{ Matrix of survival functions}
#' \item{ht}{ Matrix of hazard risk functions}
#' \item{grilleTi}{ Time grid}
#' }
#' @import stats 
#' @importFrom MASS mvrnorm
#' @export
#'
#' @author Mathilde Sautreuil
#' @seealso \code{\link{print.modSim}, \link{plot.modSim}}
#'
#'
#' @examples
#' \dontrun{
#' library(survMS)
#' ### Survival data simulated from Cox model
#' res_paramW = get_param_weib(med = 2228, mu = 2325)
#' listCoxSim_n500_p1000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500,
#'                                 p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "weibull",
#'                                 hazParams = c(res_paramW$a, res_paramW$lambda), seed = 1, d = 0)
#'print(listCoxSim_n500_p1000)
#'hist(listCoxSim_n500_p1000)
#'plot(listCoxSim_n500_p1000, ind = sample(1:500, 5))
#'plot(listCoxSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
#'
#'df_p1000_n500 = data.frame(time = listCoxSim_n500_p1000$TC,
#'                           event = listCoxSim_n500_p1000$delta,
#'                           listCoxSim_n500_p1000$Z)
#'df_p1000_n500[1:6,1:10]
#'dim(df_p1000_n500)
#' ### Survival data simulated from AFT model
#' res_paramLN = get_param_ln(var = 200000, mu = 1134)
#' listAFTSim_n500_p1000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(-1,1), n = 500,
#'                                 p = 100, pnonull = 100, betaDistr = 1, hazDistr = "log-normal",
#'                                 hazParams = c(res_paramLN$a, res_paramLN$lambda),
#'                                 Phi = 0, seed = 1, d = 0)
#' hist(listAFTSim_n500_p1000)
#' plot(listAFTSim_n500_p1000, ind = sample(1:500, 5))
# 'plot(listAFTSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
#' df_p1000_n500 = data.frame(time = listAFTSim_n500_p1000$TC,
#'                            event = listAFTSim_n500_p1000$delta,
#'                            listAFTSim_n500_p1000$Z)
#' df_p1000_n500[1:6,1:10]
#' dim(df_p1000_n500)
#' 
#' ### Survival data simulated from AH model
#' res_paramLN = get_param_ln(var=170000, mu=2325)
#' listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500, 
#'                                  p = 100, pnonull = 100, betaDistr = 1.5, hazDistr = "log-normal",
#'                                  hazParams = c(res_paramLN$a*4, res_paramLN$lambda),
#'                                  Phi = 0, seed = 1, d = 0)
#'                                  
#' print(listAHSim_n500_p1000)
#' hist(listAHSim_n500_p1000)
#' plot(listAHSim_n500_p1000, ind = sample(1:500, 5))
#' plot(listAHSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
#' }
modelSim = function(model = "cox", matDistr, matParam, n, p, pnonull, betaDistr, hazDistr, hazParams, seed, Phi = NULL, d = 0, pourc = 0.9){

  eta_i <- NULL

  testModel = c("cox", "AFT", "AFTshift", "AH")
  if(!any(testModel %in% model)){
    stop("The model must be \"cox\" or \"AFT\" or \"AFTshift\" or \"AH\"")
  }
  ## Survival Models
  TYPES_model<-c("cox", "AFT", "AFTshift", "AH")
  survmodel<-pmatch(model,TYPES_model)

  if(survmodel == 1){
    testhazdistr = c("weibull", "log-normal", "exponential", "gompertz")
    if(!any(testhazdistr %in% hazDistr)){
      stop("The distribution of the matrix must be \"weibull\" or \"log-normal\"")
    }

  }else{
    testhazdistr = c("weibull", "log-normal")
    if(!any(testhazdistr %in% hazDistr)){
      stop("The distribution of the matrix must be \"weibull\" or \"log-normal\"")
    }

  }


  testmatdistr = c("norm", "unif", "mvnorm")
  if(!any(testmatdistr %in% matDistr)){
    stop("The distribution of the matrix must be \"norm\" or \"unif\"")
  }

  if (is.null(matParam) && matDistr == "norm"){
    matParam = c(0,1)
  }else if (is.null(matParam) && matDistr == "unif"){
    matParam = c(-1,1)
  }else if (is.null(matParam) && matDistr == "mvnorm"){
    ## first arg: mu
    ## second arg: rho (coef de cor)
    matParam = c(0,0.6)
  }else if (length(matParam) != 2){
    stop("The length of \"matParam\" must be equal at 2")
  }


  ## Distribution of the matrix
  TYPES_matDistr<-c("norm", "unif", "mvnorm")
  distr<-pmatch(matDistr,TYPES_matDistr)

  ## Distribution of beta
  if (is.character(betaDistr)){
    TYPES_betaDistr<-c("norm", "unif")
    bdistr<-pmatch(betaDistr,TYPES_betaDistr)
    if (bdistr == 1){
      # beta = c(rnorm(pnonull, 0, 1), rep(0, p-pnonull))
      beta = c(rep(1, pnonull), rep(0, p-pnonull))
      beta2 = c(rnorm(pnonull, 0, 3), rep(0, p-pnonull))
    }else if (bdistr == 2){
      beta = c(rep(1, pnonull), rep(0, p-pnonull))
      beta2 = c(runif(pnonull, -1.5, 1.5), rep(0, p-pnonull))
    }
  }else if (is.numeric(betaDistr)){
    beta = c(rep(betaDistr, pnonull), rep(0, p-pnonull))
    beta2 = c(rep(2*betaDistr, pnonull), rep(0, p-pnonull))
  }else{
    stop("error")
  }

  ## Distribution of the baseline hazard
  TYPES_hazDistr<-c("weibull", "log-normal", "gompertz", "exponential")
  hdistr<-pmatch(hazDistr,TYPES_hazDistr)

  if( p == pnonull){
    pp = pnonull
  }else{
    pp = pnonull
  }

  if (is.null(Phi)){
    # hazParams[1] = hazParams[1]/3
    warning("Options \"non-linearity with interactions\" not available")
  }else{
    warning("Options \"non-linearity with interactions\" not available")
  }

  # we simulate a matrix of covariates
  set.seed(seed)
  if(distr == 1){
    Z=matrix(rnorm(n*p,matParam[1],matParam[2]),n,p)
  }else if(distr == 2){
    Z=matrix(runif(n*p,matParam[1],matParam[2]),n,p)
  }else if(distr == 3){
    mup <- c(rep(matParam[1], pnonull), rep(matParam[1]-0.5, p-pnonull))
    rhop <- matParam[2]
    Sigma1 <-  rbind(diag(floor(pnonull/2)) * (1 - rhop) + rhop, matrix(0,p-floor((pnonull/2)),floor(pnonull/2)))
    Sigma2 <-  rbind(matrix(0,floor(pnonull/2),pnonull-floor(pnonull/2)), diag(pnonull-floor(pnonull/2)) * rhop + (1-rhop), matrix(0,p-(pnonull),pnonull-floor(pnonull/2)))
    
    Sigma <- cbind(Sigma1,Sigma2, matrix(0, p, p-pnonull)) 
    diag(Sigma) <- 1
    Z <- mvrnorm(n = n, mup, Sigma)
  }else{
    stop("The matrix distribution is not correct")
  }
  Y = runif(n,0,1)

  betanorm = (1/sqrt(pp))*beta

  if(!is.na(survmodel)){

    if(survmodel == 1){
      ## cox
      if (!is.na(hdistr)){
        if (hdistr == 1){

          # print(hazParams)
          Ts = SurvTimesCoxWeib(Z, beta, Y, pp, hazParams)
          fcts = SurvFctCoxWeib(Z, beta, pp, Ts, hazParams)

          # Ts=(-(1/hazParams[2])*exp((1/sqrt(pp))*(-Z %*% beta))*log(1-Y))^(1/hazParams[1]) #(1/(p*sum(beta)))* (1/(p))* (1/(p*sum(beta)))* (1/(sqrt(5*p)))* (1/sqrt(10000*p))*
          # # print(Ts)
          # tau = max(Ts)
          # pas=100
          # grille_ti=tau*(1/pas)*c(1:(pas))
          # eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
          # h0_t = hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))
          # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
          # H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
          # F_t = 1 - exp(-H0_t*as.vector(eta_i))
          # S_t = exp(-H0_t*as.vector(eta_i))

        }
        else if (hdistr == 2){
          # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
          Ts = SurvTimesCoxLN(Z, beta, Y, pp, hazParams)
          fcts = SurvFctCoxLN(Z, beta, pp, Ts, hazParams)
        }
        else if (hdistr == 3){
          Ts = SurvTimesCoxGomp(Z, beta, Y, pp, hazParams)
          fcts = SurvFctCoxGomp(Z, beta, pp, Ts, hazParams)
        }else if (hdistr == 4){
          Ts = SurvTimesCoxExp(Z, beta, Y, pp, hazParams)
          fcts = SurvFctCoxExp(Z, beta, pp, Ts, hazParams)

        }
        else{
          stop("Distribution not defined")
        }

      }else{
        stop("Distribution not defined")
      }


    }else if(survmodel == 2){
      ## AFT
      if (!is.na(hdistr)){
        if (hdistr == 1){

          Ts = SurvTimesAFTWeib(Z, beta, Y, pp, hazParams)
          fcts = SurvFctAFTWeib(Z, beta, pp, Ts, hazParams)

          # Ts=(-(1/hazParams[2])*exp((1/sqrt(pp))*(-Z %*% beta))*log(1-Y))^(1/hazParams[1]) #(1/(p*sum(beta)))* (1/(p))* (1/(p*sum(beta)))* (1/(sqrt(5*p)))* (1/sqrt(10000*p))*
          # # print(Ts)
          # tau = max(Ts)
          # pas=100
          # grille_ti=tau*(1/pas)*c(1:(pas))
          # eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
          # h0_t = hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))
          # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
          # H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
          # F_t = 1 - exp(-H0_t*as.vector(eta_i))
          # S_t = exp(-H0_t*as.vector(eta_i))

        }
        else if (hdistr == 2){
          # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
          Ts = SurvTimesAFTLN(Z, beta, Y, pp, hazParams)
          fcts = SurvFctAFTLN(Z, beta, pp, Ts, hazParams)
        }
        else{
          stop("Distribution not defined")
        }

      }else{
        stop("Distribution not defined")
      }

    }else if(survmodel == 3){
      ## AFT shift
      if (!is.na(hdistr)){
        if (hdistr == 1){

          Ts = SurvTimesAFTshiftWeib(Z, beta, beta2, Y, pp, hazParams)
          fcts = SurvFctAFTshiftWeib(Z, beta, beta2, pp, Ts, hazParams)

          # Ts=(-(1/hazParams[2])*exp((1/sqrt(pp))*(-Z %*% beta))*log(1-Y))^(1/hazParams[1]) #(1/(p*sum(beta)))* (1/(p))* (1/(p*sum(beta)))* (1/(sqrt(5*p)))* (1/sqrt(10000*p))*
          # # print(Ts)
          # tau = max(Ts)
          # pas=100
          # grille_ti=tau*(1/pas)*c(1:(pas))
          # eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
          # h0_t = hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))
          # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
          # H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
          # F_t = 1 - exp(-H0_t*as.vector(eta_i))
          # S_t = exp(-H0_t*as.vector(eta_i))

        }
        else if (hdistr == 2){
          # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
          Ts = SurvTimesAFTshiftLN(Z, beta, beta2, Y, pp, hazParams)
          fcts = SurvFctAFTshiftLN(Z, beta, beta2, pp, Ts, hazParams)
        }
        else{
          stop("Distribution not defined")
        }

      }else{
        stop("Distribution not defined")
      }

    }else if(survmodel == 4){
      ## AH
      if (!is.na(hdistr)){
        if (hdistr == 1){

          Ts = SurvTimesAHWeib(Z, beta, Y, pp, hazParams)
          fcts = SurvFctAHWeib(Z, beta, pp, Ts, hazParams)

          # Ts=(-(1/hazParams[2])*exp((1/sqrt(pp))*(-Z %*% beta))*log(1-Y))^(1/hazParams[1]) #(1/(p*sum(beta)))* (1/(p))* (1/(p*sum(beta)))* (1/(sqrt(5*p)))* (1/sqrt(10000*p))*
          # # print(Ts)
          # tau = max(Ts)
          # pas=100
          # grille_ti=tau*(1/pas)*c(1:(pas))
          # eta_i = exp((1/sqrt(pp))*(-Z %*% beta))
          # h0_t = hazParams[1]*hazParams[2]*(grille_ti^(hazParams[1]-1))
          # h = matrix(h0_t, nrow = nrow(Z), ncol = pas, byrow = T) * as.vector(eta_i)
          # H0_t = matrix((tau/pas)*cumsum(h0_t), nrow = nrow(Z), ncol = length(grille_ti), byrow = T)
          # F_t = 1 - exp(-H0_t*as.vector(eta_i))
          # S_t = exp(-H0_t*as.vector(eta_i))

        }
        else if (hdistr == 2){
          # Ts<-exp(hazParams[1]*qnorm(1-exp((log(1-Y)/exp((1/(pp))*Z%*%beta))),mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)+hazParams[2])
          Ts = SurvTimesAHLN(Z, beta, Y, pp, hazParams)
          fcts = SurvFctAHLN(Z, beta, pp, Ts, hazParams)
        }
        else{
          stop("Distribution not defined")
        }

      }else{
        stop("Distribution not defined")
      }

    }else{
      stop("Survival models unknown")
    }
  }else{
    stop("Survival models unknown")
  }


  # censoring times
  if (d != 0){
    crate = 1/(d*mean(Ts)); #rate of censoring#
    C = rexp(n, crate);#
  }else{
    crate = 0
    C = Ts
  }

  # observed times
  TC = pmin(Ts, C) # censored observations#
  TC=as.vector(TC);#

  # censoring indicator
  delta = (Ts <= C);#
  delta=as.vector(delta)#
  crate_delta = 1-sum(delta)/n

  # tau=quantile(TC,probs=pourc)
  # tau1=quantile(TC1,probs=0.95)

  don = list(Z=Z, TC=TC, delta=delta)


  matY = matrix(0,n,n);
  for (i in 1:n) {#
    matY[i,] = (don$TC>=don$TC[i])#
  }

  vecY<-apply(matY,2,sum)


  sim <- list()
  sim$model <- model
  sim$Z <- Z
  sim$Y <- Y
  sim$TC <- TC
  sim$delta <- delta
  sim$betaNorm <- betanorm
  sim$crate <- crate
  sim$crate_delta <- crate_delta
  sim$vecY <- vecY
  # sim$tau <- tau
  sim$hazParams <- hazParams
  sim$hazDistr <- hazDistr
  sim$St <- fcts$St
  sim$ht <- fcts$ht
  sim$grilleTi <- fcts$grillet
  class(sim) <- "modSim"

  return(sim)
}

#' Print information about data simulation
#'
#' @param x output of modelSim function (must be of type modSim)
#' @param ... supplementary parameters
#'
#' @return print x
#' @export
#'
#' @examples
#' library(survMS)
#' ### Survival data simulated from AH model
#' res_paramLN = get_param_ln(var=170000, mu=2325)
#' listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500, 
#'                                  p = 100, pnonull = 100, betaDistr = 1.5, hazDistr = "log-normal",
#'                                  hazParams = c(res_paramLN$a*4, res_paramLN$lambda),
#'                                  Phi = 0, seed = 1, d = 0)
#'
#' ### Information about simulation 
#' print(listAHSim_n500_p1000)
print.modSim <- function(x, ...){

  cat("Simulated matrix of size ", dim(x$Z), "\n")
  cat("Distribution of baseline hazard function ", x$hazDistr, "\n")
  cat("Distribution parameter of baseline hazard function ", x$hazParams, "\n")
  cat("Censorship rate", x$crate, "\n")

}

#' Histogram of survival times
#'
#' @param x output of modelSim function (must be of type modSim)
#' @param ... supplementary parameters
#'
#' @return hist x
#' @export
#' @importFrom graphics hist
#'
#' @examples
#' library(survMS)
#' ### Survival data simulated from AH model
#' res_paramLN = get_param_ln(var=170000, mu=2325)
#' listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500, 
#'                                  p = 100, pnonull = 100, betaDistr = 1.5, hazDistr = "log-normal",
#'                                  hazParams = c(res_paramLN$a*4, res_paramLN$lambda),
#'                                  Phi = 0, seed = 1, d = 0)
#'                                  
#' ### Histogram of survival times 
#' hist(listAHSim_n500_p1000)
hist.modSim<- function(x, ...){
  hist(x$TC, xlab = "times", main = "Histogram of survival times")
}

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
#' res_paramW = get_param_weib(med = 2.5, mu = 1.2)
#' listCoxSimCor_n500_p1000 <- modelSim(model = "cox", matDistr = "mvnorm", 
#'                                      matParam = c(0,0.6), n = 500, 
#'                                      p = 1000, pnonull = 20, betaDistr = 1, 
#'                                      hazDistr = "weibull", 
#'                                      hazParams = c(res_paramW$a, res_paramW$lambda), 
#'                                      seed = 1, d = 0)
#' print(listCoxSimCor_n500_p1000)
#' hist(listCoxSimCor_n500_p1000)
#' # Heatmap(listCoxSimCor_n500_p1000, k = 4)
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

#' Draw Heatmap of modSim object
#'
#' @param x modSim object
#' @param k number of split for heatmap's columns
#' @param ind number of columns to keep
#' @param ... supplementary parameters
#'
#' @return draw x
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
#' draw(listCoxSimCor_n500_p1000, k = 3)
draw.modSim<- function(x, k, ind = NULL, ...){
  
  col_fun = colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
  if(!is.null(ind)){
    ind = ind
  }else{
    ind = which(x$betaNorm != 0)
  }
  rows_info <- rep("High", nrow(x$Z))
  rows_info[which(x$TC > median(x$TC))] <- "Low"
  colnames(x$Z) <- paste0("X", 1:ncol(x$Z))
  ht <- Heatmap(as.matrix(x$Z)[,ind], name = "expression", 
          row_split = rows_info, 
          # column_split = c(rep("Sign", sum(x$betaNorm[ind] != 0)),
          #                  rep("No Sign", 
          #                      ncol(x$Z[,ind]) - sum(x$betaNorm[ind] != 0))),
          col = col_fun, #row_km = 2,
          show_column_names = TRUE, column_km = k)
  draw(ht, heatmap_legend_side = "bottom")
  
}

#' Survival or hazard curves of simulated data
#'
#' @param x output of modelSim function (must be of type modSim)
#' @param ind vector (individuals to show)
#' @param ... supplementary parameters
#' @param type type of plots (survival or hazard curves)
#'
#' @return plot x
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_classic
#'
#' @examples
#' library(survMS)
#' ind = sample(1:500, 5)
#' ### Example with survival data simulated from AH model
#' res_paramLN = get_param_ln(var=170000, mu=2325)
#' listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500, 
#'                                  p = 100, pnonull = 100, betaDistr = 1.5, hazDistr = "log-normal",
#'                                  hazParams = c(res_paramLN$a*4, res_paramLN$lambda),
#'                                  Phi = 0, seed = 1, d = 0)
#' ### Two types of plot are available (survival (by default) and hazard curves)
#' ## Survival curves                                
#' plot(listAHSim_n500_p1000, ind = ind)
#' ## Hazard curves
#' plot(listAHSim_n500_p1000, ind = ind, type = "hazard")
plot.modSim <- function(x, ind, type = "surv", ...){

  testTYPE = c("surv", "hazard")
  if(!any(testTYPE %in% type)){
    stop("The type of plot must be \"surv\" or \"hazard\"")
  }
  ## Plot type
  TYPES_plot<-c("surv", "hazard")
  typeplot<-pmatch(type,TYPES_plot)


  if(typeplot == 1){
    temps <- Sr <- NULL
    ind_random = ind
    df_Ft10 = data.frame(Sr = as.vector(x$St), ind = rep(1:nrow(x$Z), length(x$grilleTi)), temps = rep(x$grilleTi, each = nrow(x$Z)), ti = rep(x$TC, length(x$grilleTi)))#hr = as.vector(h),
    p <- ggplot(aes(x = temps, y = Sr, color = as.factor(ind)), data = df_Ft10[which(df_Ft10$ind %in% ind_random),]) + geom_line() + #+ geom_point(aes(y = 0.99, x = ti, color = as.factor(ind)))
      labs(x="Times", #title="Survival curves for 5 individuals",
           y = "Survival probability", color = "Individuals")+ theme_classic()
    p

  }else if(typeplot == 2){
    temps <- hr <- NULL
    ind_random = ind
    df_Ft10 = data.frame(hr = as.vector(x$ht), ind = rep(1:nrow(x$Z), length(x$grilleTi)), temps = rep(x$grilleTi, each = nrow(x$Z)), ti = rep(x$TC, length(x$grilleTi)))#hr = as.vector(h),
    p <- ggplot(aes(x = temps, y = hr, color = as.factor(ind)), data = df_Ft10[which(df_Ft10$ind %in% ind_random),]) + geom_line() + #+ geom_point(aes(y = 0.99, x = ti, color = as.factor(ind)))
      labs(x="Times", y = "Hazard risk", color = "Individuals")+ theme_classic() #title="Hazard risk curves for 5 individuals",
    p
  }else{
    stop("wrong plot type")
  }

}



