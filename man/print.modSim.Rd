% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelSim.R
\name{print.modSim}
\alias{print.modSim}
\title{Print information about data simulation}
\usage{
\method{print}{modSim}(x, ...)
}
\arguments{
\item{x}{output of modelSim function (must be of type modSim)}

\item{...}{supplementary parameters}
}
\value{
print x
}
\description{
Print information about data simulation
}
\examples{
library(survMS)
### Survival data simulated from AH model
res_paramLN = get_param_ln(var=170000, mu=2325)
listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500, 
                                 p = 100, pnonull = 100, betaDistr = 1.5, hazDistr = "log-normal",
                                 hazParams = c(res_paramLN$a*4, res_paramLN$lambda),
                                 Phi = 0, seed = 1, d = 0)

### Information about simulation 
print(listAHSim_n500_p1000)
}
