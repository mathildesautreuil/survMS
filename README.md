survMS R package: survival Model Simulation
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- <img src="docs/reference/survMS.png" align="right" height=150/> -->

<img src="https://raw.githubusercontent.com/mathildesautreuil/survMS/master/docs/reference/survMS.png" height="2" align="right" />

<!-- [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/survMS)](https://CRAN.R-project.org/package=survMS) -->

[![](https://www.r-pkg.org/badges/version/survMS?color=bc8f8f)](https://cran.r-project.org/package=survMS)
[![](http://cranlogs.r-pkg.org/badges/grand-total/survMS?color=8fbcbc)](https://cran.r-project.org/package=survMS)
<!-- [![](http://cranlogs.r-pkg.org/badges/grand-total/survMS?color=blue)](https://cran.r-project.org/package=survMS) -->

<!-- [![R build status](https://github.com/tidyverse/dplyr/workflows/R-CMD-check/badge.svg)](https://github.com/mathildesautreuil/survMS/actions?workflow=R-CMD-check) -->
<!-- [![R build status](https://github.com/mathildesautreuil/survMS/workflows/R-CMD-check/badge.svg)](https://github.com/mathildesautreuil/survMS/actions) -->

[![R build
status](https://github.com/mathildesautreuil/survMS//workflows/R-CMD-check/badge.svg)](https://github.com/mathildesautreuil/survMS//actions)
<!-- [![](https://img.shields.io/github/last-commit/mathildesautreuil/survMS.svg)](https://github.com/mathildesautreuil/survMS/commits/master) -->

[![](https://img.shields.io/badge/survMS-website-bc8f8f.svg)](https://mathildesautreuil.github.io/survMS)

# Installation

Install from CRAN:

``` r
install.packages("survMS")
```

Or install the development version from Github:

``` r
install.packages("devtools")
devtools::install_github("mathildesautreuil/survMS")
```

And load the library:

``` r
library(survMS)
#> Loading required package: ggplot2
#> Loading required package: ComplexHeatmap
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.13.1
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
#> Loading required package: circlize
#> ========================================
#> circlize version 0.4.15
#> CRAN page: https://cran.r-project.org/package=circlize
#> Github page: https://github.com/jokergoo/circlize
#> Documentation: https://jokergoo.github.io/circlize_book/book/
#> 
#> If you use it in published research, please cite:
#> Gu, Z. circlize implements and enhances circular visualization
#>   in R. Bioinformatics 2014.
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(circlize))
#> ========================================
#> 
#> Attaching package: 'survMS'
#> The following object is masked from 'package:ComplexHeatmap':
#> 
#>     Heatmap
```

# 1 Introduction

**survMS** is an R package that generates survival data from different
survival models with varying levels of complexity. Three survival models
are implemented in the package to simulate survival data: a Cox model,
an Accelerated Failure Time (AFT) model, and an Accelerated Hazard (AH)
model. These models have their specificities. Indeed, a Cox model is a
proportional risk model, not an AFT and an AH model. A proportional risk
model means that the ratio of the risks between two individuals with
different explanatory variables does not depend on the time. In an AFT
model, the explanatory variables have an accelerating or decelerating
effect on individuals’ survival. In both models, a Cox model and an AFT
model, the survival function’s curves never intersect. When one is
interested in more complex data with intersecting survival curves, for
example, to make it more difficult for methods to predict survival, two
approaches are implemented in the package. The first approach consists
of modifying an AFT model to have intersecting survival curves. The
second approach concerns using of an AH model (for *Accelerated
Hazards*) to generate the survival data. an AH model is more flexible
than the two models mentioned above. In an AH model, the variables
accelerate or decelerate the hazard risk. The survival curves of an AH
model can therefore cross each other. The generation of survival times
carried out from the different models mentioned above. The baseline risk
function of the models is assumed to be known and follows a particular
probability distribution. For these different models, we use parametric
distribution of the baseline risk for the generation of survival time.

Note that we stayed in a linear dependency framework without
interaction. We will adapt the package in a second step to be in a
nonlinear framework to simulate interactions. The package’s
implementation is easily adapted to the nonlinear framework with
interactions by replacing the explanatory variables’ linear part with a
nonlinear function of the explanatory variables. All the functions in
the package are coded to introduce this nonlinear function of the
explanatory variables easily.

For details on methodology, see
[Vignette](https://mathildesautreuil.github.io/survMS/articles/how-to-simulate-survival-models.html).

# 2 Simulation from a Cox model

The package survMS can simulate data from a Cox model. We have chosen to
carry out this simulation to generate survival data that respects the
proportional risk hypothesis. For more details, see
[Vignette](https://mathildesautreuil.github.io/survMS/articles/how-to-simulate-survival-models.html).
Survival times can follow different distributions: Weibull, log-normal,
exponential (in progress) and Gompertz (in progress).

## 2.1 Simulation from a Cox model with baseline hazard following Weibull distribution

To simulate survival data checking the proportional risk hypothesis, we
use a Cox model with survival times following a Weibull distribution.

<!-- We recall that the generation of survival data from a Cox model based on :  -->
<!-- $$ -->
<!--     T = H_0^{-1} \left[ \frac{-\log(1-U)}{\exp(\beta^T X_{i.})}\right], -->
<!-- $$ -->
<!-- where $U \sim \mathcal{U} [0.1].$ -->
<!-- For this simulation, we consider that survival times follow a Weibull's law $\mathcal{W}(a,\lambda).$ In this case, we have the cumulative risk function expressed by :  -->
<!-- $$ -->
<!--     H_0(t) = - \log\left [1 - \Phi\left(\frac{\log t - \mu}{\sigma}\right)\right] -->
<!-- $$ -->
<!-- and survival times can therefore be simulated from :  -->
<!-- $$ -->
<!--     T = \frac{1}{\lambda^{1/a}} \left( \frac{-\log(1-U)}{\exp(\beta^T X_{i.})} \right)^{1/a} .  -->
<!-- $$ -->

``` r
res_paramW = get_param_weib(med = 1062, mu = 1134)
# c(res_paramW$a, res_paramW$lambda)
listCoxWSim_n500_p1000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500, 
                                   p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "weibull", 
                                   hazParams = c(3, 4), seed = 1, d = 0)
#> Warning in modelSim(model = "cox", matDistr = "unif", matParam = c(-1, 1), :
#> Options "non-linearity with interactions" not available
hist(listCoxWSim_n500_p1000)
```

<img src="man/figures/README-sim_cox_model-1.png" width="100%" />

## 2.2 Simulation from a Cox model with baseline hazard following log-normal distribution

To simulate survival data checking the proportional risk hypothesis, we
use a Cox model with survival times following a log-normal distribution.

``` r
res_paramLN = get_param_ln(var = 170000, mu = 2320)
listCoxSim_n500_p1000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500, 
                                  p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "log-normal", 
                                  hazParams = c(res_paramLN$a, res_paramLN$lambda), seed = 1, d = 0)
#> Warning in modelSim(model = "cox", matDistr = "unif", matParam = c(-1, 1), :
#> Options "non-linearity with interactions" not available
hist(listCoxSim_n500_p1000)
```

<img src="man/figures/README-sim_coxLN_model-1.png" width="100%" />

## 2.3 Simulation from a Cox model with a correlated covariate matrix

For this simulation, we consider the same design as 2.1 to simulate
survival times. However, the simulated covariate matrix is correlated.
The pertinent covariates are those correlated. We split this set of
pertinent covariates in two subsets: one with a correlation parameter
(*matParam\[2\]*) and the other one with (1 - *matParam\[2\]*). The
non-pertinent covariates are uncorrelated.

``` r
res_paramW = get_param_weib(med = 2.5, mu = 1.2)
listCoxSimCor_n500_p1000 <- modelSim(model = "cox", matDistr = "mvnorm", matParam = c(0,0.6), n = 500,
                                  p = 1000, pnonull = 10, betaDistr = 1, hazDistr = "weibull",
                                  hazParams = c(res_paramW$a, res_paramW$lambda), seed = 1, d = 0)
#> Warning in modelSim(model = "cox", matDistr = "mvnorm", matParam = c(0, :
#> Options "non-linearity with interactions" not available
Heatmap(listCoxSimCor_n500_p1000, ind = 1:15, k = 3)
```

<img src="man/figures/README-sim_cox_model_cor-1.png" width="100%" />

# 3 Simulation from an AFT model

From survMS package, we can simulate the data from an AFT model. We have
chosen to carry out this simulation to generate survival data that does
not respect the proportional risk hypothesis. For more details, see
[Vignette](https://mathildesautreuil.github.io/survMS/articles/how-to-simulate-survival-models.html).
Survival times can follow different distributions: Weibull, log-normal,
exponential (in progress) and Gompertz (in progress).

## 3.1 Simulation from an AFT model with baseline hazard following Weibull distribution

To simulate survival data not checking the proportional risk hypothesis,
we use an AFT model with survival times following a weibull
distribution.

``` r
res_paramW = get_param_weib(med = 1062, mu = 1134)
# c(res_paramW$a, res_paramW$lambda)
listAFTWSim_n500_p1000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(-1,1), n = 500, 
                                  p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "weibull", 
                                  hazParams = c(3, 2), seed = 1, d = 0)
#> Warning in modelSim(model = "AFT", matDistr = "unif", matParam = c(-1, 1), :
#> Options "non-linearity with interactions" not available
hist(listAFTWSim_n500_p1000)
```

<img src="man/figures/README-sim_aftW_model-1.png" width="100%" />

## 3.2 Simulation from an AFT model with baseline hazard following log-normal distribution

To simulate survival data not checking the proportional risk hypothesis,
we use an AFT model with survival times following a log-normal
distribution.

``` r
res_paramLN = get_param_ln(var = 170000, mu = 2325)
listAFTLNSim_n500_p1000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(-1,1), n = 500, 
                                  p = 10, pnonull = 10, betaDistr = 1, hazDistr = "log-normal", 
                                  hazParams = c(res_paramLN$a*3, res_paramLN$lambda), seed = 1, d = 0)
#> Warning in modelSim(model = "AFT", matDistr = "unif", matParam = c(-1, 1), :
#> Options "non-linearity with interactions" not available
hist(listAFTLNSim_n500_p1000)
```

<img src="man/figures/README-sim_aftLN_model-1.png" width="100%" />

<!-- ```{r sim_aft_model} -->
<!-- res_paramLN = get_param_ln(var = 200000, mu = 1134)#compute_param_weibull(a_list = a_list, med = 2280, moy = 2325, var = 1619996) -->
<!-- listAFTLNSim_n500_p1000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(0,1), n = 500,  -->
<!--                                   p = 100, pnonull = 100, betaDistr = 1, hazDistr = "log-normal", -->
<!--                                   hazParams = c(0.3, res_paramLN$lambda), Phi = 0, seed = 1, d = 0) -->
<!-- hist(listAFTLNSim_n500_p1000) -->
<!-- ``` -->

## 3.3 Simulation from a shifted AFT model

To obtain survival data from a model that does not respect the
proportional risk hypothesis and that their survival curves can not
cross, we have modified the AFT model by adding a term *betaDistr* in
the model. For more details, see
[Vignette](https://mathildesautreuil.github.io/survMS/articles/how-to-simulate-survival-models.html).
Survival times can follow different distributions: Weibull, log-normal,
exponential (in progress) and Gompertz (in progress).

``` r
res_paramLN = get_param_ln(var = 170000, mu = 2325)
listAFTsSim_n500_p1000 <- modelSim(model = "AFTshift", matDistr = "unif", matParam = c(-1,1), n = 500, 
                                  p = 10, pnonull = 10, betaDistr = "unif", hazDistr = "log-normal", 
                                  hazParams = c(0.3, res_paramLN$lambda), seed = 1, d = 0)
#> Warning in modelSim(model = "AFTshift", matDistr = "unif", matParam = c(-1, :
#> Options "non-linearity with interactions" not available
#> Warning in log(mat_grille_ti * (1/as.vector(eta_i)) - matrix(phi2, nrow(Z), :
#> production de NaN
hist(listAFTsSim_n500_p1000)
```

<img src="man/figures/README-sim_aft_mod_model-1.png" width="100%" />

# 4 Simulation from an AH model

From survMS package, we can also simulated survival data from another
model, an AH model. We have chosen to carry out this simulation to
generate survival data that does not respect the proportional risk
hypothesis and survival curves can cross. For more details, see
[Vignette](https://mathildesautreuil.github.io/survMS/articles/how-to-simulate-survival-models.html).

## 4.1 Simulation from an AH model with baseline hazard following Weibull distribution

``` r
res_paramW = get_param_weib(med = 1062, mu = 1134)
listAHwSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500,
                                  p = 100, pnonull = 20, betaDistr = 1, hazDistr = "weibull",
                                  hazParams = c(res_paramW$a, res_paramW$lambda),
                                  Phi = 0, seed = 1, d = 0)
#> Warning in modelSim(model = "AH", matDistr = "unif", matParam = c(-1, 1), :
#> Options "non-linearity with interactions" not available
hist(listAHwSim_n500_p1000)
```

<img src="man/figures/README-sim_ah_model-1.png" width="100%" />

## 4.2 Simulation from the an AH model with baseline hazard following log-normal distribution

``` r
res_paramLN = get_param_ln(var = 170000, mu = 2325)
listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500,
                                  p = 100, pnonull = 100, betaDistr = 1.5, hazDistr = "log-normal",
                                  hazParams = c(res_paramLN$a*4, res_paramLN$lambda), Phi = 0, 
                                 seed = 1, d = 0)
#> Warning in modelSim(model = "AH", matDistr = "unif", matParam = c(-1, 1), :
#> Options "non-linearity with interactions" not available
# summary(listAHSim_n500_p1000)
# print(listAHSim_n500_p1000)
hist(listAHSim_n500_p1000)
```

<img src="man/figures/README-sim_ahLN_model-1.png" width="100%" />

# 5 Utils functions

<!-- The package enables simulating data close to real datasets. In this goal, we have to know the mean, median, and variance of real datasets and then compute the parameters of the used distribution to have a similar distribution of survival times by giving the mean, median, and variance. -->

## 5.1 Get parameters of survival time distribution

``` r
# Weibull distribution
res_paramW = get_param_weib(med = 1062, mu = 1134)
res_paramW
#> $a
#> [1] 1.969765
#> 
#> $lambda
#> [1] 7.586963e-07
# Log-normale distribution
res_paramLN = get_param_ln(var = 600000, mu = 1134)
res_paramLN
#> $a
#> [1] 0.6188154
#> 
#> $lambda
#> [1] 6.84204
```

## 5.2 Plotting survival curves

### 5.2.1 Cox/Weibull model

``` r
set.seed(1234)
ind = sample(x = 1:500, 8)
## Cox/Weibull model
# df_p1000_n500[1:6,1:10]
plot(listCoxWSim_n500_p1000, ind = ind)
```

<img src="man/figures/README-surv_curves-1.png" width="100%" />

### 5.2.2 AFT/Log-normal model

``` r
## AFT/Weibull model
# df_p1000_n500[1:6,1:10]
plot(listAFTLNSim_n500_p1000, ind = ind)
```

<img src="man/figures/README-surv_curves_aft-1.png" width="100%" />

### 5.2.3 Shifted AFT/Log-normal model

``` r
## Shifted AFT/LN model
# df_p1000_n500[1:6,1:10]
plot(listAFTsSim_n500_p1000, ind = ind)
```

<img src="man/figures/README-surv_curves_afts-1.png" width="100%" />

### 5.2.4 AH/LN model

``` r
## Cox/Weibull model
# df_p1000_n500[1:6,1:10]
# set.seed(765)
ind = sample(x = 1:500, 8)
plot(listAHSim_n500_p1000, ind = ind)
```

<img src="man/figures/README-surv_curves_ah-1.png" width="100%" />

## 5.3 Plotting hazard curves

### 5.3.1 Cox/Weibull model

``` r
## Cox/Weibull model
# df_p1000_n500[1:6,1:10]
plot(listCoxWSim_n500_p1000, ind = ind, type = "hazard")
```

<img src="man/figures/README-haz_curves-1.png" width="100%" />

### 5.3.2 AFT/Weibull model

``` r
## AFT/Weibull model
# df_p1000_n500[1:6,1:10]
plot(listAFTWSim_n500_p1000, ind = ind, type = "hazard")
```

<img src="man/figures/README-haz_curves_aft-1.png" width="100%" />

### 5.3.3 AFT/Log-normal model

``` r
## AFT/LN model
# df_p1000_n500[1:6,1:10]
plot(listAFTLNSim_n500_p1000, ind = ind, type = "hazard")
```

<img src="man/figures/README-haz_curves_aftLN-1.png" width="100%" />

### 5.3.4 Shifted AFT/Log-normal model

``` r
## Shifted AFT/LN model
# df_p1000_n500[1:6,1:10]
plot(listAFTsSim_n500_p1000, ind = ind, type = "hazard")
```

<img src="man/figures/README-haz_curves_afts-1.png" width="100%" />

### 5.3.5 AH/LN model

``` r
## Cox/Weibull model
# df_p1000_n500[1:6,1:10]
plot(listAHSim_n500_p1000, ind = ind, type = "hazard")
```

<img src="man/figures/README-haz_curves_ah-1.png" width="100%" />

# 6 Examples

Please also see the package
[vignette](https://mathildesautreuil.github.io/survMS/articles/how-to-simulate-survival-models.html)
for more detailed examples and a description of the methodology
underpinning the **survMS** package.
