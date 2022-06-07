library(survMS)

### Survival data simulated from Cox model
res_paramW = get_param_weib(med = 2228, mu = 2325)
listCoxSim_n500_p1000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500,
                                p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "weibull",
                                hazParams = c(res_paramW$a, res_paramW$lambda), seed = 1, d = 0)
print(listCoxSim_n500_p1000)
hist(listCoxSim_n500_p1000)
plot(listCoxSim_n500_p1000, ind = sample(1:500, 5))
plot(listCoxSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")

df_p1000_n500 = data.frame(time = listCoxSim_n500_p1000$TC,
                          event = listCoxSim_n500_p1000$delta,
                          listCoxSim_n500_p1000$Z)
df_p1000_n500[1:6,1:10]
dim(df_p1000_n500)
### Survival data simulated from AFT model
res_paramLN = get_param_ln(var = 200000, mu = 1134)
listAFTSim_n500_p1000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(-1,1), n = 500,
                                p = 100, pnonull = 100, betaDistr = 1, hazDistr = "log-normal",
                                hazParams = c(res_paramLN$a, res_paramLN$lambda),
                                Phi = 0, seed = 1, d = 0)
hist(listAFTSim_n500_p1000)
plot(listAFTSim_n500_p1000, ind = sample(1:500, 5))
plot(listAFTSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
df_p1000_n500 = data.frame(time = listAFTSim_n500_p1000$TC,
                           event = listAFTSim_n500_p1000$delta,
                           listAFTSim_n500_p1000$Z)
df_p1000_n500[1:6,1:10]
dim(df_p1000_n500)

### Survival data simulated from AH model
res_paramLN = get_param_ln(var=170000, mu=2325)
listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500,
                                 p = 100, pnonull = 100, betaDistr = 1.5, hazDistr = "log-normal",
                                 hazParams = c(res_paramLN$a*4, res_paramLN$lambda),
                                 Phi = 0, seed = 1, d = 0)

print(listAHSim_n500_p1000)
hist(listAHSim_n500_p1000)
plot(listAHSim_n500_p1000, ind = sample(1:500, 5))
plot(listAHSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")

### Correlated survival data simulated from Cox model
res_paramW = get_param_weib(med = 2.5, mu = 1.2)
listCoxSimCor_n500_p1000 <- modelSim(model = "cox", matDistr = "mvnorm", matParam = c(0,0.6), n = 500,
                                  p = 1000, pnonull = 5, betaDistr = 1, hazDistr = "weibull",
                                  hazParams = c(res_paramW$a, res_paramW$lambda), seed = 1, d = 0)
print(listCoxSimCor_n500_p1000)
hist(listCoxSimCor_n500_p1000)
Heatmap(listCoxSimCor_n500_p1000)

# ### Survival data simulated from Cox model
# res_paramW = get_param_weib(med = 1062, mu = 1134)
# listCoxSim_n500_p1000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                 p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "weibull",
#                                 hazParams = c(1.5, 1), seed = 1, d = 0)
# print(listCoxSim_n500_p1000)
# hist(listCoxSim_n500_p1000)
# plot(listCoxSim_n500_p1000, ind = sample(1:500, 5))
# plot(listCoxSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
# df_p1000_n500 = data.frame(time = listCoxSim_n500_p1000$TC,
#                           event = listCoxSim_n500_p1000$delta,
#                           listCoxSim_n500_p1000$Z)
# df_p1000_n500[1:6,1:10]
# dim(df_p1000_n500)
# 
# ## COX/LN model
# listCoxSim_n500_p1000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                   p = 100, pnonull = 20, betaDistr = 1, hazDistr = "log-normal",
#                                   hazParams = c(0.25,1), seed = 1)#c(res_paramLN$a*3, res_paramLN$lambda), seed = 1, d = 0)
# print(listCoxSim_n500_p1000)
# hist(listCoxSim_n500_p1000)
# plot(listCoxSim_n500_p1000, ind = sample(1:500, 5))
# plot(listCoxSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
# 
# ## COX/exponetial model
# res_paramE = get_param_exp(med = 1062, mu = 1134)
# listCoxSim_n500_p1000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                   p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "exponential",
#                                   hazParams = res_paramE$lambda, seed = 1, d = 0)
# print(listCoxSim_n500_p1000)
# hist(listCoxSim_n500_p1000)
# plot(listCoxSim_n500_p1000, ind = sample(1:500, 5))
# 
# df_p1000_n500 = data.frame(time = listCoxSim_n500_p1000$TC,
#                            event = listCoxSim_n500_p1000$delta,
#                            listCoxSim_n500_p1000$Z)
# df_p1000_n500[1:6,1:10]
# dim(df_p1000_n500)
# 
# listCoxSim_n500_p25000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                 p = 25000, pnonull = 20, betaDistr = 1, hazDistr = "weibull",
#                                 hazParams = c(3, 4), seed = 1, d = 0)
# print(listCoxSim_n500_p25000)
# hist(listCoxSim_n500_p25000)
# plot(listCoxSim_n500_p25000, ind = sample(1:500, 5))
# plot(listCoxSim_n500_p25000, ind = sample(1:500, 5), type = "hazard")
# 
# df_p25m_n500 = data.frame(time = listCoxSim_n500_p25000$TC,
#                          event = listCoxSim_n500_p25000$delta,
#                          listCoxSim_n500_p25000$Z)
# df_p25m_n500[1:6,1:10]
# 
# ### Survival data simulated from AFT model
# res_paramLN = get_param_ln(var = 200000, mu = 1134)
# listAFTSim_n500_p1000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                 p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "log-normal",
#                                 hazParams = c(0.25, 0.5), Phi = 0, seed = 1, d = 0)
# print(listAFTSim_n500_p1000)
# hist(listAFTSim_n500_p1000)
# plot(listAFTSim_n500_p1000, ind = sample(1:500, 10), type = "surv")
# plot(listAFTSim_n500_p1000, ind = sample(1:500, 10), type = "hazard")
# df_p1000_n500 = data.frame(time = listAFTSim_n500_p1000$TC,
#                            event = listAFTSim_n500_p1000$delta,
#                            listAFTSim_n500_p1000$Z)
# df_p1000_n500[1:6,1:10]
# dim(df_p1000_n500)
# res_paramW = get_param_weib(med = 1062, mu = 1134)
# listAFTSim_n500_p1000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(0,1), n = 500,
#                                   p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "weibull",
#                                   hazParams = c(1.5, 1),
#                                   Phi = 0, seed = 1, d = 0)
# print(listAFTSim_n500_p1000)
# hist(listAFTSim_n500_p1000)
# plot(listAFTSim_n500_p1000, ind = sample(1:500, 5))
# plot(listAFTSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
# df_p1000_n500 = data.frame(time = listAFTSim_n500_p1000$TC,
#                            event = listAFTSim_n500_p1000$delta,
#                            listAFTSim_n500_p1000$Z)
# df_p1000_n500[1:6,1:10]
# dim(df_p1000_n500)
# 
# ### Survival data simulated from Shifted  AFT model
# res_paramLN = get_param_ln(var=170000, mu=2325)
# listAFTSim_n500_p1000 <- modelSim(model = "AFTshift", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                   p = 100, pnonull = 100, betaDistr = "unif", hazDistr = "log-normal",
#                                   hazParams = c(0.3, res_paramLN$lambda),
#                                   Phi = 0, seed = 100, d = 0)
# print(listAFTSim_n500_p1000)
# hist(listAFTSim_n500_p1000)
# plot(listAFTSim_n500_p1000, ind = sample(1:500, 10), type = "surv")
# plot(listAFTSim_n500_p1000, ind = sample(1:500, 10), type = "hazard")
# df_p1000_n500 = data.frame(time = listAFTSim_n500_p1000$TC,
#                            event = listAFTSim_n500_p1000$delta,
#                            listAFTSim_n500_p1000$Z)
# df_p1000_n500[1:6,1:10]
# dim(df_p1000_n500)
# 
# res_paramW = get_param_weib(med = 1062, mu = 1134)
# listAFTSim_n500_p1000 <- modelSim(model = "AFTshift", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                   p = 1000, pnonull = 20, betaDistr = "unif", hazDistr = "weibull",
#                                   hazParams = c(res_paramW$a, res_paramW$lambda),
#                                   Phi = 0, seed = 1, d = 0)
# print(listAFTSim_n500_p1000)
# hist(listAFTSim_n500_p1000)
# plot(listAFTSim_n500_p1000, ind = sample(1:500, 5))
# plot(listAFTSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
# df_p1000_n500 = data.frame(time = listAFTSim_n500_p1000$TC,
#                            event = listAFTSim_n500_p1000$delta,
#                            listAFTSim_n500_p1000$Z)
# df_p1000_n500[1:6,1:10]
# dim(df_p1000_n500)
# 
# # write.csv(df_p1000_n500, file = paste0(path, "DF_RegScreen_n500_p1000.csv"), row.names = F)
# # listAFTSim_n500_p25000 <- AFTSim(matDistr = "unif", matParam = c(-1,1), n = 500, p = 25000,
#                                  # pnonull = 20, betaDistr = 1, hazDistr = "weibull",
#                                  # hazParams = c(res_paramLN$a, res_paramLN$lambda),
#                                  # seed = 1, d = 0)
# # hist(listAFTSim_n500_p25000$TC)
# # df_p25m_n500 = data.frame(time = listAFTSim_n500_p25000$TC,
#                       # event = listAFTSim_n500_p25000$delta,
#                       # listAFTSim_n500_p25000$Z)
# # df_p25m_n500[1:6,1:10]
# 
# ### Survival data simulated from Shifted AH model
# res_paramW = get_param_weib(med = 1062, mu = 1134)
# listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                   p = 1000, pnonull = 20, betaDistr = 1, hazDistr = "weibull",
#                                   hazParams = c(res_paramW$a, res_paramW$lambda),
#                                   Phi = 0, seed = 1, d = 0)
# print(listAHSim_n500_p1000)
# hist(listAHSim_n500_p1000)
# plot(listAHSim_n500_p1000, ind = sample(1:500, 5))
# plot(listAHSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
# df_p1000_n500 = data.frame(time = listAHSim_n500_p1000$TC,
#                            event = listAHSim_n500_p1000$delta,
#                            listAHSim_n500_p1000$Z)
# df_p1000_n500[1:6,1:10]
# dim(df_p1000_n500)
# 
# res_paramLN = get_param_ln2(med = 2280, mu = 2325)
# res_paramLN = get_param_ln(var=170000, mu=2325)
# listAHSim_n500_p1000 <- modelSim(model = "AH", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                  p = 10, pnonull = 10, betaDistr = 1, hazDistr = "log-normal",
#                                  hazParams = c(res_paramLN$a*3, res_paramLN$lambda),
#                                  Phi = 0, seed = 1, d = 0)
# print(listAHSim_n500_p1000)
# hist(listAHSim_n500_p1000)
# plot(listAHSim_n500_p1000, ind = sample(1:500, 5))
# plot(listAHSim_n500_p1000, ind = sample(1:500, 5), type = "hazard")
# ## test
# 
# # write.csv(df_p1000_n500, file = paste0(path, "DF_RegScreen_n500_p1000.csv"), row.names = F)
# system.time(listAFTSim_n500_p25000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(-1,1),
#                                                n = 500, p = 25000, pnonull = 20, betaDistr = 1,
#                                                hazDistr = "weibull",
#                                                hazParams = c(res_paramLN$a, res_paramLN$lambda),
#                                                seed = 1, d = 0))
# hist(listAFTSim_n500_p25000)
# df_p25m_n500 = data.frame(time = listAFTSim_n500_p25000$TC,
# event = listAFTSim_n500_p25000$delta,
# listAFTSim_n500_p25000$Z)
# df_p25m_n500[1:6,1:10]
# 
# 
# system.time(listCoxSim_n500_p25000 <- modelSim(model = "cox", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                 p = 25000, pnonull = 20, betaDistr = 1, hazDistr = "weibull",
#                                 hazParams = c(res_paramW$a, res_paramW$lambda), seed = 1, d = 0))
# print(listCoxSim_n500_p25000)
# hist(listCoxSim_n500_p25000)
# plot(listCoxSim_n500_p25000, ind = sample(1:500, 5))
# 
# df_p25m_n500 = data.frame(time = listCoxSim_n500_p25000$TC,
#                          event = listCoxSim_n500_p25000$delta,
#                          listCoxSim_n500_p25000$Z)
# df_p25m_n500[1:6,1:10]
# 
# ### Survival data simulated from AFT model
# res_paramLN = get_param_ln(var = 200000, mu = 1134)
# system.time(listAFTSim_n500_p1000 <- modelSim(model = "AFT", matDistr = "unif", matParam = c(-1,1), n = 500,
#                                 p = 100, pnonull = 100, betaDistr = 1, hazDistr = "log-normal",
#                                 hazParams = c(res_paramLN$a, res_paramLN$lambda),
#                                 Phi = 0, seed = 1, d = 0))
# hist(listAFTSim_n500_p1000)
# df_p1000_n500 = data.frame(time = listAFTSim_n500_p1000$TC,
#                            event = listAFTSim_n500_p1000$delta,
#                            listAFTSim_n500_p1000$Z)
# df_p1000_n500[1:6,1:10]
# dim(df_p1000_n500)