
setwd("C:\\R Working Directory\\HPIF_Scenarios")

library(BVAR)
library(tibble)
library(dplyr)
library(tidyr)
library(zoo)
library(xts)
library(ggplot2)
library(plot.matrix)

source("Nelson_Siegel_Functions.R")


load('mth.series.tbl.RData')
#load('NS.params.RData')
load('NS.params.fixed.lambda.RData')
NS.params <- NS.params.fixed.lambda$Betas
load('hpi.RData')


mod.data <- mth.series.tbl[,c('date','UNRATE','NROU','NROUST',
                              'GDP','GDPC1','NGDPPOT','GDPPOT',
                              'PCETRIM12M159SFRBDAL','PCEPILFE','FEDFUNDS')]

mod.data <- mod.data %>%
              as_tibble

NS.params <- NS.params %>%
                fortify.zoo %>%
                  as_tibble

NS.params <- NS.params %>%
                mutate(Index = as.Date(Index)) %>%
                  rename(date = Index,
                         beta_0 = beta_1,
                         beta_1 = beta_2,
                         beta_2 = beta_3)

mod.data <- mod.data %>%
              full_join(NS.params, by='date')

hpi <- hpi %>%
        dplyr::select(date, hpi.sa) %>%
          as_tibble

mod.data <- mod.data %>%
              full_join(hpi, by='date')

mod.data <- mod.data %>%
              mutate(hpa.sa.1m = log(hpi.sa) - lag(log(hpi.sa)),
                     hpa.sa.12m = log(hpi.sa) - lag(log(hpi.sa), n=12),
                     beta_0.diff = beta_0 - lag(beta_0),
                     beta_1.diff = beta_1 - lag(beta_1),
                     beta_2.diff = beta_2 - lag(beta_2),
                     UNRATE.diff = UNRATE - lag(UNRATE),
                     UNRATE.ST.gap = UNRATE - NROUST,
                     UNRATE.LT.gap = UNRATE - NROU,
                     UNRATE.LT.ratio = UNRATE/NROU,
                     OUTPUT.GAP.nominal = log(NGDPPOT) - log(GDP),
                     OUTPUT.GAP.real = log(GDPPOT) - log(GDPC1)
                     )

mod.data <- mod.data %>%
              mutate(hpa.sa.gap = hpa.sa.1m - hpa.sa.12m)


mod.vars <- c('beta_0','beta_1','beta_2','UNRATE')

mod.data.sub <- mod.data %>%
                  dplyr::select(c('date',mod.vars)) %>%
                  filter(date>'1982-01-01') 
#                  fill(PCEPILFE)

mod.data.matrix <- as.data.frame(mod.data.sub[,mod.vars])

mn <- bv_minnesota(lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
                   alpha = bv_alpha(mode = 2, sd = 0.25, min = 1, max = 3),
                   var = 1e07)

soc <- bv_soc(mode = 1, sd = 1, min = 1e-04, max = 50)
sur <- bv_sur(mode = 1, sd = 1, min = 1e-04, max = 50)

priors <- bv_priors(hyper = "auto", mn = mn, sur = sur, soc = soc)

mh <- bv_metropolis(scale_hess = c(0.05, 0.0001, 0.0001),
                    adjust_acc = TRUE, acc_lower = 0.25, acc_upper = 0.45)

irfs <- bv_irf(horizon = 12, fevd = TRUE, identification = TRUE)

mod.bvar <- bvar(mod.data.matrix, lags = 1, n_draw = 50000, n_burn = 5000, n_thin = 1,
                 priors = priors, mh = mh, irf = irfs, verbose = TRUE)

print(summary(mod.bvar))

coef.mat <- summary(mod.bvar)$coef
dimnames(coef.mat)[[1]][2:4] <- dimnames(coef.mat)[[2]]

plot(coef.mat, digits=4, text.cell=list(cex=2.0),
            breaks=2,
            col=c('white','light blue'), key=NULL,
            xlab='Variable', ylab='Lagged Variable',
            main='Estimated VAR(1) Coefficients')

coef.mod.bvar <- coef(mod.bvar, conf_bands=c(0.05, 0.5, 0.95))

min.coef.mod.bvar <- as_tibble(coef.mod.bvar[1,,]) %>%
                      mutate(Variable=dimnames(coef.mod.bvar)[[2]]) %>%
                        pivot_longer(-Variable,
                                     names_to='Response',
                                     values_to='5%') %>%
                        select(Response, Variable, '5%')

mid.coef.mod.bvar <- as_tibble(coef.mod.bvar[2,,]) %>%
                      mutate(Variable=dimnames(coef.mod.bvar)[[2]]) %>%
                        pivot_longer(-Variable,
                                      names_to='Response',
                                      values_to='50%') %>%
                        select(Response, Variable, '50%')

max.coef.mod.bvar <- as_tibble(coef.mod.bvar[3,,]) %>%
                        mutate(Variable=dimnames(coef.mod.bvar)[[2]]) %>%
                        pivot_longer(-Variable,
                                      names_to='Response',
                                      values_to='95%') %>%
                                      select(Response, Variable, '95%')

conf.coef.mod.bvar <- min.coef.mod.bvar %>%
                        full_join(mid.coef.mod.bvar, by=c('Response','Variable')) %>%
                        full_join(max.coef.mod.bvar, by=c('Response','Variable')) %>%
                        arrange(Response)

print(conf.coef.mod.bvar, n=(length(mod.vars)*length(mod.vars)+length(mod.vars)))                    

print(plot(mod.bvar))

print(plot(irf(mod.bvar)))

plot(predict(mod.bvar, horizon=360))

pred.bvar <- predict(mod.bvar, horizon=360)

med.pred.bvar <- as.data.frame(pred.bvar$quants["50%",,])

colnames(med.pred.bvar) <- mod.vars

maturities <- c(3, 6, 12, 24, 36, 60, 84, 120)

lambda <- 0.0526069976895746

rates <- NS.rates(med.pred.bvar[c('beta_0','beta_1','beta_2')], lambda, maturities)

fcst.bvar <- cbind(as.data.frame(rates),
                   med.pred.bvar[,!(names(med.pred.bvar) %in% c('beta_1','beta_2','beta_3'))])

write.csv(fcst.bvar, file='fcst.bvar.csv')


# Simulate rates

rand.sim <- sample(1:45000)

x120.sim <- data.frame(sim=as.integer(),
                       month=as.integer(),
                       x120=as.numeric())

for (i in 1:1000) {
  
  betas <- pred.bvar$fcast[rand.sim[i],,]
  
  rates.sim <- NS.rates(betas, lambda, maturities)
  
  x120.sim <- rbind(x120.sim,
                    data.frame(sim=rep(i,360),
                               month=1:360,
                               x120=rates.sim[,'X120']))
  
}

x120.sim$sim <- as.factor(x120.sim$sim)

x120.quant <- quantile(x120.sim$x120[x120.sim$month==360], probs=c(0.10,0.90))

x120.outlier.sims <- x120.sim$sim[(x120.sim$x120<x120.quant[1] |
                                    x120.sim$x120>x120.quant[2]) &
                                    x120.sim$month==360]

x120.sim <- x120.sim %>%
              filter(!(sim %in% x120.outlier.sims))

rand.sim <- sample(as.integer(unique(x120.sim$sim)), 100)

x120.sim <- x120.sim %>%
              filter(sim %in% rand.sim)

x120.med <- data.frame(sim=as.factor(1),
                       month=1:360,
                       x120=rates[,"X120"])

print(ggplot(x120.sim, mapping=aes(x=month, y=x120, color=sim)) + 
               geom_line() +
               geom_line(x120.med, mapping=aes(x=month, y=x120), color='black', lwd=1.2) +
               ylab('Yield') +
               labs(title='10-Year Treasury Yield Simulations',
               subtitle='Bayesian VAR(1) forecasts of Nelson-Siegel parameters plus unemployment rate') + 
               theme(axis.title.x=element_blank(),
                     legend.position = "none"))

