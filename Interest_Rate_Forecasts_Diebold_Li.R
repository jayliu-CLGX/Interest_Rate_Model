
setwd("C:\\R Working Directory\\HPIF_Scenarios")

library(dplyr)
library(tidyr)
library(xts)
library(xtsExtra)
library(YieldCurve)
library(forecast)
library(tibble)
library(lubridate)
library(forcats)
library(ggplot2)
library(stringr)


treas.yields <- read.table('Treasury_CMT_Yields.csv', header=TRUE,
                      stringsAsFactors=FALSE, sep=",", fileEncoding="UTF-8-BOM")
treas.yields$Month <- as.Date(treas.yields$Month, '%Y-%m-%d')
rownames(treas.yields) <- treas.yields$Month


# Drop 30-year yields
treas.yields <- treas.yields %>%
                  dplyr::select(-GS30, -Month)

treas.yields <- treas.yields[complete.cases(treas.yields),]

treas.yields <- as.xts(treas.yields)


# Maturities
maturaties <- c(3, 6, 12, 24, 36, 60, 84, 120)


# Fit dynamic Nelson-Siegel model
#NS.params <- Nelson.Siegel(treas.yields, maturaties)
#save(NS.params, file='NS.params.RData')
load('NS.params.RData')


# Fit AR(1) models to the Nelson-Siegel factors
mod_beta_0 <- arima(NS.params$beta_0, order=c(1, 0, 0))
mod_beta_1 <- arima(NS.params$beta_1, order=c(1, 0, 0))
mod_beta_2 <- arima(NS.params$beta_2, order=c(1, 1, 0))


# Forecast the factors
fcst_beta_0 <- predict(mod_beta_0, n.ahead=360)[['pred']]
fcst_beta_1 <- predict(mod_beta_1, n.ahead=360)[['pred']]
fcst_beta_2 <- predict(mod_beta_2, n.ahead=360)[['pred']]


fcst.months <- seq(as.Date("2020-05-01", "%Y-%m-%d"), by="month", length.out=360)

NS.params.fcst <- data.frame(beta_0=as.numeric(fcst_beta_0),
                             beta_1=as.numeric(fcst_beta_1),
                             beta_2=as.numeric(fcst_beta_2),
                             lambda=rep(0.06,length(fcst.months)),
                             row.names=fcst.months)
NS.params.fcst <- as.xts(NS.params.fcst)


# Generate yield forecasts
maturaties <- c(3, 6, 12, 24, 36, 60, 84, 120, 360)

treas.yields.fcst <- NSrates(NS.params.fcst, maturaties)

treas.yields.fcst.df <- fortify(treas.yields.fcst)
colnames(treas.yields.fcst.df)[1] <- 'month'


treas.yields.fcst.df <- treas.yields.fcst.df %>%
                          pivot_longer(cols=starts_with('X'), 
                                       names_to='maturity', 
                                       values_to='yield')

treas.yields.fcst.df <- treas.yields.fcst.df %>%
                         filter(maturity!='X360') %>%
                          mutate(maturity=as.integer(str_remove(maturity,'X'))) %>%
                           mutate(label=as.factor(paste(maturity,'months')))

treas.yields.fcst.df$label <- factor(treas.yields.fcst.df$label, 
                                     levels=unique(treas.yields.fcst.df$label))
treas.yields.fcst.df$label <- fct_rev(treas.yields.fcst.df$label)

print(ggplot(treas.yields.fcst.df %>% arrange(maturity), 
             aes(x=month, y=yield, colour=label)) +
             ylab('Yield') +
             labs(title='Treasury Yield Forecasts, April 2020',
                  subtitle='Separate AR(1) forecasts of Nelson-Siegel parameters') + 
             theme(axis.title.x=element_blank(),
                   legend.title=element_blank()) +
             geom_line())

plot.xts(treas.yields.fcst, screens=1, auto.legend=TRUE)


# Keep yields and spreads from becoming negative
treas.yields.fcst$X3[treas.yields.fcst$X3 < 0.15] <- 0.15
treas.yields.fcst$X6[(treas.yields.fcst$X6 - treas.yields.fcst$X3) < 0.02] <- 
  treas.yields.fcst$X3[(treas.yields.fcst$X6 - treas.yields.fcst$X3) < 0.02] + 0.02
treas.yields.fcst$X12[(treas.yields.fcst$X12 - treas.yields.fcst$X6) < 0.02] <- 
  treas.yields.fcst$X6[(treas.yields.fcst$X12 - treas.yields.fcst$X6) < 0.02] + 0.02
treas.yields.fcst$X24[(treas.yields.fcst$X24 - treas.yields.fcst$X12) < 0.02] <- 
  treas.yields.fcst$X12[(treas.yields.fcst$X24 - treas.yields.fcst$X12) < 0.02] + 0.02

treas.yields.fcst.baseline <- treas.yields.fcst

treas.yields.fcst.baseline$MTG30 <- 1.5467 +
                                    1.0336 * treas.yields.fcst.baseline$X120

plot.xts(treas.yields.fcst.baseline, screens=1, auto.legend=TRUE)


# Extend Dec-2020 forecasts to generate scenarios

adverse.beg <- treas.yields.fcst.baseline["2020"]

adverse.mid <- as_tibble(treas.yields.fcst.baseline["2020-12"]) %>%
                slice(rep(row_number(), 12))
adverse.mid.rownames <- index(treas.yields.fcst.baseline["2021"])
adverse.mid <- as.data.frame(adverse.mid)
rownames(adverse.mid) <- adverse.mid.rownames
adverse.mid <- as.xts(adverse.mid)

adverse.end <- treas.yields.fcst.baseline["2021/2050-04"]
index(adverse.end) <- index(adverse.end) + months(12)

treas.yields.fcst.adverse <- do.call("rbind", list(adverse.beg, adverse.mid, adverse.end))
 
treas.yields.fcst.adverse$MTG30 <- 1.5467 +
  1.0336 * treas.yields.fcst.adverse$X120

plot.xts(treas.yields.fcst.adverse, screens=1, auto.legend=TRUE)


severe.beg <- treas.yields.fcst.baseline["2020"]

severe.mid <- as_tibble(treas.yields.fcst.baseline["2020-12"]) %>%
  slice(rep(row_number(), 36))
severe.mid.rownames <- index(treas.yields.fcst.baseline["2021/2023"])
severe.mid <- as.data.frame(severe.mid)
rownames(severe.mid) <- severe.mid.rownames
severe.mid <- as.xts(severe.mid)

severe.end <- treas.yields.fcst.baseline["2021/2050-04"]
index(severe.end) <- index(severe.end) + months(36)

treas.yields.fcst.severe <- do.call("rbind", list(severe.beg, severe.mid, severe.end))

treas.yields.fcst.severe$MTG30 <- 1.5467 +
  1.0336 * treas.yields.fcst.severe$X120

plot.xts(treas.yields.fcst.severe, screens=1, auto.legend=TRUE)

treas.yields.fcst.baseline$Month <- index(treas.yields.fcst.baseline)

write.csv(treas.yields.fcst.baseline, file='treas.yields.fcst.baseline.csv')
write.csv(treas.yields.fcst.adverse, file='treas.yields.fcst.adverse.csv')
write.csv(treas.yields.fcst.severe, file='treas.yields.fcst.severe.csv')