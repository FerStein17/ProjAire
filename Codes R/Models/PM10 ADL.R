library(forecast)
library(timeDate)
library(timeSeries)
library(fBasics)
library(fUnitRoots)
library(readxl)
library(Metrics)
library(astsa)
library(xts)
library(tseries)
library(haven)
library(zoo)
library(lmtest)
library(ggplot2)
library(forecast)
library(readxl)
install.packages("dynlm")
library(dynlm)
library(sandwich)

######################## Datos ###################################3
rm(list=ls())
install.packages("aire.zmvm")
library(aire.zmvm)

## Auto-install required R packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(aire.zmvm, dplyr, ggplot2, ggseas)

## download pollution data by station in ppb
PM10 <- get_station_data(criterion = "MAXIMOS", # Can be one of MAXIMOS (daily maximum), 
                       # MINIMOS (daily minimum), 
                       # or HORARIOS (hourly average)
                       pollutant = "PM10", 
                       year = 2017:2021) # 

PM10_max <- PM10 %>%
  group_by(date) %>%
  summarise(max = ifelse(all(is.na(value)),
                         NA,
                         base::max(value, na.rm = TRUE))) %>%
  na.omit()

PM10_max_ts<-tsdf(ts(PM10_max$max, start = c(2017,1), frequency = 365.25))

ggplot(PM10_max_ts,
       aes(x = x, y = y)) +
  geom_line(colour = "grey75", alpha = .5) +
  stat_rollapplyr(width = 30, align = "right", color = "#01C5D2") +
  xlab("date") +
  ylab("parts per billion") +
  scale_x_continuous(breaks = c(2013, 2016, 2019)) +
  ggtitle("Maximum daily ozone concentration and 30 day rolling average",
          "\nData source: SEDEMA") +
  theme_bw()







################################## ADL ########################################3

PM10<-ts(PM10_max$max[284 : 1572], start = c(2017,284), frequency = 365.25)
print(PM10)
length(PM10)
plot.ts(PM10,main = "PM10", ylab="PM10 PPB")
length(PM10)


# ¿Son white noise?   NO
acf(PM10,length(PM10)^.5,na.action = na.pass)
pacf(PM10,length(PM10)^.5,na.action = na.pass)
Box.test(PM10, lag =length(PM10)^.5 , type = c("Box-Pierce"))
Box.test(PM10, lag =length(PM10)^.5, type = c("Ljung-Box"))




setwd("/Users/dovabr/Desktop/ITAM/8ºSemestre/Macroeconometría/Proyecto/Buenos")
traffic <- read.csv('Traffic_Data_CDMX.csv')
traffic$date<-as.Date(traffic$date,"%m/%d/%Y")
traf <- ts(traffic$Traffic, start = c(2017,284), frequency = 365.25)

ADLdata <- ts.union(PM10, traf)

# estimar los adl
ADL101 <- dynlm(PM10 ~ L(PM10) +  L(traf), start =  c(2017,284), end = c(2021, 97))
coeftest(ADL101)
ADL202 <- dynlm(PM10 ~ L(PM10) +  L(traf) + L(PM10,2) +  L(traf,2), start =  c(2017,284), end = c(2021, 97))
coeftest(ADL202)
ADL303 <- dynlm(PM10 ~ L(PM10) +  L(traf) + L(PM10,2) +  L(traf,2) + L(PM10,3) +  L(traf,3), start =  c(2017,284), end = c(2021, 97))
coeftest(ADL303)
ADL404 <- dynlm(PM10 ~ L(PM10) +  L(traf) + L(PM10,2) +  L(traf,2) + L(PM10,3) +  L(traf,3) + L(PM10,4) +  L(traf,4), start =  c(2017,284), end = c(2021, 97))
coeftest(ADL404)
ADL505 <- dynlm(PM10 ~ L(PM10) +  L(traf) + L(PM10,2) +  L(traf,2) + L(PM10,3) +  L(traf,3) + L(PM10,4) +  L(traf,4) +  L(PM10,5) +  L(traf,5), start =  c(2017,284), end = c(2021, 97))
coeftest(ADL505)
ADL606 <- dynlm(PM10 ~ L(PM10) +  L(traf) + L(PM10,2) +  L(traf,2) +  L(PM10,3) + L(traf,3) + L(PM10,4) +  L(traf,4) +  L(PM10,5) +  L(traf,5) + L(PM10,6) +  L(traf,6),start =  c(2017,284), end = c(2021, 97))
coeftest(ADL606)
ADL707 <- dynlm(PM10 ~ L(PM10) +  L(traf) + L(PM10,2) +  L(traf,2) +  L(PM10,3) + L(traf,3) + L(PM10,4) +  L(traf,4) +  L(PM10,5) +  L(traf,5) + L(PM10,6) +  L(traf,6) + L(PM10,7) +  L(traf,7), start =  c(2017,284), end = c(2021, 97))
coeftest(ADL707)
ADL808 <- dynlm(PM10 ~ L(PM10) +  L(traf) + L(PM10,2) +  L(traf,2) +  L(PM10,3) + L(traf,3) + L(PM10,4) +  L(traf,4) +  L(PM10,5) +  L(traf,5) + L(PM10,6) +  L(traf,6) + L(PM10,7) +  L(traf,7) + L(PM10,8) +  L(traf,8), start =  c(2017,284), end = c(2021, 97))
coeftest(ADL808)


BIC(ADL101)
BIC(ADL202)
BIC(ADL303)
BIC(ADL404)
BIC(ADL505)
BIC(ADL606)
BIC(ADL707)
BIC(ADL808)



#Gana ADL404

grangertest(PM10~traf,order=4,data=ADLdata)

# H0: No causa
# Tráfico  SÍ granger causa PM10


resid<-residuals(ADL404)
resid<- ts(resid, start = c(2017,284), frequency = 365.25)
plot.ts(resid)
acf(resid)
pacf(resid)
Box.test(resid, lag =length(resid)^.5 , type = c("Box-Pierce"))
Box.test(resid, lag =length(resid)^.5 , type = c("Ljung-Box"))





# Predicción de 18 Abril
coeftest(ADL404)

pm10abr14<-223
pm10abr15<-122
pm10abr16<-174
pm10abr17<-233

trafabr14<-29.70463
trafabr15<-24.37246
trafabr16<-22.81552
trafabr17<-41.93305


pm10abr18<- 52.8663043 + 0.2091581*pm10abr17 + 0.3237163*trafabr17 + 0.1859682*pm10abr16 - 0.0071147*trafabr16 + 0.0831186*pm10abr15 + 0.0376447*trafabr15 + 0.1820916*pm10abr14 - 0.2877793*trafabr14
pm10abr18 




pm10abr14<-223
pm10abr15<-122
pm10abr16<-174
pm10abr17<-233
pm10abr18<-261

trafabr14<-29.70463
trafabr15<-24.37246
trafabr16<-22.81552
trafabr17<-41.93305
trafabr18<-34.50577

pm10abr19<- 52.8663043 + 0.2091581*pm10abr18 + 0.3237163*trafabr18 + 0.1859682*pm10abr17 - 0.0071147*trafabr17 + 0.0831186*pm10abr16 + 0.0376447*trafabr16 + 0.1820916*pm10abr15 - 0.2877793*trafabr15
pm10abr19 


