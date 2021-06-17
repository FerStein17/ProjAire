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
library(egcm)
#install.packages("fGarch")
library(fGarch)



######################## Datos #################################
rm(list=ls())
install.packages("aire.zmvm")
library(aire.zmvm)

## Auto-install required R packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(aire.zmvm, dplyr, ggplot2, ggseas)

## download pollution data by station in ppb
pm10 <- get_station_data(criterion = "MAXIMOS", # Can be one of MAXIMOS (daily maximum), 
                       # MINIMOS (daily minimum), 
                       # or HORARIOS (hourly average)
                       pollutant = "PM10", 
                       year = 2011:2021) # 

pm10_max <- pm10 %>%
  group_by(date) %>%
  summarise(max = ifelse(all(is.na(value)),
                         NA,
                         base::max(value, na.rm = TRUE))) %>%
  na.omit()

pm10_max_ts<-tsdf(ts(pm10_max$max, start = c(2011,1), frequency = 365.25))

ggplot(pm10_max_ts,
       aes(x = x, y = y)) +
  geom_line(colour = "grey75", alpha = .5) +
  stat_rollapplyr(width = 30, align = "right", color = "#01C5D2") +
  xlab("date") +
  ylab("µg/m³") +
  scale_x_continuous(breaks = c(2013, 2016, 2019)) +
  ggtitle("PM10",
          "\nData source: SEDEMA") +
  theme_bw()




########################### Modelo ARIMA ####################################

PM10<-ts(pm10_max$max, start = c(2011,1), frequency = 365.25)
plot.ts(PM10,ylab = "µg/m³ diarias", main = "PM10")


#raices unitarias: H0 es que sí hay raíces unitarias y es estocástica
adf.test(PM10)


PM10<-PM10[1:3750]
PM10<-ts(PM10, start = c(2011,1), frequency = 365.25)
difrems<-diff((PM10))
difrem<-ts(difrems, start = c(2011,1), frequency = 365.25)


# ¿Son WN? NO
acf(difrem,length(difrem)^.5,na.action = na.pass)
pacf(difrem,length(difrem)^.5,na.action = na.pass)
Box.test(difrem, lag =length(difrem)^.5 , type = c("Box-Pierce"))
Box.test(difrem, lag = length(difrem)^.5, type = c("Ljung-Box"))



#Arima (p,0,q) 

ar1<-arima(difrem,order =c(1,0,0),n.cond=4)
ar2<-arima(difrem,order =c(2,0,0),n.cond=4)
ar3<-arima(difrem,order =c(3,0,0),n.cond=4)
ar4<-arima(difrem,order =c(4,0,0),n.cond=4)

ma1<-arima(difrem,order =c(0,0,1),n.cond=4)
ma2<-arima(difrem,order =c(0,0,2),n.cond=4)
ma3<-arima(difrem,order =c(0,0,3),n.cond=4)
ma4<-arima(difrem,order =c(0,0,4),n.cond=4)

arma101<-arima(difrem,order =c(1,0,1),n.cond=4)
arma102<-arima(difrem,order =c(1,0,2),n.cond=4)
arma103<-arima(difrem,order =c(1,0,3),n.cond=4)
arma104<-arima(difrem,order =c(1,0,4),n.cond=4)

arma201<-arima(difrem,order =c(2,0,1),n.cond=4)
arma202<-arima(difrem,order =c(2,0,2),n.cond=4)
arma203<-arima(difrem,order =c(2,0,3),n.cond=4)
arma204<-arima(difrem,order =c(2,0,4),n.cond=4)

arma301<-arima(difrem,order =c(3,0,1),n.cond=4)
arma302<-arima(difrem,order =c(3,0,2),n.cond=4)
arma303<-arima(difrem,order =c(3,0,3),n.cond=4)
arma304<-arima(difrem,order =c(3,0,4),n.cond=4)

arma401<-arima(difrem,order =c(4,0,1),n.cond=4)
arma402<-arima(difrem,order =c(4,0,2),n.cond=4)
arma403<-arima(difrem,order =c(4,0,3),n.cond=4)
arma404<-arima(difrem,order =c(4,0,4),n.cond=4)

# BIC
BIC(ar1)
BIC(ar2)
BIC(ar3)
BIC(ar4)

BIC(ma1)
BIC(ma2)
BIC(ma3)
BIC(ma4)

BIC(arma101)
BIC(arma102)
BIC(arma103)
BIC(arma104)

BIC(arma201)
BIC(arma202)
BIC(arma203)
BIC(arma204)

BIC(arma301)
BIC(arma302)
BIC(arma303)
BIC(arma304)

BIC(arma401)
BIC(arma402)
BIC(arma403)
BIC(arma404)


# Gana ARIMA 1,1
arma101<-arima(difrem,order =c(1,0,1))

arroot<-arroots(arma101)
print(arroot)
iarroot<-arroot$roots^-1
print(iarroot)
plot.armaroots(arroot)

maroot<-maroots(arma101)
print(maroot)
imaroot<-maroot$roots^-1
print(imaroot)
plot.armaroots(maroot)


# Residuales
residualsarma<-residuals(arma101)
plot.ts(residualsarma)
acf(residualsarma)
pacf(residualsarma)
Box.test(residualsarma, lag =length(residualsarma)^.5 , type = c("Box-Pierce"))
Box.test(residualsarma, lag =length(residualsarma)^.5 , type = c("Ljung-Box"))




#pronostico arma
predarma<-predict(arma101,n.ahead=10)
predarma

# 3751 es Jueves
summary.lm(regpm10,robust=TRUE)
tiempopred<-c(3751:3760)
tiempopred2<-tiempopred^2
predtend<-2.120e+02 - 3.695e-02*tiempopred + 6.243e-06*tiempopred2
predest<-c(8.693e+00, 7.175e+00, -3.466e+00, -2.884e+01, 0, 8.517e+00, 6.620e+00, 8.693e+00, 7.175e+00, -3.466e+00)


pronostico<-predarma$pred+predtend+predest
pronostico







garch10 <- garchFit(residualsarma~garch(1,0), data = residualsarma)
garch20 <- garchFit(residualsarma~garch(2,0), data = residualsarma)
garch30 <- garchFit(residualsarma~garch(3,0), data = residualsarma)
garch40 <- garchFit(residualsarma~garch(4,0), data = residualsarma)
garch50 <- garchFit(residualsarma~garch(5,0), data = residualsarma)
garch60 <- garchFit(residualsarma~garch(6,0), data = residualsarma)
garch70 <- garchFit(residualsarma~garch(7,0), data = residualsarma)

garch11 <- garchFit(residualsarma~garch(1,1), data = residualsarma)
garch22 <- garchFit(residualsarma~garch(2,2), data = residualsarma)
garch33 <- garchFit(residualsarma~garch(3,3), data = residualsarma)
garch44 <- garchFit(residualsarma~garch(4,4), data = residualsarma)
garch55 <- garchFit(residualsarma~garch(5,5), data = residualsarma)
garch66 <- garchFit(residualsarma~garch(6,6), data = residualsarma)
garch77 <- garchFit(residualsarma~garch(7,7), data = residualsarma)

garch12 <- garchFit(residualsarma~garch(1,2), data = residualsarma)
garch13 <- garchFit(residualsarma~garch(1,3), data = residualsarma)
garch14 <- garchFit(residualsarma~garch(1,4), data = residualsarma)
garch15 <- garchFit(residualsarma~garch(1,5), data = residualsarma)
garch16 <- garchFit(residualsarma~garch(1,6), data = residualsarma)
garch17 <- garchFit(residualsarma~garch(1,7), data = residualsarma)

garch21 <- garchFit(residualsarma~garch(2,1), data = residualsarma)
garch23 <- garchFit(residualsarma~garch(2,3), data = residualsarma)
garch24 <- garchFit(residualsarma~garch(2,4), data = residualsarma)
garch25 <- garchFit(residualsarma~garch(2,5), data = residualsarma)
garch26 <- garchFit(residualsarma~garch(2,6), data = residualsarma)
garch27 <- garchFit(residualsarma~garch(2,7), data = residualsarma)

garch31 <- garchFit(residualsarma~garch(3,1), data = residualsarma)
garch32 <- garchFit(residualsarma~garch(3,2), data = residualsarma)
garch34 <- garchFit(residualsarma~garch(3,4), data = residualsarma)
garch35 <- garchFit(residualsarma~garch(3,5), data = residualsarma)
garch36 <- garchFit(residualsarma~garch(3,6), data = residualsarma)
garch37 <- garchFit(residualsarma~garch(3,7), data = residualsarma)

garch41 <- garchFit(residualsarma~garch(4,1), data = residualsarma)
garch42 <- garchFit(residualsarma~garch(4,2), data = residualsarma)
garch43 <- garchFit(residualsarma~garch(4,3), data = residualsarma)
garch45 <- garchFit(residualsarma~garch(4,5), data = residualsarma)
garch46 <- garchFit(residualsarma~garch(4,6), data = residualsarma)
garch47 <- garchFit(residualsarma~garch(4,7), data = residualsarma)

garch51 <- garchFit(residualsarma~garch(5,1), data = residualsarma)
garch52 <- garchFit(residualsarma~garch(5,2), data = residualsarma)
garch53 <- garchFit(residualsarma~garch(5,3), data = residualsarma)
garch54 <- garchFit(residualsarma~garch(5,4), data = residualsarma)
garch56 <- garchFit(residualsarma~garch(5,6), data = residualsarma)
garch57 <- garchFit(residualsarma~garch(5,7), data = residualsarma)

garch61 <- garchFit(residualsarma~garch(6,1), data = residualsarma)
garch62 <- garchFit(residualsarma~garch(6,2), data = residualsarma)
garch63 <- garchFit(residualsarma~garch(6,3), data = residualsarma)
garch64 <- garchFit(residualsarma~garch(6,4), data = residualsarma)
garch65 <- garchFit(residualsarma~garch(6,5), data = residualsarma)
garch67 <- garchFit(residualsarma~garch(6,7), data = residualsarma)

garch71 <- garchFit(residualsarma~garch(7,1), data = residualsarma)
garch72 <- garchFit(residualsarma~garch(7,2), data = residualsarma)
garch73 <- garchFit(residualsarma~garch(7,3), data = residualsarma)
garch74 <- garchFit(residualsarma~garch(7,4), data = residualsarma)
garch75 <- garchFit(residualsarma~garch(7,5), data = residualsarma)
garch76 <- garchFit(residualsarma~garch(7,6), data = residualsarma)

summary(garch10)
summary(garch20)
summary(garch30)
summary(garch40)
summary(garch50)
summary(garch60)
summary(garch70)

summary(garch11)
summary(garch22)
summary(garch33)
summary(garch44)
summary(garch55)
summary(garch66)
summary(garch77)

summary(garch12)
summary(garch13)
summary(garch14)
summary(garch15)
summary(garch16)
summary(garch17)

summary(garch21)
summary(garch23)
summary(garch24)
summary(garch25)
summary(garch26)
summary(garch27)

summary(garch31)
summary(garch32)
summary(garch34)
summary(garch35)
summary(garch36)
summary(garch37)

summary(garch41)
summary(garch42)
summary(garch43)
summary(garch45)
summary(garch46)
summary(garch47)

summary(garch51)
summary(garch52)
summary(garch53)
summary(garch54)
summary(garch56)
summary(garch57)

summary(garch61)
summary(garch62)
summary(garch63)
summary(garch64)
summary(garch65)
summary(garch67)

summary(garch71)
summary(garch72)
summary(garch73)
summary(garch74)
summary(garch75)
summary(garch76)


#el mejor es el garch52
garch52 <- garchFit(residualsarma~garch(5,2), data = residualsarma)
summary(garch52)
predictse<-predict(garch52)
predictse

modelo<-garchFit(resid~arma(1,1)+garch(5,2), data=residualsarma)
summary(modelo)

se<-predictse$standardDeviation[1:10]
upper<-pronostico+1.96*se
print(upper)

lower<-pronostico-1.96*se
print(lower)

pronosticofinal<-pronostico
upperfinal<-upper
lowerfinal<-lower

#hacerlo serie de tiempo
pronosticots<-ts(pronosticofinal, start = c(2011,3751), frequency = 365.25)


lowerts<-ts(lowerfinal, start = c(2011,3751), frequency = 365.25)
upperts<-ts(upperfinal, start = c(2011,3751), frequency = 365.25)


#grafica con intervalos
ts.plot(PM10,pronosticots,upperts,lowerts,col=1:4,xlim=c(2021,2021.4),ylab="PM10")

#rmse
PM10<-ts(pm10_max$max[1:3760], start = c(2011,1), frequency = 365.25)
plot.ts(PM10,main="remesas")
observaciones<-PM10[3751:3760]
observaciones
rmse(pronosticofinal,observaciones)

#68.26971















################### Predicción vs datos guardados######################

predicciones<-predict(arma101,n.ahead=10)
pronostico<-predicciones$pred
pronostico2<-(cumsum(pronostico))+(189)
pronosticots<-ts(pronostico2,start=c(2011,3751),frequency=365.25)

# Intervalos de pronostico 
upper<-pronostico2+1.96*predicciones$se
lower<-pronostico2-1.96*predicciones$se
lowerts<-ts(lower, start = c(2011,3751), frequency = 365.25)
upperts<-ts(upper, start = c(2011,3751), frequency = 365.25)
ts.plot(PM10,pronosticots,upperts,lowerts,col=1:4,xlim=c(2011,2022),ylab="PM10 PPB")
ts.plot(PM10,pronosticots,upperts,lowerts,col=1:4,xlim=c(2021,2022),ylab="PM10 PPB")


# RMSE
PM10<-ts(pm10_max$max, start = c(2011,1), frequency = 365.25)
plot.ts(PM10,main = "PM10")
observaciones<-PM10[3751:3760]
rmse(pronosticots,observaciones)

# RMSE = 75.01651




