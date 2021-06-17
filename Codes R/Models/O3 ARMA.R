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
library(tseries)
#install.packages("fGarch")
library(fGarch)

######################## Datos ###################################3
rm(list=ls())
install.packages("aire.zmvm")
library(aire.zmvm)

## Auto-install required R packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(aire.zmvm, dplyr, ggplot2, ggseas)

## download pollution data by station in ppb
O3 <- get_station_data(criterion = "MAXIMOS", # Can be one of MAXIMOS (daily maximum), 
                       # MINIMOS (daily minimum), 
                       # or HORARIOS (hourly average)
                       pollutant = "O3", 
                       year = 2011:2021) # 

O3_max <- O3 %>%
  group_by(date) %>%
  summarise(max = ifelse(all(is.na(value)),
                         NA,
                         base::max(value, na.rm = TRUE))) %>%
  na.omit()

O3_max_ts<-tsdf(ts(O3_max$max, start = c(2011,1), frequency = 365.25))

ggplot(O3_max_ts,
       aes(x = x, y = y)) +
  geom_line(colour = "grey75", alpha = .5) +
  stat_rollapplyr(width = 30, align = "right", color = "#01C5D2") +
  xlab("date") +
  ylab("parts per billion") +
  scale_x_continuous(breaks = c(2013, 2016, 2019)) +
  ggtitle("Maximum daily ozone concentration and 30 day rolling average",
          "\nData source: SEDEMA") +
  theme_bw()




########################### Modelo ARMA ####################################

O3<-ts(O3_max$max[1:3760], start = c(2011,1), frequency = 365.25)
print(O3)
length(O3)
plot.ts(O3,main = "O3", ylab="PPB")
# Omitiré un año 14 días
O3<-O3[1:3750]
O3<-ts(O3, start = c(2011,1), frequency = 365.25)
length(O3)


# ¿Son white noise?   NO
acf(O3,length(O3)^.5,na.action = na.pass)
pacf(O3,length(O3)^.5,na.action = na.pass)
Box.test(O3, lag =length(O3)^.5 , type = c("Box-Pierce"))
Box.test(O3, lag =length(O3)^.5, type = c("Ljung-Box"))


# ¿Tendencia? NO
time<-(1:3750)
time2<-time^2

# Lineal
reg1<-lm(O3~time,data=O3)
summary.lm(reg1,robust=TRUE)
resid1<-residuals(reg1)
plot.ts(resid1,type="l",main="residuales")
trend1<-fitted.values(reg1)
ts.plot(O3, trend1, gpars = list(col = c("black", "blue")))

# Cuadrática
reg2<-lm(O3~time+time2,data=O3)
summary.lm(reg2,robust=TRUE)
resid2<-residuals(reg2)
plot.ts(resid2,type="l",main="residuales")
trend2<-fitted.values(reg2)
ts.plot(O3, trend2, gpars = list(col = c("black", "blue")))

# Exponencial
lnlogO3<-log(O3)
reg3<-lm(lnlogO3~time,data=lnlogO3)
summary.lm(reg3,robust=TRUE)
b0<-exp(4.553e+00)
b1<- 1.790e-06
reg4<-nls(O3~a*exp(b*time),start=list(a=b0,b=b1))
summary(reg4)
resid4<-residuals(reg4)
plot.ts(resid4,type="l",main="residuales")
trend3<-fitted.values(reg4)
ts.plot(O3, trend3, gpars = list(col = c("black", "blue")))

BIC(reg1)
BIC(reg2)
BIC(reg4)

# Sin tendencia


datos<- read_excel("/Users/dovabr/Desktop/ITAM/8ºSemestre/Macroeconometría/Proyecto/Datos temperatura.xlsx")
dd1<-ts(datos$d1, start = c(2011,1), frequency = 365.25)
dd2<-ts(datos$d2, start = c(2011,1), frequency = 365.25)
dd3<-ts(datos$d3, start = c(2011,1), frequency = 365.25)
dd4<-ts(datos$d4, start = c(2011,1), frequency = 365.25)
dd5<-ts(datos$d5, start = c(2011,1), frequency = 365.25)
dd6<-ts(datos$d6, start = c(2011,1), frequency = 365.25)
dd7<-ts(datos$d7, start = c(2011,1), frequency = 365.25)
covid<-ts(datos$covid, start = c(2011,1), frequency = 365.25)

dataO3<-ts.intersect(O3,dd1,dd2,dd3,dd4,dd5,dd6,dd7,covid)

regO3<-lm(O3~dd2+dd3+dd4+dd5+dd6+dd7+covid,data=dataO3)
summary.lm(regO3,robust=TRUE)

#wald
waldtest(regO3,reg2)

# Sí hay estacionalidad

# ¿White Noise? NO
summary.lm(regO3,robust=TRUE)
resid<-residuals(regO3)
plot.ts(resid,type="l",main="residuales")
acf(resid,length(resid)^.5,na.action = na.pass)
pacf(resid,length(resid)^.5,na.action = na.pass)
Box.test(resid, lag =length(resid)^.5 , type = c("Box-Pierce"))
Box.test(resid, lag = length(resid)^.5, type = c("Ljung-Box"))


#AR
#Arima (p,0,q) 

ar1<-arima(resid,order =c(1,0,0),n.cond=4)
ar2<-arima(resid,order =c(2,0,0),n.cond=4)
ar3<-arima(resid,order =c(3,0,0),n.cond=4)
ar4<-arima(resid,order =c(4,0,0),n.cond=4)

#ARMA(0,q) con el n.con para que sea comparable el numero de obs 
ma1<-arima(resid,order =c(0,0,1),n.cond=4)
ma2<-arima(resid,order =c(0,0,2),n.cond=4)
ma3<-arima(resid,order =c(0,0,3),n.cond=4)
ma4<-arima(resid,order =c(0,0,4),n.cond=4)

arma101<-arima(resid,order =c(1,0,1),n.cond=4)
arma102<-arima(resid,order =c(1,0,2),n.cond=4)
arma103<-arima(resid,order =c(1,0,3),n.cond=4)
arma104<-arima(resid,order =c(1,0,4),n.cond=4)

arma201<-arima(resid,order =c(2,0,1),n.cond=4)
arma202<-arima(resid,order =c(2,0,2),n.cond=4)
arma203<-arima(resid,order =c(2,0,3),n.cond=4)
arma204<-arima(resid,order =c(2,0,4),n.cond=4)

arma301<-arima(resid,order =c(3,0,1),n.cond=4)
arma302<-arima(resid,order =c(3,0,2),n.cond=4)
arma303<-arima(resid,order =c(3,0,3),n.cond=4)
arma304<-arima(resid,order =c(3,0,4),n.cond=4)

arma401<-arima(resid,order =c(4,0,1),n.cond=4)
arma402<-arima(resid,order =c(4,0,2),n.cond=4)
arma403<-arima(resid,order =c(4,0,3),n.cond=4)
arma404<-arima(resid,order =c(4,0,4),n.cond=4)

#BIC
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


#Gana ARMA 2,1
arma201<-arima(resid,order =c(2,0,1))
arma201

#nos quedamos con el arma 2,1
arroot<-arroots(arma201)
print(arroot)
iarroot<-arroot$roots^-1
print(iarroot)
plot.armaroots(arroot)

maroot<-maroots(arma201)
print(maroot)
imaroot<-maroot$roots^-1
print(imaroot)
plot.armaroots(maroot)

# ¿Son White Noise? Sí
residualesarma<-residuals(arma201)
plot.ts(residualesarma)
acf(residualesarma)
pacf(residualesarma)
Box.test(residualesarma, lag =length(residualesarma)^.5 , type = c("Box-Pierce"))
Box.test(residualesarma, lag =length(residualesarma)^.5 , type = c("Ljung-Box"))


fittedarma<-fitted.values(arma201)
ts.plot(resid, fittedarma, gpars = list(col = c("black", "blue")))


#pronostico arma
predarma<-predict(arma201,n.ahead=10)
predarma

# 3751 es Jueves
summary.lm(regO3,robust=TRUE)
predtend<- 97.124 
predest<-c( 3.414, 3.547,  3.889, 2.195, 0, 2.919, 2.698, 3.414, 3.547,  3.889)

pronostico<-predarma$pred+predtend+predest
pronostico

#Intervalos de pronostico 
upper<-pronostico+1.96*predarma$se
print(upper)

lower<-pronostico-1.96*predarma$se
print(lower)


#hacerlo serie de tiempo
pronosticots<-ts(pronostico, start = c(2011,3751), frequency = 365.25)
lowerts<-ts(lower, start = c(2011,3751), frequency = 365.25)
upperts<-ts(upper, start = c(2011,3751), frequency = 365.25)


#grafica con intervalos
ts.plot(O3,pronosticots,upperts,lowerts,col=1:4,xlim=c(2021,2021.4),ylab="remesas")

#rmse
O3<-ts(O3_max$max[1:3760], start = c(2011,1), frequency = 365.25)
plot.ts(O3,main="O3")
observaciones<-O3[3751:3760]
observaciones
rmse(pronostico,observaciones)

# RMSE = 28.31494
















############### Predicción ####################3

O3<-ts(O3_max$max[1:3760], start = c(2011,1), frequency = 365.25)


# ¿Tendencia? 
time<-(1:3760)
time2<-time^2

datos<- read_excel("/Users/dovabr/Desktop/ITAM/8ºSemestre/Macroeconometría/Proyecto/Datos temperatura.xlsx")
dd1<-ts(datos$d1, start = c(2011,1), frequency = 365.25)
dd2<-ts(datos$d2, start = c(2011,1), frequency = 365.25)
dd3<-ts(datos$d3, start = c(2011,1), frequency = 365.25)
dd4<-ts(datos$d4, start = c(2011,1), frequency = 365.25)
dd5<-ts(datos$d5, start = c(2011,1), frequency = 365.25)
dd6<-ts(datos$d6, start = c(2011,1), frequency = 365.25)
dd7<-ts(datos$d7, start = c(2011,1), frequency = 365.25)
covid<-ts(datos$covid, start = c(2011,1), frequency = 365.25)

dataO3<-ts.intersect(O3,time,time2,dd1,dd2,dd3,dd4,dd5,dd6,dd7,covid)

regO3<-lm(O3~time+time2+dd2+dd3+dd4+dd5+dd6+dd7+covid,data=dataO3)
summary.lm(regO3,robust=TRUE)


# ¿White Noise? NO
summary.lm(regO3,robust=TRUE)
resid<-residuals(regO3)
plot.ts(resid,type="l",main="residuales")
acf(resid,length(resid)^.5,na.action = na.pass)
pacf(resid,length(resid)^.5,na.action = na.pass)
Box.test(resid, lag =length(resid)^.5 , type = c("Box-Pierce"))
Box.test(resid, lag = length(resid)^.5, type = c("Ljung-Box"))


#AR
#Arima (p,0,q) 

ar1<-arima(resid,order =c(1,0,0),n.cond=4)
ar2<-arima(resid,order =c(2,0,0),n.cond=4)
ar3<-arima(resid,order =c(3,0,0),n.cond=4)
ar4<-arima(resid,order =c(4,0,0),n.cond=4)

#ARMA(0,q) con el n.con para que sea comparable el numero de obs 
ma1<-arima(resid,order =c(0,0,1),n.cond=4)
ma2<-arima(resid,order =c(0,0,2),n.cond=4)
ma3<-arima(resid,order =c(0,0,3),n.cond=4)
ma4<-arima(resid,order =c(0,0,4),n.cond=4)

arma101<-arima(resid,order =c(1,0,1),n.cond=4)
arma102<-arima(resid,order =c(1,0,2),n.cond=4)
arma103<-arima(resid,order =c(1,0,3),n.cond=4)
arma104<-arima(resid,order =c(1,0,4),n.cond=4)

arma201<-arima(resid,order =c(2,0,1),n.cond=4)
arma202<-arima(resid,order =c(2,0,2),n.cond=4)
arma203<-arima(resid,order =c(2,0,3),n.cond=4)
arma204<-arima(resid,order =c(2,0,4),n.cond=4)

arma301<-arima(resid,order =c(3,0,1),n.cond=4)
arma302<-arima(resid,order =c(3,0,2),n.cond=4)
arma303<-arima(resid,order =c(3,0,3),n.cond=4)
arma304<-arima(resid,order =c(3,0,4),n.cond=4)

arma401<-arima(resid,order =c(4,0,1),n.cond=4)
arma402<-arima(resid,order =c(4,0,2),n.cond=4)
arma403<-arima(resid,order =c(4,0,3),n.cond=4)
arma404<-arima(resid,order =c(4,0,4),n.cond=4)

#BIC
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


#Gana ARMA 2,1
arma201<-arima(resid,order =c(2,0,1))
arma201

residualesarma<-residuals(arma201)
plot.ts(residualesarma)
acf(residualesarma)
pacf(residualesarma)
Box.test(residualesarma, lag =length(residualesarma)^.5 , type = c("Box-Pierce"))
Box.test(residualesarma, lag =length(residualesarma)^.5 , type = c("Ljung-Box"))


fittedarma<-fitted.values(arma201)
ts.plot(resid, fittedarma, gpars = list(col = c("black", "blue")))



#pronostico arma
predarma<-predict(arma201,n.ahead=10)
predarma

# 3761 es Sábado
summary.lm(regO3,robust=TRUE)
tiempopred<-c(3761:3770)
tiempopred2<-tiempopred^2
predtend<- 2.124e+02 - 3.803e-02*tiempopred + 6.608e-06*tiempopred2
predest<-c(-3.273e+00, -2.892e+01, 0, 8.860e+00, 6.654e+00, 8.609e+00, 7.274e+00, -3.273e+00, -2.892e+01, 0)


pronostico<-predarma$pred+predtend+predest
pronostico







garch10 <- garchFit(residualesarma~garch(1,0), data = residualesarma)
garch20 <- garchFit(residualesarma~garch(2,0), data = residualesarma)
garch30 <- garchFit(residualesarma~garch(3,0), data = residualesarma)
garch40 <- garchFit(residualesarma~garch(4,0), data = residualesarma)
garch50 <- garchFit(residualesarma~garch(5,0), data = residualesarma)
garch60 <- garchFit(residualesarma~garch(6,0), data = residualesarma)
garch70 <- garchFit(residualesarma~garch(7,0), data = residualesarma)

garch11 <- garchFit(residualesarma~garch(1,1), data = residualesarma)
garch22 <- garchFit(residualesarma~garch(2,2), data = residualesarma)
garch33 <- garchFit(residualesarma~garch(3,3), data = residualesarma)
garch44 <- garchFit(residualesarma~garch(4,4), data = residualesarma)
garch55 <- garchFit(residualesarma~garch(5,5), data = residualesarma)
garch66 <- garchFit(residualesarma~garch(6,6), data = residualesarma)
garch77 <- garchFit(residualesarma~garch(7,7), data = residualesarma)

garch12 <- garchFit(residualesarma~garch(1,2), data = residualesarma)
garch13 <- garchFit(residualesarma~garch(1,3), data = residualesarma)
garch14 <- garchFit(residualesarma~garch(1,4), data = residualesarma)
garch15 <- garchFit(residualesarma~garch(1,5), data = residualesarma)
garch16 <- garchFit(residualesarma~garch(1,6), data = residualesarma)
garch17 <- garchFit(residualesarma~garch(1,7), data = residualesarma)

garch21 <- garchFit(residualesarma~garch(2,1), data = residualesarma)
garch23 <- garchFit(residualesarma~garch(2,3), data = residualesarma)
garch24 <- garchFit(residualesarma~garch(2,4), data = residualesarma)
garch25 <- garchFit(residualesarma~garch(2,5), data = residualesarma)
garch26 <- garchFit(residualesarma~garch(2,6), data = residualesarma)
garch27 <- garchFit(residualesarma~garch(2,7), data = residualesarma)

garch31 <- garchFit(residualesarma~garch(3,1), data = residualesarma)
garch32 <- garchFit(residualesarma~garch(3,2), data = residualesarma)
garch34 <- garchFit(residualesarma~garch(3,4), data = residualesarma)
garch35 <- garchFit(residualesarma~garch(3,5), data = residualesarma)
garch36 <- garchFit(residualesarma~garch(3,6), data = residualesarma)
garch37 <- garchFit(residualesarma~garch(3,7), data = residualesarma)

garch41 <- garchFit(residualesarma~garch(4,1), data = residualesarma)
garch42 <- garchFit(residualesarma~garch(4,2), data = residualesarma)
garch43 <- garchFit(residualesarma~garch(4,3), data = residualesarma)
garch45 <- garchFit(residualesarma~garch(4,5), data = residualesarma)
garch46 <- garchFit(residualesarma~garch(4,6), data = residualesarma)
garch47 <- garchFit(residualesarma~garch(4,7), data = residualesarma)

garch51 <- garchFit(residualesarma~garch(5,1), data = residualesarma)
garch52 <- garchFit(residualesarma~garch(5,2), data = residualesarma)
garch53 <- garchFit(residualesarma~garch(5,3), data = residualesarma)
garch54 <- garchFit(residualesarma~garch(5,4), data = residualesarma)
garch56 <- garchFit(residualesarma~garch(5,6), data = residualesarma)
garch57 <- garchFit(residualesarma~garch(5,7), data = residualesarma)

garch61 <- garchFit(residualesarma~garch(6,1), data = residualesarma)
garch62 <- garchFit(residualesarma~garch(6,2), data = residualesarma)
garch63 <- garchFit(residualesarma~garch(6,3), data = residualesarma)
garch64 <- garchFit(residualesarma~garch(6,4), data = residualesarma)
garch65 <- garchFit(residualesarma~garch(6,5), data = residualesarma)
garch67 <- garchFit(residualesarma~garch(6,7), data = residualesarma)

garch71 <- garchFit(residualesarma~garch(7,1), data = residualesarma)
garch72 <- garchFit(residualesarma~garch(7,2), data = residualesarma)
garch73 <- garchFit(residualesarma~garch(7,3), data = residualesarma)
garch74 <- garchFit(residualesarma~garch(7,4), data = residualesarma)
garch75 <- garchFit(residualesarma~garch(7,5), data = residualesarma)
garch76 <- garchFit(residualesarma~garch(7,6), data = residualesarma)

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


#el mejor es el garch74
garch74 <- garchFit(residualesarma~garch(7,4), data = residualesarma)
summary(garch74)
predictse<-predict(garch74)
predictse

modelo<-garchFit(resid~arma(2,1)+garch(7,4), data=resid)
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
pronosticots<-ts(pronosticofinal, start = c(2011,3761), frequency = 365.25)


lowerts<-ts(lowerfinal, start = c(2011,3761), frequency = 365.25)
upperts<-ts(upperfinal, start = c(2011,3761), frequency = 365.25)


#grafica con intervalos
ts.plot(O3,pronosticots,upperts,lowerts,col=1:4,xlim=c(2021,2021.4),ylab="remesas")



