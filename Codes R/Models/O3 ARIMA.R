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


########################### Modelo ARIMA ####################################

O3<-ts(O3_max$max, start = c(2011,1), frequency = 365.25)
print(O3)
length(O3)
plot.ts(O3,main = "PPB diarias")

#raices unitarias: H0 es que sí hay raíces unitarias y es estocástica
adf.test(O3)


O3<-O3[1:3750]
O3<-ts(O3, start = c(2011,1), frequency = 365.25)
length(O3)

difo3<-diff((O3))
print(difo3)
length(difo3)
plot.ts(difo3,main="Dif O3 PPB")
difrem<-ts(difo3, start = c(2011,1), frequency = 365.25)

#checamos si es WN
acf(difrem,length(difrem)^.5,na.action = na.pass)
pacf(difrem,length(difrem)^.5,na.action = na.pass)
Box.test(difrem, lag =length(difrem)^.5 , type = c("Box-Pierce"))
Box.test(difrem, lag = length(difrem)^.5, type = c("Ljung-Box"))


#AR
#Arima (p,0,q) 

ar1<-arima(difrem,order =c(1,0,0),n.cond=4)
ar2<-arima(difrem,order =c(2,0,0),n.cond=4)
ar3<-arima(difrem,order =c(3,0,0),n.cond=4)
ar4<-arima(difrem,order =c(4,0,0),n.cond=4)

#ARMA(0,q) con el n.con para que sea comparable el numero de obs 
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

# Gana ARIMA 1,1
arma101<-arima(difrem,order =c(1,0,1))
arma101

#nos quedamos con el arma101
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


#checamos sus residuales
residualsarma<-residuals(arma101)
plot.ts(residualsarma)
acf(residualsarma)
pacf(residualsarma)
Box.test(residualsarma, lag =length(residualsarma)^.5 , type = c("Box-Pierce"))
Box.test(residualsarma, lag =length(residualsarma)^.5 , type = c("Ljung-Box"))





predicciones<-predict(arma101,n.ahead=10)
pronostico<-predicciones$pred
pronostico
pronostico2<-(cumsum(pronostico))+(92)
pronosticots<-ts(pronostico2,start=c(2011,3751),frequency=365.25)
pronosticots

#Intervalos de pronostico 
upper<-pronostico2+1.96*predicciones$se
print(upper)
lower<-pronostico2-1.96*predicciones$se
print(lower)
#hacerlo serie de tiempo
lowerts<-ts(lower, start = c(2011,3751), frequency = 365.25)
upperts<-ts(upper, start = c(2011,3751), frequency = 365.25)


#grafica con intervalos
ts.plot(O3,pronosticots,upperts,lowerts,col=1:4,xlim=c(2011,2022),ylab="O3 PPB")
ts.plot(O3,pronosticots,upperts,lowerts,col=1:4,xlim=c(2021,2022),ylab="O3 PPB")

#error#rmse
O3<-ts(O3_max$max[1:3760], start = c(2011,1), frequency = 365.25)
length(O3)
plot.ts(O3,main="O3 PPB")
observaciones<-O3[3751:3760]
rmse(pronosticots,observaciones)


# RMSE =  27.2247
# GANA
















O3<-ts(O3_max$max, start = c(2011,1), frequency = 365.25)
O3<-O3[1:3760]
O3<-ts(O3, start = c(2011,1), frequency = 365.25)


difo3<-diff((O3))
print(difo3)
length(difo3)
plot.ts(difo3,main="Dif O3 PPB")
difrem<-ts(difo3, start = c(2011,1), frequency = 365.25)

#checamos si es WN
acf(difrem,length(difrem)^.5,na.action = na.pass)
pacf(difrem,length(difrem)^.5,na.action = na.pass)
Box.test(difrem, lag =length(difrem)^.5 , type = c("Box-Pierce"))
Box.test(difrem, lag = length(difrem)^.5, type = c("Ljung-Box"))


#AR
#Arima (p,0,q) 

ar1<-arima(difrem,order =c(1,0,0),n.cond=4)
ar2<-arima(difrem,order =c(2,0,0),n.cond=4)
ar3<-arima(difrem,order =c(3,0,0),n.cond=4)
ar4<-arima(difrem,order =c(4,0,0),n.cond=4)

#ARMA(0,q) con el n.con para que sea comparable el numero de obs 
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

# Gana ARIMA 1,1
arma101<-arima(difrem,order =c(1,0,1))
arma101

#nos quedamos con el arma101
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


#checamos sus residuales
residualsarma<-residuals(arma101)
plot.ts(residualsarma)
acf(residualsarma)
pacf(residualsarma)
Box.test(residualsarma, lag =length(residualsarma)^.5 , type = c("Box-Pierce"))
Box.test(residualsarma, lag =length(residualsarma)^.5 , type = c("Ljung-Box"))




predicciones<-predict(arma101,n.ahead=10)
pronostico<-predicciones$pred
pronostico
pronostico2<-(cumsum(pronostico))+(129)
pronosticots<-ts(pronostico2,start=c(2011,3761),frequency=365.25)
pronosticots

#Intervalos de pronostico 
upper<-pronostico2+1.96*predicciones$se
print(upper)
lower<-pronostico2-1.96*predicciones$se
print(lower)
#hacerlo serie de tiempo
lowerts<-ts(lower, start = c(2011,3761), frequency = 365.25)
upperts<-ts(upper, start = c(2011,3761), frequency = 365.25)


#grafica con intervalos
ts.plot(O3,pronosticots,upperts,lowerts,col=1:4,xlim=c(2011,2022),ylab="O3 PPB")
ts.plot(O3,pronosticots,upperts,lowerts,col=1:4,xlim=c(2021,2021.4),ylab="O3 PPB")


