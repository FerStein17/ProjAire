rm(list=ls())

library(dynlm)
library(zoo)
library(lmtest)
library(sandwich)
library(dplyr)

if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github('diegovalle/aire.zmvm')

if (!require("pacman")) install.packages("pacman")
pacman::p_load(aire.zmvm, dplyr, ggplot2, ggseas)

## download pollution data by station in ppb
o3 <- get_station_data(criterion = "MAXIMOS", # Can be one of MAXIMOS (daily maximum), 
                       # MINIMOS (daily minimum), 
                       # or HORARIOS (hourly average)
                       pollutant = "O3", # Can be one of "SO2", "CO", "NOX", "NO2", "NO", "O3", 
                       # "PM10", "PM25", "WSP", "WDR", "TMP", "RH"
                       year = 2017:2021) # A numeric vector, the earliest year allowed is 1986



o3_max <- o3 %>%
  group_by(date) %>%
  summarise(max = ifelse(all(is.na(value)),
                         NA,
                         base::max(value, na.rm = TRUE))) %>%
  na.omit()

# We only need the data from 11-Oct-2017
o3_max = o3_max[284:1572,]


#We now extract data for PM10
pm10 <- get_station_data(criterion = "MAXIMOS", # Can be one of MAXIMOS (daily maximum), 
                       # MINIMOS (daily minimum), 
                       # or HORARIOS (hourly average)
                       pollutant = "PM10", # Can be one of "SO2", "CO", "NOX", "NO2", "NO", "O3", 
                       # "PM10", "PM25", "WSP", "WDR", "TMP", "RH"
                       year = 2017:2021) # A numeric vector, the earliest year allowed is 1986



pm10_max <- pm10 %>%
  group_by(date) %>%
  summarise(max = ifelse(all(is.na(value)),
                         NA,
                         base::max(value, na.rm = TRUE))) %>%
  na.omit()

# We only need the data from 11-Oct-2017
pm10_max = pm10_max[284:1572,]

names_pm10 <-c("date", "pm10_max")
names_o3 <-c("date","o3_max")

colnames(pm10_max)<-names_pm10
colnames(o3_max)<-names_o3


#Now we need Traffic Data
setwd("C:/Users/ferna/OneDrive/Escritorio/Decimo Semestre/Macroeconometria Aplicada/Proyecto")
traffic <- read.csv('Traffic_Data_CDMX.csv')


traffic$date<-as.Date(traffic$date,"%m/%d/%Y")
traffic <- traffic[,1:2]

#Turn them into TS
max_o3_ts <- ts(o3_max$o3_max, start = c(2017,10), frequency = 365.25)
max_pm10_ts <- ts(pm10_max$pm10_max, start = c(2017.10), frequency = 365.25)
traffic_ts <- ts(traffic$Trafico, start = c(2017,10), frequency = 365.25)

ts.plot(max_o3_ts,traffic_ts,col=1:2,xlim=c(2017,2021),ylab="Trafico y O3")

ts.plot(max_pm10_ts,traffic_ts,col=1:2,xlim=c(2017,2021),ylab="Trafico y O3")

