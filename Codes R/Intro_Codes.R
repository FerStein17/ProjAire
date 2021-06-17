rm(list=ls())

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
                       year = 2012:2021) # A numeric vector, the earliest year allowed is 1986



o3_max <- o3 %>%
  group_by(date) %>%
  summarise(max = ifelse(all(is.na(value)),
                         NA,
                         base::max(value, na.rm = TRUE))) %>%
  na.omit()

#Change to IMECA to be consistent

o3_max$max <- convert_to_imeca(o3_max$max, "O3")

# ozone threshold for declaring a 'smog alert' and the dates during which they
# were valid
# source: http://www.aire.cdmx.gob.mx/descargas/ultima-hora/calidad-aire/pcaa/pcaa-modificaciones.pdf
contingencia_levels <- data.frame(imeca = c(
                                          180, 180, 150, 150),
                                  start = c(2011.5795,
                                            2012.6052,  2016.291, 2016.4986),
                                  end = c(
                                          2012.6025,    2016.2883, 2016.4959, Inf))
max_daily_df <- tsdf(ts(o3_max$max, start = c(2012,1), frequency = 365.25))

contingencia <- o3_max
contingencia$date <- max_daily_df$x
contingencia$contingencia <- case_when(
  contingencia$date > 2012.6052 & contingencia$max > 180 ~ TRUE,
  contingencia$date > 2016.291 & contingencia$max > 150 ~ TRUE,
  TRUE~ FALSE
)

ggplot(max_daily_df,
       aes(x = x, y = y)) +
  geom_line(colour = "grey75", alpha = .5) +
  stat_rollapplyr(width = 30, align = "right", color = "#01C5D2") +
  geom_segment(data = contingencia_levels,
               aes(x=start, y=imeca, xend=end, yend=imeca), color="darkred",
               linetype = 2)  +
  geom_point(data=filter(contingencia, contingencia == TRUE),
             aes(x=date, y=max), color = "#999999",
             size = 2, shape = 21, fill = "firebrick3" ) +
  xlab("date") +
  ylab("IMECA") +
  scale_x_continuous(breaks = c(2012, 2014,2016, 2018,2020)) +
  ggtitle("Concentracion máxima diaria de ozono",
          subtitle = paste0("Líneas Rojas: Valores necesarios para activar Fase I ",
                            "\nLíneas Azules: Promedio de 30 días lines are ",
                            "Puntos Rojos: Indican cuando el ozono rebasa el límite" ,
                            "\nElaboración Propia con Datos de SEDEMA")) +
  theme_bw()


##########
#  PM10  #
##########

pm10 <- get_station_data(criterion = "MAXIMOS", # Can be one of MAXIMOS (daily maximum), 
                       # MINIMOS (daily minimum), 
                       # or HORARIOS (hourly average)
                       pollutant = "PM10", # Can be one of "SO2", "CO", "NOX", "NO2", "NO", "O3", 
                       # "PM10", "PM25", "WSP", "WDR", "TMP", "RH"
                       year = 2012:2021) # A numeric vector, the earliest year allowed is 1986



pm10_max <- pm10 %>%
  group_by(date) %>%
  summarise(max = ifelse(all(is.na(value)),
                         NA,
                         base::max(value, na.rm = TRUE))) %>%
  na.omit()

#Change to IMECA to be consistent

pm10_max$max <- convert_to_imeca(pm10_max$max, "PM10")

contingencia_levels <- data.frame(imeca = c(
  176, 176, 150, 150),
  start = c(2011.5795,
            2012.6052,  2016.291, 2016.4986),
  end = c(
    2012.6025,    2016.2883, 2016.4959, Inf))
max_daily_df <- tsdf(ts(pm10_max$max, start = c(2012,1), frequency = 365.25))

contingencia <- pm10_max
contingencia$date <- max_daily_df$x
contingencia$contingencia <- case_when(
  contingencia$date > 2012.6052 & contingencia$max > 176 ~ TRUE,
  contingencia$date > 2016.291 & contingencia$max > 150 ~ TRUE,
  TRUE~ FALSE
)

ggplot(max_daily_df,
       aes(x = x, y = y)) +
  geom_line(colour = "grey75", alpha = .5) +
  stat_rollapplyr(width = 30, align = "right", color = "#01C5D2") +
  geom_segment(data = contingencia_levels,
               aes(x=start, y=imeca, xend=end, yend=imeca), color="darkred",
               linetype = 2)  +
  geom_point(data=filter(contingencia, contingencia == TRUE),
             aes(x=date, y=max), color = "#999999",
             size = 2, shape = 21, fill = "firebrick3" ) +
  xlab("date") +
  ylab("IMECA") +
  scale_x_continuous(breaks = c(2012, 2014,2016, 2018,2020)) +
  ggtitle("Concentracion máxima diaria de PM10",
          subtitle = paste0("Líneas Rojas: Valores necesarios para activar Fase I ",
                            "\nLíneas Azules: Promedio de 30 días lines are ",
                            "Puntos Rojos: Indican cuando el ozono rebasa el límite" ,
                            "\nElaboración Propia con Datos de SEDEMA")) +
  theme_bw()


