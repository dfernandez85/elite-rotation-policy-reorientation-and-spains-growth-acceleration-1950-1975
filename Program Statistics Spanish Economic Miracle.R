#Program statistics Spanish Economic Miracle
remove(list=ls())

# Load data
library(rio)
dat <- import("https://www.rug.nl/ggdc/historicaldevelopment/maddison/data/mpd2018.dta")

# Cleaning data
# Year as Date
library(lubridate)
library(zoo)
dat$year <- as.Date(as.yearmon(dat$year))

# Selecting countries
countries <- c("Austria", "Belgium", 
               "Denmark", "Finland", "France", "Germany",
               "Greece", "Italy","Netherlands", "Norway",
               "Portugal", "Spain" , "Sweden", "Switzerland",
               "United Kingdom")
               
dat <- dat[(dat$country) %in% countries, ]
dat$country <- as.character(dat$country)

# Selecting years and variables
dat <- dat[which(dat$year >= "1850-01-01"), c(2:3, 5:6)]

# Subsetting dataset only Spain
dat1 <- dat[which(dat$country == "Spain"),]

# Creating YoY growth moving average for Spain
library(dplyr)
library(pracma)
library(quantmod)
dat1$growth <- Delt(dat1$rgdpnapc, k = 1)
dat1$mm <- movavg(dat1$growth, n = 10, type = "s")

# Creating dummy variables for historical periods
library(tidyr)
library(dplyr)
dat1$regimes <- 0
dat1$regimes[which(dat1$year > 1867 & dat1$year < 1876)] <- 1
dat1$regimes[which(dat1$year > 1875 & dat1$year < 1924)] <- 2
dat1$regimes[which(dat1$year > 1923 & dat1$year < 1931)] <- 3
dat1$regimes[which(dat1$year > 1930 & dat1$year < 1936)] <- 4
dat1$regimes[which(dat1$year > 1935 & dat1$year < 1940)] <- 5
dat1$regimes[which(dat1$year > 1939 & dat1$year < 1976)] <- 6
dat1$regimes[which(dat1$year > 1978 & dat1$year <= 2016)] <- 7

dat1$regimes <- as.factor(dat1$regimes)
dat1$regimes <- recode(dat1$regimes, 
                       "1" = "Democratic Sexennium",
                       "2" = "Bourbon Restoration",
                       "3" = "Primo Rivera Dictatorship",
                       "4" = "Second Republic",
                       "5" = "Civil War",
                       "6" = "Francoism",
                       "7" = "Current Democracy"
                       )

# Economic growth different regimes
Dates <- data.frame(start = as.Date(c("1868-01-01", "1876-01-01", 
                                      "1923-01-01", "1931-01-01", 
                                      "1936-01-01", "1940-01-01",
                                      "1978-01-01")),
                    end =   as.Date(c("1875-01-01", "1922-01-01", 
                                      "1930-01-01", "1935-01-01",
                                      "1939-01-01", "1975-01-01", 
                                      "2016-01-01")),
                    regimes = factor(c("Democratic Sexennium",
                                "Bourbon Restoration",
                                "Primo Rivera Dictatorship",
                                "Second Republic",
                                "Civil War",
                                "Francoism",
                                "Current Democracy"))
                            )

Dates$regimes<-  factor(Dates$regimes, levels = Dates$regimes)         

# Creating plots
# Plot Growth per capita GDP Spain from 1868 to 2016
library(ggplot2)
library(scales)

plot1 <- ggplot() +
  geom_rect(data = Dates, alpha = 0.9,
            aes(xmin = start,
                xmax = end,
                ymin = -Inf,
                ymax = Inf,
                fill = regimes)) +
  geom_line(data = dat1, aes(year, mm), lwd = 1.2) +
  labs(title = "Figure 1 Economic Growth in Spain during Different Regimes",
       y = "Growth in real GDP per capita (10-year moving average)",
       x = "Year") +
  theme_minimal() + 
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 0.5,
                                   vjust = 0.5),
        axis.text = element_text(face = "bold")) +
  scale_y_continuous(breaks=seq(-0.03, 0.06, 0.01),
                     minor_breaks = NULL,
                     labels = percent) +
  scale_x_date(limits = as.Date(c("1868-01-01", "2016-01-01")),
               expand = c(0, 0),
               date_breaks = "5 years",
               minor_breaks = NULL,
               date_labels = "%Y") +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             colour = "black") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()
        )
                     
# Plot GDP per Capita Spain vs Rest of Western Europe
# Create income (weighted mean by population) for Western Europe (excluding Spain)

income <- dat %>%
  group_by(year) %>% 
  summarise(Western_Europe = weighted.mean(rgdpnapc[which(country != "Spain")], 
                                               pop[which(country != "Spain")], 
                                               na.rm = TRUE))

# Create data frame Spanish Income, Western Europe Income
dat2 <- data.frame(dat1$year, dat1$rgdpnapc, income$Western_Europe)
colnames(dat2) <- c("year", "Spain", "Western_Europe")
dat2 <- dat2 %>%
  mutate(Share = Spain / Western_Europe)

# Create plot

plot2 <- ggplot() +
  geom_rect(data = Dates, alpha = 0.9,
            aes(xmin = start,
                xmax = end,
                ymin = -Inf,
                ymax = Inf,
                fill = regimes)) +
  geom_line(data = dat2, aes(year, Share), lwd = 1.2) +
  labs(title = "Figure 2 GDP per Capita (Spain / Rest of Western Europe)",
       y = "Spain income as % of the rest of Western Europe",
       x = "Year") +
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 0.5,
                                   vjust = 0.5),
        axis.text = element_text(face = "bold")) +
  scale_y_continuous(breaks=seq(0.0, 1, 0.1),
                     minor_breaks = NULL,
                     labels = percent) +
  scale_x_date(limits = as.Date(c("1868-01-01", "2016-01-01")),
               expand = c(0, 0),
               date_breaks = "5 years",
               minor_breaks = NULL,
               date_labels = "%Y") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_blank()
  )

ggsave("Figure_1.png", plot1, width = 20, 
       height = 11, units = "cm", dpi = 600)
ggsave("Figure_2.png", plot2, width = 20, 
       height = 11, units = "cm", dpi = 600)
