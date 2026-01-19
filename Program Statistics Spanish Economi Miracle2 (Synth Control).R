#Program synthetic control Spanish Economic Miracle
remove(list=ls())

# Load data
library(pwt9)
dat1 <- force(pwt9.1)

# Cleaning data
countries <- c("Argentina", "Australia", "Austria", "Belgium", 
               "Bolivia (Plurinational State of)",
               "Brazil", "Canada", "Colombia", "Costa Rica", "Denmark",
               "Ecuador", "El Salvador", "Finland", "France", "Germany",
               "Guatemala", "Honduras", "Ireland", "Israel", "Italy",
               "Mexico", "Netherlands", "New Zealand", "Nicaragua", "Norway",
               "Panama", "Peru", "Spain", "Sweden", "Switzerland",
               "Trinidad and Tobago", "Turkey", "United Kingdom",
               "United States of America", "Uruguay", 
               "Venezuela (Bolivarian Republic of)")

library(textclean)
dat <- dat1[(dat1$country) %in% countries, ]
dat$country <- as.character(dat$country)
dat$country <- mgsub(dat$country, 
                     "Bolivia (Plurinational State of)", 
                     "Bolivia")
dat$country <- mgsub(dat$country, 
                     "Venezuela (Bolivarian Republic of)", 
                             "Venezuela")

# Create variable GDP per Capita PPP (using rgdpo)
dat$gdpcap <- dat$rgdpo / dat$pop

# Create variable GDP per Capita (using rgdpna) (use this or the other)
dat$gdpcap <- dat$rgdpna / dat$pop

# Create variable capital per capita (using rkna)
dat$rknacapita <- dat$rnna / dat$pop

# Eliminate variables and observations not used in analysis
dat <- dat[which(dat$year < 1976) , c(1:3, 7, 10, 40:45, 53:54)]

# Create average of all countries included in sample
library(dplyr)
avg <- dat %>%
  group_by(year) %>%
  summarise_at(vars(-c(1:2)), funs(mean(.)))
avg$country <- rep("Average", 26)
avg$isocode <- rep("AVG", 26)
dat <- rbind(dat, avg)
dat <- dat[c(937:962, 1:936), ]

# Create variable Region (numerical)
library(dplyr)
dat$region <- case_when(dat$isocode == "AVG" ~ 1,
                        dat$isocode == "ARG" ~ 2,
                        dat$isocode == "AUS" ~ 3,
                        dat$isocode == "AUT" ~ 4,
                        dat$isocode == "BEL" ~ 5,
                        dat$isocode == "BOL" ~ 6,
                        dat$isocode == "BRA" ~ 7,
                        dat$isocode == "CAN" ~ 8,
                        dat$isocode == "CHE" ~ 9,
                        dat$isocode == "COL" ~ 10,
                        dat$isocode == "CRI" ~ 11,
                        dat$isocode == "DEU" ~ 12,
                        dat$isocode == "DNK" ~ 13,
                        dat$isocode == "ECU" ~ 14,
                        dat$isocode == "ESP" ~ 15,
                        dat$isocode == "FIN" ~ 16,
                        dat$isocode == "FRA" ~ 17,
                        dat$isocode == "GBR" ~ 18,
                        dat$isocode == "GTM" ~ 19,
                        dat$isocode == "HND" ~ 20,
                        dat$isocode == "IRL" ~ 21,
                        dat$isocode == "ISR" ~ 22,
                        dat$isocode == "ITA" ~ 23,
                        dat$isocode == "MEX" ~ 24,
                        dat$isocode == "NIC" ~ 25,
                        dat$isocode == "NLD" ~ 26,
                        dat$isocode == "NOR" ~ 27,
                        dat$isocode == "NZL" ~ 28,
                        dat$isocode == "PAN" ~ 29,
                        dat$isocode == "PER" ~ 30,
                        dat$isocode == "SLV" ~ 31,
                        dat$isocode == "SWE" ~ 32,
                        dat$isocode == "TTO" ~ 33,
                        dat$isocode == "TUR" ~ 34,
                        dat$isocode == "URY" ~ 35,
                        dat$isocode == "USA" ~ 36,
                        dat$isocode == "VEN" ~ 37)

#Getting matrices ready
library(Synth)
library(MSCMT)
library(ggplot2)

datprep <- listFromLong(dat, unit.variable="region", 
                       time.variable="year",unit.names.variable="country")

treatment.identifier <- "Spain"
controls.identifier  <- setdiff(colnames(datprep[[1]]), 
                                c(treatment.identifier, "Average"))

times.dep  <- cbind("gdpcap" = c(1950,1959))
times.pred <- cbind("pop" = c(1950,1959), 
                    "hc" = c(1950,1959),
                    "gdpcap" = c(1950,1959),
                    "csh_c" = c(1950,1959),
                    "csh_i" = c(1950,1959),
                    "csh_g" = c(1950,1959),
                    "csh_x" = c(1950,1959),
                    "csh_m" = c(1950,1959),
                    "rknacapita" = c(1950,1959))
agg.fns <- rep("mean", ncol(times.pred))

res <- mscmt(datprep, treatment.identifier, controls.identifier, 
             times.dep, times.pred, agg.fns, seed = 1, 
             outer.optim = "DEoptim",verbose = FALSE)

# Create table characteristics Spain vs Synthetic Spain
library(reactable)
dattable <- res$predictor.table
colnames(dattable) <- c("Spain", "Synthetic Spain", "Sample")
reactable(dattable)

# Create Plots Spain vs Synthetic Spain
library(scales)
plot3 <- ggplot(res, type="comparison", ylab = "Real GDP per capita", xlab = "Year", 
       main = "Figure 3 GDP per Capita Spain and Synthetic Spain", 
       labels = c("Spain","Synthetic Spain")) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 0.5,
                                   vjust = 0.5),
        axis.text = element_text(face = "bold")) +
  scale_x_date(limits = as.Date(c("1950-01-01", "1975-01-01")),
               expand = c(0, 0),
               date_breaks = "5 years",
               minor_breaks = NULL,
               date_labels = "%Y") +
  scale_y_continuous(breaks=seq(2500, 17500, 2500),
                     minor_breaks = NULL,
                     labels = dollar)

plot4 <- ggplot(res, type="gaps", ylab = "Differential GDP per capita", xlab = "Year", 
                main = "Figure 4 GDP per Capita Spain over Synthetic Spain", 
                labels = c("Spain","Synthetic Spain")) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 0.5,
                                   vjust = 0.5),
        axis.text = element_text(face = "bold")) +
  scale_x_date(limits = as.Date(c("1950-01-01", "1975-01-01")),
               expand = c(0, 0),
               date_breaks = "5 years",
               minor_breaks = NULL,
               date_labels = "%Y") +
  scale_y_continuous(breaks=seq(0, 5000, 1000),
                     minor_breaks = NULL,
                     labels = dollar)

ggsave("Figure_3.png", plot3, width = 20, 
       height = 11, units = "cm", dpi = 600)

ggsave("Figure_4.png", plot4, width = 20, 
       height = 11, units = "cm", dpi = 600)

# Placebo test
library(parallel)
cl <- makeCluster(detectCores())
resplacebo <- mscmt(datprep, treatment.identifier, controls.identifier,
                    times.dep, times.pred, agg.fns, cl=cl, placebo=TRUE,
                    seed=1, outer.optim="DEoptim")
stopCluster(cl)

plot5 <- ggplot(resplacebo, exclude.ratio=5, ratio.type="mspe",
                ylab = "Differential GDP per capita", xlab = "Year", 
                main = "Figure 5 Placebo Test") +
          theme(axis.ticks = element_blank(),
                axis.text.x = element_text(angle = 90, 
                                           hjust = 0.5,
                                           vjust = 0.5),
                axis.text = element_text(face = "bold")) +
          scale_x_date(limits = as.Date(c("1950-01-01", "1975-01-01")),
               expand = c(0, 0),
               date_breaks = "5 years",
               minor_breaks = NULL,
               date_labels = "%Y") +
          scale_y_continuous(breaks=seq(-5000, 5000, 1000),
                minor_breaks = NULL,
                labels = dollar)

ggsave("Figure_5.png", plot5, width = 20, 
       height = 11, units = "cm", dpi = 600)

# Test p-value
pvalue(resplacebo, exclude.ratio=5, ratio.type="mspe", alternative="less")
ggplot(resplacebo, exclude.ratio=5, ratio.type="mspe", 
       type="p.value",alternative="less", size = 3)
did(resplacebo, range.post=c(1970,1990), exclude.ratio=5, alternative="less")

# Post-pre mspe ratio
library(tibble)
ppratio <- as.data.frame(ppratio(resplacebo, type = "mspe"))
colnames(ppratio) <- c("ppratio")
ppratio <- rownames_to_column(ppratio, var = "Country")

# Plot Post-pre mspe ratio
plot6 <- ggplot() + 
  geom_col(data = ppratio, aes(ppratio, y = reorder (Country, ppratio)), 
             fill = "black") +
  geom_col(data = ppratio[1,], aes(x = ppratio, y = reorder (Country, ppratio)), 
             fill = "red") +
  labs(x = NULL, 
       y = "Country", 
       title = "Figure 6 post-pre MSPE ratio") + 
  scale_x_continuous(breaks=seq(0, 600, 100),
                     minor_breaks = NULL)

ggsave("Figure_6.png", plot6, width = 20, 
       height = 11, units = "cm", dpi = 600)
