# APC Simulated Data Generation
# Contact: Matthew Prest, mp4090@columbia.edu
rm(list = ls()) # Clearing Environment and Plots
gc() # Freeing unused memory
#Importing Packages,  moved all import statements here for clarity.
library(ggplot2) # For efficient plotting
library(dplyr)  # For dataframe manipulation
library(APCtools)  # For Hexamap
library(infotheo)  # For Mutual Information
library(stringr)  # For string manipulation
cat("\014") # Clearing Console

# Importing SEER Data, starting with Incidence by Year & Age
# df <- read.csv(file.choose())  # Use this line to navigate to file path if needed
df <- read.csv("W:\\Matt P\\Projects\\202305 APC\\S8 GC Incidence Age Year.csv")

# Cleaning SEER Data, getting year and age as integer values only
df <- df %>% filter(!(Year %in% c("1975-2019", "Statistic could not be calculated.", "Rates are per 100,000 and age-adjusted to the 2000 US Std Population (single ages to 84 - Census P25-1130) standard.")))
df <- df %>% filter(!(Age %in% c("85+ years", "Unknown", "", "  ~")))
df$Age <- as.integer(str_remove(df$Age, " years"))
df$Year <- as.integer(df$Year)
# Removing commas from the rate, count and pop variables, and coercing into floats
df$Rate <- as.numeric(gsub(",", "", df$Rate))
df$Count <- as.numeric(gsub(",", "", df$Count))
df$Pop <- as.numeric(gsub(",", "", df$Pop))
# Adding birth year
df$Cohort <- df$Year - df$Age
# Only keeping, age, period, cohort and rate variables
df <- df %>% select(Age, Year, Cohort, Rate, Pop)
# Renaming all lowercase so that the hexamap function can work properly
names(df) <- c('age','period','cohort','rate','pop')
df$pop <- NULL

p_min <- min(df$period)
c_min <- min(df$cohort)

# Variables to test:
form <- c('Linear', 'Quadratic', 'Sinusoidal')
predictors <- c('A','P','C','AP','AC','PC','APC')
noise_lvl <- c(0.01, 0.1, 0.5)
vars <- expand.grid(form, predictors, noise_lvl)

df$rate <- 0.0
df <- merge(df, vars, all = T)
names(df) <- c("age", "period", "cohort", "rate", "form", "predictor", "noise_lvl")
df <- df %>% select(age, period, cohort, form, predictor, noise_lvl, rate)
rm(vars, form, noise_lvl, predictors)
df$predictor <- as.character(df$predictor)


# Building simulated data
df$rate <- ifelse(df$form == 'Linear' & str_detect(df$predictor, 'A'), df$rate + df$age, df$rate)
df$rate <- ifelse(df$form == 'Linear' & str_detect(df$predictor, 'P'), df$rate + df$period-p_min, df$rate)
df$rate <- ifelse(df$form == 'Linear' & str_detect(df$predictor, 'C'), df$rate + df$cohort-c_min, df$rate)

df$rate <- ifelse(df$form == 'Quadratic' & str_detect(df$predictor, 'A'), df$rate + df$age^2, df$rate)
df$rate <- ifelse(df$form == 'Quadratic' & str_detect(df$predictor, 'P'), df$rate + (df$period-p_min)^2, df$rate)
df$rate <- ifelse(df$form == 'Quadratic' & str_detect(df$predictor, 'C'), df$rate + (df$cohort-c_min)^2, df$rate)

df$rate <- ifelse(df$form == 'Sinusoidal' & str_detect(df$predictor, 'A'), df$rate + 5+sin(df$age), df$rate)
df$rate <- ifelse(df$form == 'Sinusoidal' & str_detect(df$predictor, 'P'), df$rate + 5+sin(df$period-p_min), df$rate)
df$rate <- ifelse(df$form == 'Sinusoidal' & str_detect(df$predictor, 'C'), df$rate + 5+sin(df$cohort-c_min), df$rate)

# Getting group means to set the noise scales too
grpd <- df %>% group_by(form, predictor) %>% summarise(g_mean = mean(rate))
df <- merge(df, grpd, by=c('form','predictor'), all = T)
df$rate <- df$rate + df$noise_lvl*df$g_mean*runif(length(df$rate))
df$g_mean <- NULL
rm(grpd, c_min, p_min)
df <- df %>% select(form, predictor, noise_lvl, age, period, cohort, rate)

write.csv(df, "W:\\Matt P\\Projects\\202305 APC\\Simulated_APC_Data.csv", row.names = F)


















# GAM without covariates
# mod1 <- gam(rate ~ te(age, cohort, bs = "ps", k = c(8,8)), data = df)
# 
# plot_partialAPCeffects(mod1, dat = df, variable = 'age')







