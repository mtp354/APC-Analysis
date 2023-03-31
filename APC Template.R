# APC Template Code
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
df <- read.csv("W:\\Matt P\\Projects\\202201 CISNET\\Uterine\\APC\\S8 GC Incidence Age Year.csv")

# Cleaning SEER Data, getting year and age as integer values only
df <- df %>% filter(!(Year %in% c("1975-2019", "Statistic could not be calculated.", "Rates are per 100,000 and age-adjusted to the 2000 US Std Population (single ages to 84 - Census P25-1130) standard.")))
df <- df %>% filter(!(Age %in% c("85+ years", "Unknown", "", "  ~")))
df$Age <- as.integer(str_remove(df$Age, " years"))
df$Year <- as.integer(df$Year)
df$Rate <- as.numeric(df$Rate)
# Adding birth year
df$Cohort <- df$Year - df$Age
# Only keeping, age, period, cohort and rate variables
df <- df %>% select(Age, Year, Cohort, Rate)
names(df) <- c('age','period','cohort','rate')

# Now Creating Hexamap plots:
rmax <- max(df$rate)  # Getting the max value to set the color range, can manually set if desired
plot_APChexamap(dat = df, y_var = "rate", color_range = c(0,rmax), color_vec = hcl.colors(rmax), label_size = 0.9)

# Now Calculating Mutual Information Scores for comparative effect strength
# 1st step, discretize data
df <- as.matrix(df)  # Convert to matrix
df <- discretize(df)  # Discretizing

# Mutual information matrix:
info <- round(mutinformation(df), 7)

# Converting to % contribution
info <- info[4,1:3]/info[4,4]
info <- 100*info/sum(info)

# Converting to Dataframe for export
info <- data.frame(info)
names(info) <- c("Variable")

# For exporting to excel via clipboard:
write.table(info, file = "clipboard", sep = "\t", row.names = T, col.names = T)






