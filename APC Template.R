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
# Adding an minimum age cutoff if desired
age_min = 0
df <- df %>% filter(Age > age_min)
# Renaming all lowercase so that the hexamap function can work properly
names(df) <- c('age','period','cohort','rate','pop')

# Now Creating Hexamap plots: (Make sure the plot window is sufficiently large in RStudio)
rmax <- max(df$rate)  # Getting the max value to set the color range, can manually set if desired
plot_APChexamap(dat = df, y_var = "rate", color_range = c(0,rmax), color_vec = hcl.colors(rmax), label_size = 0.9)

# Creating Standard line graphs:
# Incidence by age
df %>% select(age, period, rate, pop) %>% group_by(age) %>% summarise(rate = weighted.mean(rate, pop)) %>% ggplot(aes(x = age, y = rate)) + geom_line() +
  ggtitle("Incidence by Age") + xlab("Age (Years)") + ylab("Incidence per 100k") + ylim(0, rmax)

# Incidence by period
df %>% select(age, period, rate, pop) %>% group_by(period) %>% summarise(rate = weighted.mean(rate, pop)) %>% ggplot(aes(x = period, y = rate)) + geom_line() +
  ggtitle("Incidence by Period") + xlab("Period (Calendar Years)") + ylab("Incidence per 100k") + ylim(0, rmax)

# incidence by age for 5 groups of cohorts
df %>% select(age, cohort, rate) %>% mutate(cohort = cut(cohort, 5, dig.lab=10)) %>% group_by(age, cohort) %>% summarise(rate = mean(rate)) %>%
  ggplot(aes(x = age, y = rate, colour=cohort)) + geom_line(linewidth=0.8, alpha=0.8) + 
  ggtitle("Incidence by Age, stratified by Cohort") + xlab("Age (Years)") + ylab("Incidence per 100k") + ylim(0, rmax)

# Defining APC Mutual Information score function
score <- function(df){
  test1 <- df  # Saving unbinned df
  df <- discretize(df, disc = 'equalwidth', nbins = min(length(unique(df$age)), length(unique(df$period))))
  test2 <- df  # Saving binned df
  rinfo <- entropy(df$rate)  # self info
  remaining <- condentropy(df$rate, data.frame(df$age, df$period, df$cohort))  # Info after conditioning
  
  # MI for each predictor
  a_mi <- min(mutinformation(df$rate, df$age), condinformation(df$rate, df$age, df$cohort))
  p_mi <- min(mutinformation(df$rate, df$period), condinformation(df$rate, df$period, df$cohort))
  c_mi <- min(mutinformation(df$rate, df$cohort), condinformation(df$rate, df$cohort, df$age), condinformation(df$rate, df$cohort, df$period))
  mi <- c(a_mi, p_mi, c_mi)
  
  # Calculating Uncertainty
  test1$id <- 1:nrow(test1)
  test2$id <- 1:nrow(test2)
  test1 <- merge(test1, test2, by='id', all = T) %>% select(rate.x, rate.y)
  names(test1) <- c('unbinned','binned')
  test1 <- test1 %>% group_by(binned) %>% summarise(unbinned = mean(unbinned)) %>% 
    mutate(bin_diff = binned - lag(binned)) %>% mutate(rate_diff = unbinned - lag(unbinned)) %>% filter(!is.na(bin_diff))
  test1 <- mean(test1$rate_diff/test1$bin_diff)
  
  return(round(c(100*(rinfo - remaining)/rinfo, test1, 100*mi/sum(mi)), 3))
}


# Calculating and converting to Dataframe for export
info <- data.frame(c('Total Information Contained (%)', 'Uncertainty (Incidence per 100k)', 'Attribution to Age (%)', 'Attribution to Period (%)', 'Attribution to Cohort (%)'), score(df))
names(info) <- c("Variable","Score")

# For exporting to excel via clipboard:
write.table(info, file = "clipboard", sep = "\t", row.names = F, col.names = T)





















