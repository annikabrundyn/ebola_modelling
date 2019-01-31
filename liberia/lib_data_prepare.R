###### Prepare Liberia data ------

# clear environment and set working directory
rm(list = ls())
wd = "~/Liberia code"
setwd(wd)

# install required libraries
library(ggplot2)
library(grid)
library(gridExtra)

# read in the raw data
full_data_raw = read.csv("CDCallcasecounts.csv", sep = ";")

full = full_data_raw
colnames(full)[1] = "Date"
full$Date = as.Date(full$Date, format = "%Y/%m/%d")
full = full[-266,] #remove empty final row

# only look at dataframe columns for Liberia - column 4 and 5
lib.cases = full[,c(1,4,5)]

# remove any missing observations
lib.cases = na.omit(lib.cases)

# rename column headings
colnames(lib.cases)[2] = "Cases"
colnames(lib.cases)[3] = "Deaths"

# plot raw LIB case data
lib_raw_cases = ggplot(lib.cases, aes(Date, Cases)) +
  geom_point(shape = 1) + 
  ggtitle("Cumulative Cases") + 
  theme_minimal()

# plot raw LIB death data
lib_raw_deaths = ggplot(lib.cases, aes(Date, Deaths)) +
  geom_point(shape = 1) + 
  ggtitle("Cumulative Deaths") + 
  theme_minimal()

# arrange plots of raw data in a grid format
title1=textGrob("Liberia Raw Data", gp=gpar(fontface="bold", fontsize = 16, col = "tomato"))
grid.arrange(lib_raw_cases, lib_raw_deaths, nrow = 1, ncol = 2, 
             top = title1)

# remove duplicate dates
dup.dates = which( duplicated( lib.cases$Date))
lib.data = lib.cases[-dup.dates,]

# order data from earliest date to latest date
lib.data = lib.data[ order(lib.data$Date), ]

# start from date when cases are first > 50
row.lib = min(which(lib.data$Cases >= 50))      
start.lib = lib.data$Date[row.lib]  
lib.data = lib.data[c(row.lib:nrow(lib.data)),]

# convert dates to days after initial starting date
lib.data$Day = as.numeric(difftime(lib.data$Date, start.lib-1, units = "days"))  
lib.data = lib.data[,c(1,4,2,3)]
rownames(lib.data) = 1:nrow(lib.data)

# identify position of dips in cumulative cases
lib.data$inc = c(NA, lib.data$Cases[2:nrow(lib.data)] - lib.data$Cases[1:(nrow(lib.data)-1)] )
lib.data[c(125, 126, 128, 129, 130, 132, 133, 134, 135), 2] = 10672
neg.inc = which(lib.data$inc < 0)
lib.data = lib.data[-c(39, 118),]
rownames(lib.data) = 1:nrow(lib.data)

#plot(Deaths ~ Date, data = lib.data, type="l")
#plot(Deaths ~ Date, data = lib.data , subset = seq(1,nrow(lib.data),5), type="l")
#plot(Cases ~ Date, data = lib.data , type="l")
#plot(inc ~ Date, data = lib.data , type="l")

# smooth the death data to only select every 5th point
lib.data$Deaths[-seq(1, nrow(lib.data), 5)] = NA

# plot non missing death indices
death_ind <- ! is.na(lib.data$Deaths)
plot(Deaths ~ Date, data = lib.data , subset = death_ind, type="l")

# delete a few inconsistent observations - clear errors
lib.data = lib.data[-c(123:133),]

# plot cleaned LIB case data
lib_clean_cases = ggplot(lib.data, aes(Date, Cases)) +
  geom_point(shape = 1) + 
  ggtitle("Cumulative Cases") + 
  theme_minimal()

# plot cleaned LIB death data
lib_clean_deaths = ggplot(lib.data, aes(Date, Deaths)) +
  geom_point(shape = 1, size = 1.6) + 
  ggtitle("Cumulative Deaths") + 
  theme_minimal()

# arrange plots of cleaned data in grid format
title1=textGrob("Liberia Cleaned Data", gp=gpar(fontface="bold", fontsize = 16, col = "tomato"))
grid.arrange(lib_clean_cases, lib_clean_deaths, nrow = 1, ncol = 2, 
             top = title1)

# save liberia data in rds and csv format
saveRDS(lib.data, file = "lib_data_23Sept.rds")
write.csv(lib.data, "lib_data_23Sept.csv")
