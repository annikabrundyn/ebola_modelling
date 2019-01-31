###### Prepare Sierra Leone data ------

# clear environment and set working directory
rm(list = ls())
wd = "~/Ebola/final"
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

# only look at dataframe columns for Sierra Leone - column 6 and 7
sl.data.raw = full[,c(1,6,7)]

sl.cases = sl.data.raw

# remove any missing observations
sl.cases = na.omit(sl.cases)

# rename column headings
colnames(sl.cases)[2] = "Cases"
colnames(sl.cases)[3] = "Deaths"

# plot raw SL case data
sl_raw_cases = ggplot(sl.cases, aes(Date, Cases)) +
  geom_point(shape = 1) + 
  ggtitle("Cumulative Cases") + 
  theme_minimal()
#sl_raw_cases

# plot raw SL death data
sl_raw_deaths = ggplot(sl.cases, aes(Date, Deaths)) +
  geom_point(shape = 1) + 
  ggtitle("Cumulative Deaths") + 
  theme_minimal()
#sl_raw_deaths

# arrange plots of raw data in a grid format
title1=textGrob("Sierra Leone Raw Data", gp=gpar(fontface="bold", fontsize = 16, col = "#00AFBB"))
grid.arrange(sl_raw_cases, sl_raw_deaths, nrow = 1, ncol = 2, 
             top = title1)

# remove duplicate dates
dup.dates = which( duplicated( sl.cases$Date))
sl.data = sl.cases[-dup.dates,]

# order data from earliest date to latest date
sl.data = sl.data[ order(sl.data$Date), ]

# start from date when cases are first > 50
row.sl = min(which(sl.data$Cases >= 50))      
start.sl = sl.data$Date[row.sl]     # SL: starting date "2014-06-02", observation 18
sl.data = sl.data[c(row.sl:nrow(sl.data)),]

# convert dates to days after initial starting date
sl.data$Day = as.numeric(difftime(sl.data$Date, start.sl-1, units = "days"))  
sl.data = sl.data[,c(1,4,2,3)]
rownames(sl.data) = 1:nrow(sl.data)

# identify position of dips in cumulative cases
sl.data$inc = c(NA, sl.data$Cases[2:nrow(sl.data)] - sl.data$Cases[1:(nrow(sl.data)-1)] )
neg.inc.sl = which(sl.data$inc < 0)
sl.data = sl.data[-c(5, 44, 45, 76, 141, 170),]
rownames(sl.data) = 1:nrow(sl.data)

# smooth the death data to only select every 14th point
sl.data$Deaths[-seq(1, nrow(sl.data), 14)] = NA

# plot non missing death indices
death_ind <- ! is.na(sl.data$Deaths)
plot(Deaths ~ Date, data = sl.data , subset = death_ind, type="l")

# plot cleaned SL case data
sl_clean_cases = ggplot(sl.data, aes(Date, Cases)) +
  geom_point(shape = 1) + 
  ggtitle("Cumulative Cases") + 
  theme_minimal()

# plot cleaned SL death data
sl_clean_deaths = ggplot(sl.data, aes(Date, Deaths)) +
  geom_point(shape = 1, size = 1.7) + 
  ggtitle("Cumulative Deaths") + 
  theme_minimal()

# arrange plots of cleaned data in grid format
title1=textGrob("Sierra Leone Cleaned Data", gp=gpar(fontface="bold", fontsize = 16, col = "#00AFBB"))
grid.arrange(sl_clean_cases, sl_clean_deaths, nrow = 1, ncol = 2, 
             top = title1)

# save sierra leone data in rds and csv format
saveRDS(sl.data, file = "sl_data_23Sept.rds")
write.csv(sl.data, "sl_data_23Sept.csv")
