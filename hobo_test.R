# Combine multiple raw HOBO data files into a single .csv file 
# Author: Bryan Urban
# Email: burban@fraunhofer.org
# Date: 2015-04-30
# https://gist.github.com/bjurban/ef74a5accf42c43480a1

## SETUP --------------------

# install these two packages first if you don't have them:
# install.packages("data.table")
# install.packages("lubridate")
library(data.table)
library(lubridate)

# change these to match the folders containing the data
raw_dir <- "C:/Users/Tyson/Documents/HOBOware/Calibration"
out_dir <- "C:/Users/Tyson/Documents/HOBOware"


## LOAD RAW DATA ------------
# get file names: 
pattern = ".*csv$" # for identifying files to read
fns <- list.files(raw_dir, pattern=pattern, full.names=TRUE)

# load data into lists
read_and_label <- function(x,...){
  z <- fread(x,...)
  
  # add file name without the extension as id column
  pattern <- "(.*\\/)([^.]+)(\\.csv$)"
  z$ids <- sub(pattern, "\\2", x)
  z
}

# reads columns 2 and 3 (timestamp and temperature) into a list of data.table
all_data <- 
  lapply(fns, function(x,...) {try(read_and_label(x,...))},
         select=2:4, header=FALSE, skip=2
  ) 

## PROCESS RAW DATA ---------
# drop errors, merge into one large data.table, name columns, parse timestamp
all_data <- all_data[sapply(all_data, is.data.table)]
all_data <- rbindlist(all_data)

setnames(all_data, c("ts", "temp", "lux", "ids"))
all_data[, ts:=floor_date(mdy_hms(ts), "minute")] # floor ts to nearest minute

# drop rows with missing data
all_data <- all_data[complete.cases(all_data),]

## SAVE PROCESSED DATA ------
# write data in one big file
write.csv(all_data, paste(out_dir, "all_data.csv", sep="/"), 
          row.names=FALSE)

library(dplyr)
all_data <- all_data %>% 
  group_by(ts) %>% 
  mutate(meantemp = mean(temp),
         meanlux = mean(lux)) %>% 
  ungroup() %>% 
  mutate(vartemp = temp - meantemp,
         varlux = lux - meanlux)

# look at mean/var over time to see sensor bias/gc mean differences
library(ggplot2)
plt <- ggplot(all_data, aes(x = ts, y = vartemp, group = ids, color = ids)) +
  geom_line() +
  theme_bw()
plt

meandat <- all_data %>% 
  filter(ts < as.Date("2018-06-15"), ids == "GC1") %>% 
  tidyr::gather(var, val, meantemp:meanlux)

vardat <-  all_data %>% 
  filter(ts < as.Date("2018-06-15")) %>% 
  tidyr::gather(var, val, vartemp:varlux)

plt <- ggplot(meandat, aes(x = ts, y = val, group = var)) +
  geom_line() +
  facet_wrap(~var, ncol = 1, scales = "free") +
  theme_bw()
plt

plt <- ggplot(vardat, aes(x = ts, y = val, group = ids, color = ids)) +
  geom_line() +
  facet_wrap(~var, ncol = 1, scales = "free") +
  theme_bw()
plt

vardat %>% 
  group_by(ids, var) %>% 
  summarise(meanval = mean(val))
