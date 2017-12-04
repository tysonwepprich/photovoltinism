# GDD from Canadian sites for Galerucella
# using Daymet instead of PRISM
# output should be simple ID, DOY, GDD data.frame for analysis

fs <- list.files("data/daymet", full.names = TRUE)

outdf <- list()
for (f in 1:length(fs)){
  
  all_content = readLines(fs[f])
  skip_lead = all_content[-c(1:6)]
  data = read.csv(textConnection(skip_lead), header = TRUE, stringsAsFactors = FALSE)
  
  data$date <- strptime(paste(data$year, data$yday, sep = " "), format = "%Y %j")
  
  # gdd for each day
  LDT <- 10
  UDT <- 37.8
  
  data$GDD <- TriDD(data$tmax..deg.c., data$tmin..deg.c., LDT, UDT)
  data$site <- stringr::str_split_fixed(string = list.files("data/daymet")[f], 
                                        pattern = "_", n = 2)[, 1]
  data <- data[, c("year", "yday", "GDD", "site")]
  outdf[[f]] <- data
}

outdf <- bind_rows(outdf)
names(outdf) <- c("Year", "DOY", "GDD", "ID")

outdf$ID <- factor(outdf$ID)
levels(outdf$ID) <- c("Delta, BC", "Edmonton, AB", "North Bay, ON", "Winnipeg, MB")

saveRDS(outdf, "data/gdd_canada_sites.RDS")
