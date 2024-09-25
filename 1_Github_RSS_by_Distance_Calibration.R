###########################################################################################################################################################
##
##  Kristina L Paxton
##
##  April 1 2021
##
##    Code to examine the relationship between RSS values and distance for a given study area using an exponential decay model
##      -- Analysis in Ecology and Evolution Paper is based on data from a test with of a radio transmitter at 135 known locations distributed throughout the Guam Network 
##         -- At each test location a transmitter was held stationary for a 5-minute time period 
##            -- Removed first and last minute of each test to ensure that times matched between tests and the node network
##            -- For each test, calculated an average RSS value for the 3-min middle time period individually for each node that detected the transmitter
##        
##
##    1st step: Data preparation - isolates raw RSS data from node network that is associated with test data 
##      -- Creates the data file that was published with this paper at: https://doi.org/10.5066/P94LQWIE 
##            -- Data associated with this test is coded as 'A' in the column 'DataSet'
##    2nd step: Exponential Decay Model - uses the dataset created in step 1 to examine the relationship between RSS values and distance for a node network
##
##
##    Files Needed
##        1. TestInfo.csv: file with test information that is saved in the working directory defined below
##            - For each unique test location there should a row of information associated with each of the 3 inner minutes of that test
##                - meaning that each test will have 3 rows of information with each row identifying one of the 3 inner minutes of the test (see TestInfo_Example.csv)
##            - Columns
##                -- TagId: unique identifier of the transmitter used during the specified test
##                -- TestId: unique identifier given to each of the unique locations where a test was conducted
##                -- Date: date when the specified test was conducted
##                -- Start.Time: time when the specified test was started
##                -- End.Time: time when the specified test was ended
##                -- Min: 1-minute time period of the specified test
##                -- Hour: hour of the the specified minute of the test
##                -- TestUTMx: Easting location of the specified test 
##                -- TestUTMy: Northing location of the specified test 
##
##       2. BeepData.rds: file with the raw RSS values collected by the node network during the time period of the test. The file should be saved in the working directory defined below
##            -- output file from Import_beep.data.Github.R script (see BeepData_Example.rds)
##            -- If using a file created from another source the following columns are needed in the specified format: 
##                -- TagId: Factor identifying the unique code of a tag
##                -- Time.local: POSIXct value in the format: "2021-12-29 09:02:51" that identifies the unique datetime stamp (IN THE LOCAL TIME ZONE OF THE STUDY) when the tag was detected by the specified node
##                    ***** If you don't have a column with the local time, but only UTC time - then change lines 42-44 in Functions_RSS.Based.Localizations.R from 'Time.local' to 'Time' ************
##                    ***** BUT be sure that Dates and Times in TestInfo.csv are also in UTC time and not local time *********
##                -- NodeId: Factor identifying the node in the network that detected the specified tag
##                -- TagRSSI: Integer indicating the RSS value (in dB) of the signal of the specified tag received by the specified node
##
##       3. Nodes.csv: file with a list of all nodes in your network and their UTM locations. The file should be saved in the working directory defined below
##            - Columns
##                -- NodeId: unique identifier of given to each node 
#                 -- NodeUTMx: Easting location of the specified node
##                -- NodeUTMy: Northing location of the specified node 
##            - Other columns can also be in the file for your reference 
##            - If Node names are being converted to scientific notation in Excel open the file with a text editor (e.g. BBedit) to change names to the correct format and save the file  
## 
##       4. Functions_RSS.Based.Localizations.R - R script with functions needed to run the code below - file should be saved in the working directory defined below    
##      
##  Important Output
##      1. Calibration Dataset that contains the average RSS values for a given node that detected the signal of a test transmitter at a known location
##      2. Model coefficients from an exponential decay model that show the relationship between RSS and Distance for a node network
##            exponential model formula: avgRSS ~ a * exp(-S * distance) + K
##                  a = intercept
##                  S = decay factor
##                  K = horizontal asymptote
##
##  
##########################################################################################################################################################

# Packages needed
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(data.table)

# Reset R's brain - removes all previous objects
rm(list=ls())

node_file <- function(health) {
  if (nrow(health) < 1) stop("no node health data!")
  #health$timediff <- as.integer(health$time - health$recorded_at)
  #health <- health[health$timediff == 0,]
  health <- aggregate(health[,c("latitude", "longitude")],list(health$node_id), mean, na.rm=TRUE)
  if (any(is.na(health))) {health <- health[-which(is.na(health$latitude) | is.na(health$longitude)),]}
  #
  colnames(health)[colnames(health)=="latitude"] <- "node_lat"
  colnames(health)[colnames(health)=="longitude"] <- "node_lng"
  colnames(health)[colnames(health)=="Group.1"] <- "node_id"
  return(health)
}

out <- function(x, contents, timezone) {
  x <- which(names(contents) == x)
  timecol <- contents[, x]
  if (is.character(timecol)) {
    DatePattern <- "[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}[T, ][[:digit:]]{2}:[[:digit:]]{2}:[[:digit:]]{2}(.[[:digit:]]{3})?[Z]?"
    exactDatePattern <- "^[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}[T, ][[:digit:]]{2}:[[:digit:]]{2}:[[:digit:]]{2}(.[[:digit:]]{3})?[Z]?$"
    brokenrow <- grep(exactDatePattern, timecol, invert = TRUE) # find row that has a date embedded in a messed up string (i.e. interrupted rows)
    timecol[brokenrow] <- substring(timecol[brokenrow], regexpr(DatePattern, timecol[brokenrow]))
    timecol[brokenrow[which(regexpr(DatePattern, timecol[brokenrow]) < 0)]] <- NA
    newtimecol <- as.POSIXct(timecol, tz = timezone)
  } else {
    newtimecol <- timecol
  }
  return(newtimecol)
}

#data.setup(test, tag_col="Tag.Id", tagid="0C5F5CED", time_col="Time..UTC.",timezone="UTC",x="Longitude",y="Latitude", node_ids=node_ids, loc_precision=4)

tag_col <- "Tag.Id"
tagid = "0C5F5CED"
time_col="Time..UTC."
timezone="UTC"
x="Longitude"
y="Latitude"
loc_precision=4

data.setup <- function(test, tag_col, tagid, time_col, timezone, x, y, node_ids = NA, loc_precision = NA, latlon=T) {
  test$tag_id <- test[,tag_col]
  test <- test[test$tag_id %in% tagid,]
  outtime <- lapply(time_col, out, contents=test, timezone=timezone)
  test[time_col] <- outtime
  if(length(time_col) < 2) {
    test$time <- test[,time_col]
    if(is.na(loc_precision)) {
      if(any(colnames(test) == "beep_number")) {
        test$TestId <- seq(1:nrow(test))
        test.info <- test
       #paste(format(test$Time, "%Y-%m-%d %H:%M"), test$beep_number, sep="_")
        test.info$Start.Time <- test.info$time - 1
        test.info$Stop.Time <- test.info$time + 1
        test.info$lat <- test[,y]
        test.info$lon <- test[,x]
      }
    } else {
      multfactor = 10^loc_precision
      test$lat <- format(trunc(test[,y]*multfactor)/multfactor, nsmall=loc_precision)
      test$lon <- format(trunc(test[,x]*multfactor)/multfactor, nsmall=loc_precision)
      test$id <- paste(test$lat, test$lon, sep="_")
      test <- test %>%
        mutate(c_diff = ifelse(id != lag(id), 1, 0))
      test$c_diff[1] <- 0
      test$TestId <- cumsum(test$c_diff)
      test.info <- setDT(test)[, .(Start.Time = min(time), Stop.Time = max(time)), by = TestId]
      test.info$id <- test$id[match(test.info$TestId, test$TestId)]
      testid <- test[!duplicated(test$id),]
      test.info$TestId <- testid$TestId[match(test.info$id, testid$id)]
      test.info$lat <- as.numeric(test$lat[match(test.info$TestId, test$TestId)])
      test.info$lon <- as.numeric(test$lon[match(test.info$TestId, test$TestId)])
    }
  } else {
    test.info <- test
    test.info$lat <- test[,y]
    test.info$lon <- test[,x]}
  
  test.UTM <- test.info %>%
    dplyr::group_by(TestId) %>%
    dplyr::slice_head(n=1)
  
  #if(length(tagid) > 1) {
    df1 <- test %>%
      group_by(TestId) %>%
      summarise(tag_id = list(unique(tagid))) %>%
      unnest(tag_id)
    test.info <- left_join(test.info, df1) 
  #}
  
  start <- min(test.info$Start.Time)
  end <- max(test.info$Stop.Time)
  
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = "/home/jess/Documents/radio_projects/data/radio_projects/office/my-db.duckdb", read_only = TRUE)
  testdata <- tbl(con, "blu") |>
    filter(time >= start & time <= end) |>
    filter(tag_id %in% tagid) |>
    collect()
  
  #testdata$syncid <- paste(format(testdata$time, "%Y-%m-%d %H:%M"), testdata$sync, sep="_")
  
  start_buff = start - 2*60*60
  end_buff = end + 2*60*60
  
  nodes <- tbl(con, "node_health") |>
    filter(time >= start_buff  & time <= end_buff) |>
    collect()
  
  DBI::dbDisconnect(con)
  
  test.dat <- setDT(testdata)[test.info,  TestId := +(i.TestId), on = .(tag_id, time > Start.Time, time < Stop.Time), by = .EACHI]
  test.dat <- test.dat[!is.na(test.dat$TestId),]
  
  summary.test.tags <- test.dat %>%
    dplyr::group_by(node_id, TestId) %>%
    dplyr::summarise(avgRSS = mean(tag_rssi),
                     sdRSS = sd(tag_rssi),
                     n.det = n())
  
  if(any(!is.na(node_ids))) {nodes <- nodes[nodes$node_id %in% node_ids,]}
  
  nodes <- node_file(nodes)
  
  dst <- raster::pointDistance(test.UTM[,c("lon", "lat")], nodes[,c("node_lng", "node_lat")], lonlat = latlon, allpairs = T)
  dist_df <- data.frame(dst, row.names = test.UTM$TestId)
  colnames(dist_df) <- nodes$node_id
  dist_df$TestId <- as.integer(rownames(dist_df))
  
  dist.gather <- dist_df %>%
    tidyr::gather(key = "node_id", value = "distance", -TestId)
  
  summary.dist <- summary.test.tags %>%
    dplyr::left_join(dist.gather) 
  
  # Add UTMs of nodes and test locations to dataset
  combined.data <- summary.dist %>%
    dplyr::left_join(nodes[, c("node_id", "node_lat", "node_lng")]) %>%
    dplyr::left_join(test.UTM[, c("TestId", "lat", "lon")]) 
return(combined.data)}

## Set by User
# Working Directory - Provide/path/on/your/computer/where/master/csv/file/of/nodes/is/found/and/where/Functions_CTT.Network.R/is/located
working.directory <- "/home/jess/Documents/radio_projects/EcolEvol.Manuscript_Optimizing.Trilateration"

# Directory for Output data - Provide/path/where/you/want/output/data/to/be/stored/
outpath <- "/home/jess/Documents/radio_projects/paxton/"

## Bring in functions 
setwd(working.directory)
source("4_Functions_RSS.Based.Localizations.R")

#re-sample data for calibration 
#create time window by reducing location precision
#or can input data with TestId column (user-defined window)

#INPUT
#test DATA FRAME looks like...
#tag_col, time_col, x, y, latlon=T
#must pass Start.Time, Stop.Time, TestId if you pass your own
#tagid = character or vector of tag ids
#time_col = character or vector of time columns
#timezone only applies to if your time column(s) are characters

#OUTPUT: combined/matched up data set

mytest <- read.csv("~/Downloads/cal_20m_up.csv")

node_ids = c(
  "B25AC19E",
  "44F8E426",
  "FAB6E12",
  "1EE02113",
  "565AA5B9",
  "EE799439",
  "1E762CF3",
  "A837A3F4",
  "484ED33B"
)

testout <- data.setup(mytest, tag_col="Tag.Id", tagid="0C5F5CED", time_col="Time..UTC.",timezone="UTC",x="Longitude",y="Latitude", node_ids=node_ids) #, loc_precision=4


## Bring in 3 Needed files - Test Information, RSS values, and Node Information - change file names in " " as needed
#test.info_k <- read.csv("Test.Info_Example.csv", header = T)

#test <- test[-which(test$Tag.Type=="userLocation"),]



#testdata$TestId <- 0L

str(test.info) # check that data imported properly

#beep.dat <- readRDS("BeepData_Example.rds") 
#str(beep.dat) # check that data imported properl

#nodes <- read.csv("Nodes_Example.csv", header = T)
#str(nodes)



# rearrange data



## Combine distances with summary data 



#################################################################
#   Step 1
# Run function to isolate raw RSS data from a node network that
# is associated with test data  
#################################################################

    ## Variables to define for function
        ## TEST.TYPE = category indicating the type of dataset
            ## "Calibration" - test dataset used for calibrating the relationship between RSS and Distance (purpose of this R script)
            ## "LocError" - test data used to determine localization error associated with RSS-based localization estimates for a node network
        ## DATE.FORMAT = format of the date column in Test.Info.csv file using R standard date expressions (e.g., "%m-%d-%y", "%Y-%m-%d")
        ## TIMEZONE = Time zone where data was collected, use grep("<insert location name here>", OlsonNames(), value=TRUE) to find valid time zone name


    ## Output in R environment
        ## combined.data - dataframe that contains the average RSS values for a given node associated with each unique test
            ## Columns:
                ## NodeId - unique identifier of the specified node
                ## TestId - unique identifier of the specified test
                ## avgRSS - average RSS value across a 3-min time period for the given node and test
                ## sdRSS - standard deviation of RSS values across the 3-min time period for the given node and test
                ## distance - true distance (in meters) between the specified node and test location
                ## NodeUTMx - Easting location of the specified node
                ## NodeUTMy - Northing location of the specified node
                ## TestUTMx - Easting location of the specified test
                ## TestUTMy - Northing location of the specified test

    ## Output saved
        ## Calibration_Dataset.csv - .csv file of the dataframe 'combine.data' saved in the folder specified by the outpath


##* Define Variables - replace values below with user specified values **## 
#TEST.TYPE <- "Calibration"
#DATE.FORMAT <- "%m/%d/%y"
#TIME.ZONE <- "Pacific/Guam"


# Combine RSS data from a node network with test information into a dataframe
#combined.data <- data.setup(TEST.TYPE, DATE.FORMAT, TIME.ZONE)

atomic <- function(test, nodes, testdata) { #make sure these are the right data frames
dst <- raster::pointDistance(test[,c("Longitude", "Latitude")], nodes[,c("node_lng", "node_lat")], lonlat = T, allpairs = T)
dist_df <- data.frame(dst, row.names = test$syncid)
colnames(dist_df) <- nodes$node_id
dist_df$syncid <- rownames(dist_df)

# rearrange data
dist.gather <- dist_df %>%
  tidyr::gather(key = "node_id", value = "distance", -syncid)

summary.dist <- testdata %>% 
  dplyr::left_join(dist.gather) 

# Add UTMs of nodes and test locations to dataset
combined.data <- summary.dist %>%
  dplyr::left_join(nodes[, c("node_id", "node_lat", "node_lng")]) %>%
  dplyr::left_join(test[, c("syncid", "lat", "lon")]) 
}
##########################################################
#    Step 2
# Exponential Decay Function to Examine Relationship 
# between distance and Tag RSS values 
##########################################################


## Visualize data

  # Plot of the relationship between RSS and distance
ggplot(data = combined.data, aes(x = distance, y = avgRSS, color = node_id)) +
  geom_point(size = 2)



## Preliminary Exponential Decay Model to determine starting values for the final model 

    # SSasymp - self start for exponential model to finding the starting values of the data
    # Asym - horizontal asymptote (when large values ) - y values decay to this value 
    # R0 - numeric value when avgRSS (i.e., response variable) = 0 
    # lrc - natural logarithm of the rate constant (rate of decay)


  # Preliminary Model
#combined.data <- combined.data[!is.na(combined.data$distance),]
exp.mod <- nls(avgRSS ~ SSasymp(distance, Asym, R0, lrc), data = combined.data)
  # Summary of Model
summary(exp.mod)
  # rate of decay
exp(coef(exp.mod)[["lrc"]])




## Final Exponential Decay Model with user provided self-starting values 
## based on visualization of the data and values in the Preliminary Model Output

    # exponential model formula: avgRSS ~ a * exp(-S * distance) + K
      # a = intercept
      # S = decay factor
      # K = horizontal asymptote


##  ***** Variables to define for final model below - replace values below with values from exp.mod ****  ## 
a <- coef(exp.mod)[["R0"]]
S <- exp(coef(exp.mod)[["lrc"]])
K <- coef(exp.mod)[["Asym"]]
  
  # Final Model
nls.mod <- nls(avgRSS ~ a * exp(-S * distance) + K, start = list(a = a, S = S, K= K), 
               data = combined.data)
  # Model Summary
summary(nls.mod)
  # Model Coefficients 
    # **** you will use these coefficient values to estimate distance in Github_TestDataset_LocalizationError.R and Github_Simulations.R scripts *****
coef(nls.mod)


## Check the fit of the model and get predicted values

  # Get residuals and fit of model and add variables to main table
combined.data$E <- residuals(nls.mod)
combined.data$fit <- fitted(nls.mod)

  # Plot residuals by fit or distance
ggplot(combined.data, aes(x = distance, y = E, color = node_id)) +
         geom_point(size = 2)
ggplot(combined.data, aes(x = fit, y = E, color = node_id)) +
  geom_point(size = 2)

  # Get model predictions
combined.data$pred <- predict(nls.mod)


## Save Final Dataset with Residuals and Predictions
#write.csv(combined.data, paste0(outpath, "Calibration_Dataset_withResiduals_Predictions.csv"),
#          row.names = F)


## Plot with predicted line

#pdf(paste0(outpath, "Relationship_Distance~RSSI.pdf"), width = 8, height = 6)

ggplot(combined.data, aes(x = distance, y = avgRSS)) + 
  geom_point() +
  geom_line(aes(y = pred), color = "#689FBB", lwd = 1.25) +
  scale_y_continuous(name = "RSS (dB)") +
  scale_x_continuous(name = "Distance (m)") +
  theme_classic()

#dev.off()




