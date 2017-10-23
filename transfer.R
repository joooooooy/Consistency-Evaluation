doc <- "jp.1.2017-07-17"

file.address <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/Simu July (lamba.ctr pwr hr)/lambda/Intern project/",
                   doc,
                   ".csv",
                   sep = "")
  

test <- read.csv(file=file.address,
                 header = FALSE,
                 skip=1)
test <- test[,2:2001]
HR.class <- test

test <- read.csv(file=file.address,
                 header = FALSE,
                 skip=1)
test <- test[,2:2001]
byes.min <- test

test <- read.csv(file=file.address,
                 header = FALSE,
                 skip=1)
test <- test[,2:2001]
byes.mean <- test

test <- read.csv(file=file.address,
                 header = FALSE,
                 skip=1)
test <- test[,2:2001]
jp <- test


last <- 103
mydate <- strptime("07/24/2017:11:53:00",format = "%m/%d/%Y:%H:%M:%S")
duration <- 3600*last
mydate+duration

