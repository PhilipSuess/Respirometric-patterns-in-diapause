#### extract DGE parameters ####

#### load libraries ####
library(ggplot2)
library(tidyverse)
library("dplyr")
library("zoo")
library(hrbrthemes) # themes for graphs
library(socviz) # %nin
library(stringi)

####load data table with individual ID####
ID <- "dia_0_1"
stwd <- "Chapters/Chapter III respirometry/Data"
setwd(stwd)
resptab <- read.table(paste(ID,".txt",sep = ""))

#### sd for finding the open phase ####
sdofcom <- sd(na.omit(resptab$co2))

#### reduce size of respiratory file to decrease computing time ####
resptab_op <- resptab
resptab_op <- resptab_op %>%
  slice(which(row_number() %% 10 == 1))
resptab_op <- resptab_op %>%
  dplyr::mutate(rollpeak = zoo::rollmean(co2, k = 60, fill = NA))

#### check the fit of standard deviation threshold for open phase visually ####
ggplot(subset(resptab_op))+
  geom_point(aes(timeh, rollpeak))+
  geom_hline(yintercept = 2*sdofcom)+
  geom_hline(yintercept = 2.5 * sdofcom)+
  theme_bw()

#### everything that is higher than 2* the standard deviation but smaller than 2.5* the standard deviation is thought to be
#### involved in a peak 
resptab_op$isop <- NA

for(i in 1:length(resptab_op$time)){
  if(isTRUE(resptab_op$rollpeak[i] > (2*sdofcom) & resptab_op$rollpeak[i] < (2.5 * sdofcom))==T){
    resptab_op$isop[i] <- 1
  }
  else{
    resptab_op$isop[i] <- 0
  }
}
#check what counts as open phase and compare to the actual traces 
ggplot(resptab_op)+
  geom_point(aes(time, isop))+
  theme_bw()


#### get end time of the open phases ####
timeofpeaks <- resptab_op$time[resptab_op$isop == 1]
timeofpeaks <- na.omit(timeofpeaks)

# reduce the found times for the open phase to only one time per open phase
end <- 1
for(i in 1:length(timeofpeaks)){
  if(isTRUE(timeofpeaks[i] +  3600 > timeofpeaks[i+1]) == F){
    end[i] <- timeofpeaks[i]
  }
  else{
    end[i] <- NA
  }
}
end <- na.omit(end)
allpeaks <- end[end>3600] # avoid peaks in the beginning due to handling
resptab$peak <- 0
resptab$peak[resptab$time %in% end] <- 1 # add the found times of the open phase to the original datafile
#check visually for open phase and found times
ggplot(subset(resptab))+
  geom_line(aes(time,co2))+
  geom_vline(xintercept = resptab$time[resptab$peak == 1], col = "red")+
  theme_bw()

#### get the maximum value of the open phase and timing ####
allpeaks <- c(3600,allpeaks)

max <- 1:(length(allpeaks)-1)
timepeak <- 1:(length(allpeaks)-1)

for(i in 1:(length(allpeaks)-1)){
  max[i] <- max(resptab$co2[resptab$time > allpeaks[i] & resptab$time < allpeaks[i+1]])
  timepeak[i] <- min(resptab$time[resptab$time > allpeaks[i] & resptab$time < allpeaks[i+1] & resptab$co2 == max[i]])
  
}


#### create a table with the maximum values ####
maxtable <- as.data.frame(matrix(NA,length(max),2)) # table with the length of open phases
colnames(maxtable) <- c("max","time") # name the columns
maxtable$max <- max # add the maximum values of an open phase
maxtable$time <- timepeak # add the timing of the maximum value of the open phase

mxtbl <- subset(maxtable, time > 3600) # only take open phases that are not immediately at the start or the end
mxtbl <- subset(mxtbl, time < length(resptab$time) - 3600) # only take open phases that are not immediately at the start or the end

resultstable <- as.data.frame(matrix(NA,length(mxtbl$max),11)) #create a results table with the length of the found open phases

#### name the table and fill with the max amplitude and the timing ####
colnames(resultstable) <- c("peak_nr","peak_amplitude","time_max_peak","peak_length","peak_area","period","cf_phase","cf_max","cf_area", "closed_phase","cf_spike_nr")
resultstable$peak_nr <- 1:(length(mxtbl$max)) # number the open phases
resultstable <- as.numeric(resultstable) # make everything numeric

resultstable$peak_amplitude <- mxtbl$max # maximum of open phase added to the results table
resultstable$time_max_peak <- mxtbl$time # timing of the open phase added to the results table

#### decrease file size to limit the computing time ####
resptab <- resptab %>%
  dplyr::mutate(rollpeak = zoo::rollmean(co2, k = 600, fill = NA))

#### find beginning and end of a peak ####
pknr <- seq(1,length(allpeaks),2)
end <- 1

for(i in 1:nrow(resultstable)){
  time <- resultstable$time_max_peak[i]
  CO_2 <- resultstable$peak_amplitude[i]
  while(CO_2 > 0.2*sdofcom){
    
    CO_2 <- resptab$rollpeak[resptab$time == time]
    
    time = time +1
  } 
  end[i] <- time
}


start <- 1
for(i in 1:nrow(resultstable)){
  time <- resultstable$time_max_peak[i]
  CO_2 <- resultstable$peak_amplitude[i]
  while(CO_2 > 0.2 * sdofcom){
    CO_2 <- resptab$co2[resptab$time == time]
    time = time - 1
  }
  start[i] <- time
}

ggplot(resptab)+
  geom_line(aes(timeh, co2))+
  geom_vline(xintercept = start/3600, col = "green")+
  geom_vline(xintercept = end/3600, col = "red")+
  geom_hline(yintercept = 0.2 * sdofcom, col = "red")+
  theme_bw()

#### add it to the original file ####
resptab$peak <- 0
resptab$peak[resptab$time %in% end] <- 1
resptab$peak[resptab$time %in% start] <- 1

#### check if all peaks are marked in the original file ####

allpeaks <- resptab$time[resptab$peak == 1]

ggplot(resptab)+
  geom_line(aes(time, co2))+
  geom_vline(xintercept = allpeaks, col = "red")+
  theme_bw()

#### get the amount of CO2 produced per open phase ####
pknr <- seq(1,length(allpeaks),2)
area <- 1
for(i in pknr){
  area[i] <- sum(resptab$co2[resptab$time > allpeaks[i] & resptab$time < allpeaks[i+1]])
}

area <- na.omit(area)
resultstable$peak_area <- area / 60

#### get peak length ####
peak_length <- 1
for(i in pknr){
  peak_length[i] <- allpeaks[i+1] - allpeaks[i]
}


peak_length <- na.omit(peak_length)

resultstable$peak_length <- peak_length

#### get period ####

period <- 1
for(i in pknr){
  period[i] <- allpeaks[i+3] - allpeaks[i +1]  
}

period <- na.omit(period)

resultstable$period <- append(NA,period)

#### cf phase ####

cf_phase <- 1
for(i in pknr){
  cf_phase[i] <- allpeaks[i+2] - allpeaks[i+1]  
}

cf_phase <- na.omit(cf_phase)


resultstable$cf_phase <- append(NA,cf_phase)

#### cf max ####

cf_max <- 1
for(i in pknr){
  cf_max[i] <- max(resptab$co2[resptab$time < allpeaks[i+2] & resptab$time > allpeaks[i+1]])  
}

cf_max <- na.omit(cf_max)


resultstable$cf_max <- append(NA,cf_max)

#### cf area ####
cf_area <- 1
for(i in pknr){
  cf_area[i] <- sum(resptab$co2[resptab$time < allpeaks[i+2] & resptab$time > allpeaks[i+1]])  
}

cf_area <- na.omit(cf_area)


resultstable$cf_area <- append(NA,cf_area) / 60

####closed phase ####

pknr <- seq(1,length(allpeaks),2)
startc <- 1
for(i in 1:nrow(resultstable)){
  time <- resultstable$time_max_peak[i]
  CO_2 <- resultstable$peak_amplitude[i]
  while(CO_2 > 0.015*sdofcom){
    
    CO_2 <- resptab$rollpeak[resptab$time == time]
    
    time = time +1
  } 
  startc[i] <- time
}



pknr <- seq(1,length(allpeaks),2)
endc <- 1
for(i in 1:length(startc)){
  time <- startc[i]
  CO_2 <- resptab$co2[resptab$time == time]
  while(CO_2 < 0.1*sdofcom){
    CO_2 <- resptab$co2[resptab$time == time]
    time = time +1
  } 
  endc[i] <- time
}


ggplot(resptab_op)+
  geom_vline(xintercept = endc, col = "red")+
  geom_vline(xintercept = startc, col = "green")+
  geom_line(aes(time,co2))+
  geom_hline(yintercept = 0.02, col = "magenta")+
  theme_bw()


closed_length <- endc - startc

resultstable$closed_phase <- closed_length

#### closed phase peaks ####

cf_peaks <- 1
for(i in pknr){
  cf_peaks <- append(cf_peaks,resptab$time[resptab$time > allpeaks[i+1] & resptab$co2 > 0.2 * sdofcom])
}

cf_peaks <- unique(cf_peaks)

cf_peak_times <- 1
for(i in 1:length(cf_peaks)){
  if(isTRUE(cf_peaks[i] + 300 > cf_peaks[i+1])==F){
    cf_peak_times[i] <- cf_peaks[i]
  }
  else{
    cf_peak_times[i] <- NA
  }
}
cf_peak_times <- na.omit(cf_peak_times)

cf_peak_nr <- 1
for(i in pknr){
  cf_peak_nr[i] <- length(resptab$time[resptab$time > allpeaks[i+1] + 60 & resptab$time < allpeaks[i+2] - 60 & resptab$time %in% cf_peak_times])
}

cf_peak_nr <- na.omit(cf_peak_nr)

cf_peak_nr[-length(cf_peak_nr)]

resultstable$cf_spike_nr <- append(NA,cf_peak_nr[-length(cf_peak_nr)])


resultstable$id <- ID

ggplot(resptab)+
  geom_line(aes(timeh, co2))+
  geom_vline(xintercept = start/3600, col = "green")+
  geom_vline(xintercept = end/3600, col = "red")+
  geom_hline(yintercept = 0.2 * sdofcom, col = "red")+
  theme_bw()

#### save results ####

write.table(resultstable, file = paste(ID,"results.txt",sep="_"))


