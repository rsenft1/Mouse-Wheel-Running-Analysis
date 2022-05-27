#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# Step 1: Setup ----------------------------------------------
######### USER: make a new R project directory and put this script inside. Also add the "scripts" folder into the project directory

## Load necessary libraries
#This package-checking chunk of code is from vikram: https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
packages = c("tidyr", "here",
             "ggplot2", "grid", "stringr", "matrixStats", "reshape2", "logr", "dplyr")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

## Make directories
dir_scripts <- here("scripts") #where the scripts are
dir_rawData <- here("data", "01_Raw_data") #where the organized data is
dir_orgData <- here("data", "02_Organized_data") #where the organized data is
dir_data_clean <- here("data", "03_Clean_data") #where the cleaned data is
dir_data_indSpeed <- here(dir_data_clean, "01_IndSpeed")
dir_data_AvgSpeedPerDay <- here(dir_data_clean, "02_AvgSpeedPerDay")
dir_data_PercentRunningPerDay <- here(dir_data_clean, "03_PercentRunningPerDay")
figDir <- here("figs")
resultsDir <- here("results")
docDir <- here("doc")
for(d in list(figDir, resultsDir, docDir, dir_data_clean, dir_rawData, dir_orgData, dir_data_AvgSpeedPerDay, dir_data_PercentRunningPerDay, dir_data_indSpeed)){dir.create(d)}
options(stringsAsFactors=FALSE)
#open the log: 
logpath <- here(docDir, "log.txt")
lf <- log_open(logpath)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# Step 2: Organize, clean, and load data ---------------------
######### USER: add .log files into the data > Raw_data folder made by Step 1. 

## organize .log files into 1 folder per animal
source(here(dir_scripts, "01_organize_data_files.R")) 

## load in functions necessary for loading and cleaning data:
source(here(dir_scripts, "02_load_and_clean_data_functions.R")) 
source(here(dir_scripts, "02_load_and_clean_data.R")) 

#list of animals to be tested. Must also be the name of the folder storing the files for each mouse.
micelist <- list.dirs(here(dir_orgData),full.names=FALSE) 
micelist <- micelist[2:length(micelist)] #list.dirs will give the wd as the first element in the list. This gets rid of that. 
log_print(paste(length(micelist), "mice found for experiment"))
myTime <- getTime(here(dir_orgData,micelist[1]))
log_print(paste("Time for experiment:", myTime))
genotypes <- getGenotypes(dir_orgData)
loadAndCleanData(dir_orgData, dir_data_indSpeed, micelist, myTime, toFilter=TRUE, dayofChange=14, newlightsOnTime=7, newlightsOffTime= 12, oldlightsOnTime = 7, oldlightsOffTime=19)

#Open all the files you created as dataframes in r:
for(mouse in micelist){
  name <- paste0(mouse, "dt")
  temp <- read.csv(here(dir_data_indSpeed, paste0(name, ".csv")), header=TRUE, sep=",")
  assign(name, temp)
}

# Get light cycle information
df = get(paste0(micelist[1],"dt"))
lightCycle <- getlightCycle(get(paste0(micelist[1],"dt")))

# This is needed for graphing day/night cycles: (df can be any data frame from the experiment)
cycleMappedOn <- df$time[lightCycle$start.points]
cycleMappedOff <- df$time[lightCycle$end.points]
cycleMapped <- data.frame(start = cycleMappedOn, end = cycleMappedOff) #NOTE could add this into the function if we only want to use this index with time in days

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# Step 3: Summarize data by genotype ---------------------
# summarize individual data
df <- get(paste0(micelist[1],"dt"))
hourly_avg_speed <- df %>%
  group_by(expDay, hour) %>% 
  summarise(avgspeed=mean(speed), time=sum(time), hour=mean(hour), activeFactor=round(mean(activeFactor))) %>%  
  na.omit()

hourly_avg_speed$row_number <- as.integer(row.names(hourly_avg_speed))  

#smooth
ggplot(data=hourly_avg_speed, aes(x=hour, y=avgspeed))+
  geom_point()+
  geom_smooth(method="loess")

#spline
spline.df <- as.data.frame(spline(hourly_avg_speed$hour, hourly_avg_speed$avgspeed, method='periodic'))

ggplot(data=hourly_avg_speed, aes(x=hour, y=avgspeed))+
  geom_point()+
  geom_line(data = spline.df, aes(x = x, y = y))


ggplot(data=hourly_avg_speed, aes(x=as.factor(hour), y=avgspeed))+
  geom_point()+
  geom_boxplot()

ggplot(data=hourly_avg_speed, aes(x=row_number, y=avgspeed))+
  geom_area(fill="cornflowerblue")+
  geom_line()+
  theme_bw()

#Scaling function for scaling an animal's speed to be between 0 and 1 for each animal
scaleCol <- function(x){(x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))}
rowSems <- function(mat){rowSds(mat, na.rm=TRUE)/sqrt(rowSums(!is.na(mat)))} 
source(here(dir_scripts, "03_summarize_by_genotype.R"))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# Step 4: Individual activity graph ---------------------
source(here(dir_scripts, "g_smooth_ind.R"))
for (i in 1:length(micelist)){
  percentDone <- (i-1)/length(micelist)*100
  print(paste0("Graphing data for Mouse ",i, " ... ", round(percentDone), "% complete"))
  smoothIndActivity(i, micelist, genotypes, cycleMapped, downSample = 15)
}

log_close()
