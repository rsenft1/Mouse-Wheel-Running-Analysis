
# Introduction ------------------------------------------------------------

# The purpose of this script is to process .log files output by the raspberry Pi
# Running wheel system...Input so far is a directory containing folders (1 per mouse)
# each containing all log files for that mouse. 
# Script is based off a script received from Jin Hyung Cho(Jesse Gray's lab) written by Jesse Gray and Jin
# March 2, 2021: I modified it significantly to make certain plots with functions and add analysis capabilities. 


# Step 1: Setup (user input required) ---------------------------------------------
#
#Would you like to filter your data?
toFilter = FALSE

#This package-checking chunk of code is from vikram: https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## First specify the packages of interest
packages = c("tidyr", "here",
             "ggplot2", "grid", "stringr", "matrixStats", "reshape2")

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
 options(stringsAsFactors=FALSE)
 library(here)
 #Use "here" package to set wd
 wd = here::i_am("Edit-R_runningwheel_RS.R")
 list.dirs(wd)
 DataFolder = "Organized_Data"
 micelist <- list.dirs(here(DataFolder),full.names=FALSE) #list of animals to be tested. Must also be the name of the 
 # folder storing the files for each mouse.
 micelist <- micelist[2:length(micelist)] #list.dirs will give the wd as the first element in the
 # list. This gets rid of that. 
 print(paste(as.character(length(micelist)), "mice found"))
 #dir.create(paste(wd,"/analysisOutput/",sep=""))
 #savePath <- paste(wd,"/analysisOutput/",sep="")
 dir.create(here("Figs"))


# Step 2: Load Functions for data ---------------------------------------------------------------
 
# 1. smoothing function for plotting the data
 smo.dist <- function( v, t ) {
  if( !t ) return(v)
  if( !(t%%2) ) t <- t+1                  # t must be odd, t=2*u+1
  u <- (t-1)/2
  n <- 3*t                                # { left tail, plateau, right tail }
  L <- length(v)
  v <- c(rep(v[1],t+u),v,rep(v[L],t+u))   # has length L+2*(t+u) = L+n-1
  d <- c(1:t,rep(t+1,t),t:1)              # this is dist. [norm. via /(2*t*(t+1))]
  w <- rep(0,L)
  for ( i in 1:n ) w <- w + d[i] * v[i:(i+L-1)]
  return( w / (2*t*(t+1)) )
 }

 # 2. function for loading the data
 getData <- function(directory,binSize){ # binSize in seconds
   
   setwd(directory)
   filelist = list.files(pattern = ".*.log")
   
   #assuming tab separated values with a header    
   z = unlist(lapply(filelist, function(x)read.table(x, header=F,skip=7))) #get rid of 6 lines
   #z = unlist(lapply(filelist, function(x)read.table(x, header=T))) 
   # conversion factor for meter/min
   z <- (z* 2*pi*0.0625 *60) # not bin size dependent 	
   
   # create time bins time expressed in days
   myTime <- (1:length(z)) / 3600 * binSize / 24 # 10 sec bins, time in days
   
   data.frame(time=myTime, speed=z, timefactor = factor(ceiling(myTime))) #basic df
 }
 addTimeInfo <- function(df, binSize=10, lightsOnTime=7, lightsOffTime= 19){
    #now add in time information
    #Current time
    currentHour <- as.numeric(substring(myTime,1,2))
    currentMinute <- as.numeric(substring(myTime,4,5))
    nRows <- length(df$time)
    #initialize variables that will be hours and whether it's the animal's active cycle
    Hours <- replicate(nRows,0)
    activeFactor <- replicate(nRows,0)
    minute <- 60/binSize #num of entries for each minute
    byMin <- seq(from=1,to=nRows, by=60/binSize)
    for (k in byMin){
       Hours[k:(k+minute-1)] <- currentHour
       currentMinute <- currentMinute+1
       currentMinute <- currentMinute%%60
       if (currentMinute==0){
          currentHour <- currentHour+1
       }
       currentHour <- currentHour%%24
    }
    #add hours to df
    Hours <- head(Hours, nRows) #fix extra entries added at the end by going in minutes
    df$hour <- Hours
    for (h in 1:length(Hours)){
       if ((lightsOnTime <= Hours[h])&(Hours[h]< lightsOffTime)){
          activeFactor[h] <- FALSE
       }
       else{
          activeFactor[h] <- TRUE
       }
    }
    #add active factor to df
    df$activeFactor <- activeFactor
    return (df)
 }
 addTimeChange <- function(df, dayofChange, binSize=10, newlightsOnTime=7, newlightsOffTime= 12, oldlightsOnTime = 7, oldlightsOffTime=19){
   #Assumes the time change was made during the LIGHT cycle
   #now add in time information for the cycle change. default to short photoperiod
   Hours <- df$hour
   startHour <- df$hour[1]
   activeFactor <- df$activeFactor
   day <- df$timefactor
   activeStart<- FALSE
   for (h in 1:length(Hours)){
     if (as.numeric(day[h])==dayofChange){
       if (Hours[h]==newlightsOffTime){ #if lights are on
         activeStart <- TRUE
       }
       if ((Hours[h]>=newlightsOffTime)&(Hours[h]<oldlightsOffTime)){ #if lights are off
         activeFactor[h] <- activeStart       
       }
     }
     else if (as.numeric(day[h])>dayofChange){
       if ((newlightsOnTime <= Hours[h])&(Hours[h]< newlightsOffTime)){
         activeFactor[h] <- FALSE
       }
       else{
         activeFactor[h] <- TRUE
       }
     }
   }
   #update active factor in df
   df$activeFactor <- activeFactor
   return (df)
 }
# 3. function to find start time
 getTime <- function(directory){
   setwd(directory)
   filelist = list.files(pattern = ".*.log")
   myTime <- scan(filelist[1], '', skip = 2, nlines = 1, sep = '\n')
   return (myTime)
 }
 #4. function to get boundaries of on/off light cycle
 #MUST ADD SOMETHING TO BE ABLE TO MODIFY LIGHT CYCLE PARTWAY THROUGH AN EXPERIMENT
 getlightCycle <- function(df){
   #use activeFactor to found boundaries for graphing output later:
   activeInfo <- rle(df$activeFactor)
   activeLengths <- unlist(activeInfo[1])
   whichCycle <- unlist(activeInfo[2])
   cycleBorders <- cumsum(activeLengths)
   activeStart <- cycleBorders[whichCycle==0]+1 # 1 after the inactive cycle ends
   activeEnd <- cycleBorders[whichCycle==1]
   #plot.data <- data.frame(start.points=lightsOff, end.points=lightsOn)
   if(whichCycle[1]==1){ #if we start in the active cycle
     activeStart <- append(0,activeStart) #first start should be beginning of the recording
   }
   if(whichCycle[length(whichCycle)]==0){ #if we end in the light cycle
     activeStart <- head(activeStart,-1) #remove the last 'start' which will be the last time point
   }
   plot.data <- data.frame(start.points=activeStart, end.points=activeEnd)
 }
 
#4. function to get animals' genotypes
 getGenotypes <- function(directory){
    setwd(directory)
    mice = list.files()
    genotypes=replicate(length(mice),"N/A")
    for (m in 1:length(mice)){
       setwd(here(directory, mice[m]))
       filelist = list.files(pattern = ".*.log")
       geno <- scan(filelist[1], '', skip = 5, nlines = 1, sep = '\n')
       if (length(geno)>0){
          genotypes[m] <- geno
       }
    }
    return (genotypes)
 }


# Step 3. Load data into data frames --------------------------------------
 dayofChange=14
 genotypes <- getGenotypes(here(DataFolder))
 #genotypes[3:4] <- 'Pet1-cre(Tg) Vglut3(fl/fl)'
 #genotypes[5] <- "Vglut3(fl/fl)"

 #get Time from first file and Genotypes 
 direct <- here(DataFolder,micelist[1])
 myTime <- getTime(direct)
 
 # Loading the data (RS: don't need to chop off because I tested all wheels before)
 ptm <- proc.time()
 for (i in 1:(length(micelist))){
   #direct <- paste(wd,"/",micelist[i],"/", sep="")
   direct <- here(DataFolder,micelist[i])
   name <- paste(micelist[i],"dt", sep="")
   print(paste0("Loading data for Mouse ",i))
   temp <- getData(direct,binSize=10)
   temp$smoothSpeed <- smo.dist(temp$speed, 200)
   temp <- addTimeInfo(temp)
   #temp <- addTimeChange(temp, dayofChange, binSize=10, newlightsOnTime=7, newlightsOffTime= 12, oldlightsOnTime = 7, oldlightsOffTime=19)
   #head(temp <- temp[-c(1:40),]) #getting rid of few lines of the data because I was testing the wheel. Chops off a couple of minutes at the beginning of the file to account for spinning the wheel.
   assign(name,temp)
 }
 ptm - proc.time()

  filter.data <- function(name, thres=20){
   df <- get(paste0(name,"dt"))
   temp <- df
   thres= thres*6 #convert meters to "speed"
   for (d in 1:max(as.numeric(temp$timefactor))){
     if (sum(temp$speed[temp$timefactor==d])<thres){
       temp$speed[temp$timefactor==d]<- NA
       temp$smoothSpeed[temp$timefactor==d]<- NA
     }
   }
   return(temp)
  }
 if (toFilter==TRUE){
    for (mouse in micelist){
       temp <- filter.data(mouse)
       assign(paste0(mouse,"dt"),temp)
    } 
 }
 
# Step 4. Light cycle information --------------------------------------
 
#get borders of light cycle from one df (same for all)
 df = get(paste0(micelist[1],"dt"))
 lightCycle <- getlightCycle(get(paste0(micelist[1],"dt")))
# This is needed for graphing day/night cycles: (df can be any data frame from the experiment)
 cycleMappedOn <- df$time[lightCycle$start.points]
 cycleMappedOff <- df$time[lightCycle$end.points]
 cycleMapped <- data.frame(start = cycleMappedOn, end = cycleMappedOff) #NOTE could add this into the function if we only want to use this index with time in days
 

# Step 5. Load Graphing Functions----------------------------------------------------------------
 
 ## GRAPH 1: Average speed for duration of experiment for individual mice
 
 smoothActivity <- function(mouseNum, micelist, genotypes, cycleMapped, downSample=5, subset=FALSE){
   #name for mouse and genotype
   myName <- micelist[mouseNum]
   myData <- get(paste0(myName,"dt"))
   if (subset!=FALSE){
     myData <- myData[myData$timefactor%in%subset, ]
     myData[with(myData, timefactor %in% subset), ] # can have multiple conditions after %in%
   }
   myGeno <- genotypes[mouseNum]
   myGenofileName <- str_replace_all(myGeno,"/","-")
   grob <- grobTree(textGrob(myName, x=0.05,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13, fontface="bold")))
   grobGeno <- grobTree(textGrob(myGeno, x=0.05,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=13, fontface="italic")))
   # the plot itself:
   p <- ggplot(data=cycleMapped, aes(xmin=start, xmax=end, ymin=0, ymax=Inf))+
     geom_rect(fill='gray', alpha=0.6)+
     geom_area(data=myData[seq(1, nrow(myData), downSample), ], inherit.aes = FALSE, aes(x=time, y=smoothSpeed), fill="cornflowerblue", alpha=0.5, color="black")+
     #geom_line(data=myData, inherit.aes = FALSE, aes(x=time, y=smoothSpeed, group=1),size=0.5)+
     scale_y_continuous(limits = c(0, 25), expand = c(0, 0))+
     scale_x_continuous(limits=c(min(myData$time), max(myData$time)), expand = c(0, 0), breaks=seq(1,max(as.numeric(myData$timefactor)),1))+
     xlab("Days")+
     ylab("Mean wheel speed (meter/min)")+
     annotation_custom(grob)+
     annotation_custom(grobGeno)+
     theme_classic()
   ggsave(
      filename = paste0(myName,"_",myGenofileName,".pdf"),
      plot = p,
      device = "pdf",
      path = here("Figs"),
      scale = 1,
      width = 11,
      height = 4,
      units = "in",
      dpi = 300,
      limitsize = TRUE
   )
   #ggsave(path = here("figs"), filename = "fig1.png")
   #pdf(paste0(name,".pdf"))
   #print(p)
   #dev.off()
 }

 
 #Scaling function for scaling an animal's speed to be between 0 and 1 for each animal

 scaleCol <- function(x){(x-min(x))/(max(x)-min(x))}
 
 
 #The following code generates a data frame object for each unique genotype
 # Run before running either of the group-based graphing functions
 rowNum <- length(df$speed)
 allGeno = unique(genotypes)
 for (geno in allGeno){
   myMice <- micelist[genotypes==geno]
   speedMatrix <- matrix(data=NA, nrow=rowNum, ncol=length(myMice))
   colnames(speedMatrix) <- myMice
   Each.averageSpeedSmooth <- speedMatrix
   mouseNum <- 0
   for (mouse in myMice){
     df <- get(paste0(mouse,"dt"))
     speedMatrix[,mouse] <- df$speed
     Each.averageSpeedSmooth[,mouse]<-df$smoothSpeed
     mouseNum <- mouseNum +1
   }
   #smooth and scale speed for each mouse in matrix
   #Each.averageSpeedSmooth <- apply(speedMatrix, 2, smo.dist, 200)
   Each.averageSpeedScaledSmooth <- apply(Each.averageSpeedSmooth, 2, scaleCol)
   #derive average speed for raw data and for smoothed
   avgSpeedRaw <- rowMeans(speedMatrix, na.rm=TRUE)
   avgSpeedSmooth <- rowMeans(Each.averageSpeedSmooth, na.rm=TRUE)
   avgSpeedSmoothScaled <- rowMeans(Each.averageSpeedScaledSmooth, na.rm=TRUE)
   assign(paste0(geno,"_speeds"), speedMatrix)
   assign(paste0(geno,"_speedsScaled"), avgSpeedSmoothScaled)
   
   #get sd and from sd, the sem of the smoothed and scaled data
   sdSpeed <- rowSds(Each.averageSpeedSmooth, na.rm=TRUE)
   semSpeed <- sdSpeed/sqrt(length(myMice))
   
   sdSpeedScaled <- rowSds(Each.averageSpeedScaledSmooth, na.rm=TRUE)
   semSpeedScaled <- sdSpeedScaled/sqrt(length(myMice))
   
   temp <- data.frame(time = df$time, timefactor = df$timefactor,avg.raw.speed = avgSpeedRaw, avg.smoothed.speed = avgSpeedSmooth, avg.scaled.smoothed.speed = avgSpeedSmoothScaled, sem.smooth = semSpeed, sem.scaled.smooth = semSpeedScaled)
   assign(geno, temp)
 }
 ## GRAPH 2: Average speed for duration of experiment for mice grouped by genotype 
#Function for graphing average of raw speed for two different genotypes
 smoothActivityGroups <- function(geno1, geno2, cycleMapped){
    #name for mouse and genotype
    dataOne = get(geno1)
    dataTwo = get(geno2)
    grobOne <- grobTree(textGrob(geno1, x=0.05,  y=0.95, hjust=0, gp=gpar(col="deepskyblue", fontsize=13, fontface="bold")))
    grobTwo <- grobTree(textGrob(geno2, x=0.05,  y=0.9, hjust=0, gp=gpar(col="tomato", fontsize=13, fontface="bold")))
    #writing this for two groups only at the moment
    fileName <- str_replace_all(paste0("Avg_speed",geno1,"_vs_",geno2),"/","-")
    # the plot itself:
    p <- ggplot(data=cycleMapped, aes(xmin=start, xmax=end, ymin=0, ymax=Inf))+
       #geom_area(data=dataTwo, inherit.aes = FALSE, aes(x=time, y=semSpeed), fill="cornflowerblue", alpha=0.5)+
       geom_rect(fill='gray', alpha=0.6)+
       geom_ribbon(data = dataOne, inherit.aes = FALSE,aes(y = avg.smoothed.speed, x = time, ymin = (avg.smoothed.speed - sem.smooth), ymax = (avg.smoothed.speed + sem.smooth)), fill = "deepskyblue", alpha=0.5)+
       geom_ribbon(data = dataTwo, inherit.aes = FALSE,aes(y = avg.smoothed.speed, x = time, ymin = (avg.smoothed.speed - sem.smooth), ymax = (avg.smoothed.speed + sem.smooth)), fill = "tomato", alpha=0.4)+
       geom_line(data = dataOne, inherit.aes = FALSE, aes(x=time, y=avg.smoothed.speed, group=1),size=0.5, color="dodgerblue4")+
       geom_line(data = dataTwo, inherit.aes = FALSE, aes(x=time, y=avg.smoothed.speed, group=1),size=0.5, linetype = "31",color="orange4")+
       ylim(0,25)+
       xlab("Days")+
       ylab("Mean wheel speed (meter/min)")+
       annotation_custom(grobOne)+
       annotation_custom(grobTwo)+
       theme_classic()
    ggsave(
       filename = paste0(fileName,".pdf"),
       plot = p,
       device = "pdf",
       path = here("Figs"),
       scale = 1,
       width = 11,
       height = 4,
       units = "in",
       dpi = 600,
       limitsize = TRUE
    )
 }
 #Function for graphing average of scaled (0 to 1) speed for two different genotypes
 smoothScaledActivityGroups <- function(geno1, geno2, cycleMapped){
    #name for mouse and genotype
    dataOne = get(geno1)
    dataTwo = get(geno2)
    grobOne <- grobTree(textGrob(geno1, x=0.05,  y=0.95, hjust=0, gp=gpar(col="deepskyblue", fontsize=13, fontface="bold")))
    grobTwo <- grobTree(textGrob(geno2, x=0.05,  y=0.9, hjust=0, gp=gpar(col="tomato", fontsize=13, fontface="bold")))
    #writing this for two groups only at the moment
    fileName <- str_replace_all(paste0("Avg_speed_all_points_",geno1,"_vs_",geno2),"/","-")
    # the plot itself:
    p <- ggplot(data=cycleMapped, aes(xmin=start, xmax=end, ymin=0, ymax=Inf))+
       geom_rect(fill='gray', alpha=0.6)+
       theme_classic()+
       geom_ribbon(data = dataOne, inherit.aes = FALSE,aes(y = avg.scaled.smoothed.speed, x = time, ymin = (avg.scaled.smoothed.speed - sem.scaled.smooth), ymax = (avg.scaled.smoothed.speed + sem.scaled.smooth)), fill = "deepskyblue", alpha=0.5)+
       geom_ribbon(data = dataTwo, inherit.aes = FALSE,aes(y = avg.scaled.smoothed.speed, x = time, ymin = (avg.scaled.smoothed.speed - sem.scaled.smooth), ymax = (avg.scaled.smoothed.speed + sem.scaled.smooth)), fill = "tomato", alpha=0.4)+
       geom_line(data = dataOne, inherit.aes = FALSE, aes(x=time, y=avg.scaled.smoothed.speed, group=1),size=0.5, color="dodgerblue4")+
       geom_line(data = dataTwo, inherit.aes = FALSE, aes(x=time, y=avg.scaled.smoothed.speed, group=1),size=0.5, linetype = "31",color="orange4")+
       xlab("Days")+
       annotation_custom(grobOne)+
       annotation_custom(grobTwo)+
       scale_y_continuous("Mean wheel speed (meter/min)", seq(0,1.1, by=0.2))
       
    ggsave(
       filename = paste0(fileName,"_scaled.pdf"),
       plot = p,
       device = "pdf",
       path = here("Figs"),
       scale = 1,
       width = ,
       height = 4,
       units = "in",
       dpi = 600,
       limitsize = TRUE
    )
 }

 ## GRAPH 3: Average speed per day by genotype (relies on the Graph 2 for grouping code chunk to be run)
 avgSpeedPerDayGroup <- function(geno1, geno2){ #note, the genotypes must be strings that match the df name for each genotype.
    #graph average as well as individuals
    df1 <- get(geno1)
    df2 <- get(geno2)
    mat1 <- get(paste0(geno1, "_speeds")) # speed matrix
    mat2 <- get(paste0(geno2, "_speeds"))
    fileName <- str_replace_all(paste0("Avg_speed_per_day",geno1,"_vs_",geno2),"/","-")
    timefactor <- df1$timefactor
    meandf1 <- setNames(aggregate(mat1,by=list(timefactor), mean), c("Day", colnames(mat1)))
    meandf2 <- setNames(aggregate(mat2,by=list(timefactor), mean), c("Day", colnames(mat2)))
    meandf1 <- melt(meandf1 ,  id.vars = 'Day', variable.name = 'Mouse', value.name = 'Average.Speed')
    meandf2 <- melt(meandf2 ,  id.vars = 'Day', variable.name = 'Mouse', value.name = 'Average.Speed')
    f1 <- function(x) c(Mean = mean(x), SEM = sd(x)/sqrt(length(x)))
    group1ALL <- do.call(data.frame, aggregate(Average.Speed~Day, meandf1, f1))
    group2ALL <- do.call(data.frame, aggregate(Average.Speed~Day, meandf2, f1))
    p <- ggplot(group1ALL, aes(Day, Average.Speed.Mean)) + 
       geom_line(data=meandf1, aes(Day, Average.Speed, color = geno1, group=Mouse), size=1, alpha=0.1)+
       geom_line(data=meandf2, aes(Day, Average.Speed, color = geno2, group=Mouse), size=1, alpha=0.1)+
       geom_line(aes(color = geno1), group=1, size=1)+
       geom_errorbar(aes(ymin=Average.Speed.Mean-Average.Speed.SEM, ymax=Average.Speed.Mean+Average.Speed.SEM), width=.2,color="turquoise4",
                     position=position_dodge(.9))+
       geom_line(data=group2ALL, aes(color = geno2), group=2, size=1)+
       geom_errorbar(data=group2ALL, aes(ymin=Average.Speed.Mean-Average.Speed.SEM, ymax=Average.Speed.Mean+Average.Speed.SEM), width=.2, color="tomato4",
                     position=position_dodge(.9))+
       
       theme_classic()+
       theme(legend.position=c(0.2,0.95),
             legend.background = element_rect(fill=alpha('white', 0)),
             legend.title = element_blank())
    ggsave(
       filename = paste0(fileName,".pdf"),
       plot = p,
       device = "pdf",
       path = here("Figs"),
       scale = 1,
       width = 6,
       height = 4,
       units = "in",
       dpi = 600,
       limitsize = TRUE
    )
 }
 
 ## GRAPH 4: % Time running by genotype (relies on the Graph 2 for grouping code chunk to be run)
 PercentTimeRunGroup <- function(geno1, geno2){
    df1 <- get(geno1)
    df2 <- get(geno2)
    mat1 <- get(paste0(geno1, "_speeds")) # speed matrix
    mat2 <- get(paste0(geno2, "_speeds"))
    fileName <- str_replace_all(paste0("Percent_time_running_per_day",geno1,"_vs_",geno2),"/","-")
    timefactor <- df1$timefactor
    meandf1 <- setNames(aggregate(mat1,by=list(timefactor), function(x){sum(x>0 )/8640*100}), c("Day", colnames(mat1)))
    meandf2 <- setNames(aggregate(mat2,by=list(timefactor), function(x){sum(x>0 )/8640*100}), c("Day", colnames(mat2)))
    meandf1 <- melt(meandf1 ,  id.vars = 'Day', variable.name = 'Mouse', value.name = 'PercentRunning')
    meandf2 <- melt(meandf2 ,  id.vars = 'Day', variable.name = 'Mouse', value.name = 'PercentRunning')
    f1 <- function(x) c(Mean = mean(x), SEM = sd(x)/sqrt(length(x)))
    group1ALL <- do.call(data.frame, aggregate(PercentRunning~Day, meandf1, f1))
    group2ALL <- do.call(data.frame, aggregate(PercentRunning~Day, meandf2, f1))
    p <- ggplot(group1ALL, aes(Day, PercentRunning.Mean)) + 
       geom_line(data=meandf1, aes(Day, PercentRunning, color = geno1, group=Mouse), size=1, alpha=0.1)+
       geom_line(data=meandf2, aes(Day, PercentRunning, color = geno2, group=Mouse), size=1, alpha=0.1)+
       geom_line(aes(color = geno1), group=1, size=1)+
       geom_errorbar(aes(ymin=PercentRunning.Mean-PercentRunning.SEM, ymax=PercentRunning.Mean+PercentRunning.SEM), width=.2,color="turquoise4",
                     position=position_dodge(.9))+
       geom_line(data=group2ALL, aes(color = geno2), group=2, size=1)+
       geom_errorbar(data=group2ALL, aes(ymin=PercentRunning.Mean-PercentRunning.SEM, ymax=PercentRunning.Mean+PercentRunning.SEM), width=.2, color="tomato4",
                     position=position_dodge(.9))+
       theme_classic()+
       theme(legend.position=c(0.2,0.95),
             legend.background = element_rect(fill=alpha('white', 0)),
             legend.title = element_blank())
    ggsave(
       filename = paste0(fileName,".pdf"),
       plot = p,
       device = "pdf",
       path = here("Figs"),
       scale = 1,
       width = 6,
       height = 4,
       units = "in",
       dpi = 600,
       limitsize = TRUE
    )
 }
 
   #GRAPH 5: Bouts of running distribution
 BoutLengthGroup <- function(geno1, geno2){
   #avg and median graphs for bout length comparison across 2 genotypes
    maxBoutLength <- 30 #min. consider bouts of activity between 0 and 30 min in length. if >30 min expected for a single bout (unlikely), then need to change this number
    fileName <- str_replace_all(paste0("BoutsRunning_",geno1,"_vs_",geno2),"/","-")
    bins <- seq(0,maxBoutLength, by=1) 
    df1 <- get(geno1)
    df2 <- get(geno2)
    mat1 <- get(paste0(geno1, "_speeds")) # speed matrix (must be raw speeds, not smoothed)
    mat2 <- get(paste0(geno2, "_speeds"))
    densityMatrix1 <- matrix(data=NA, nrow=maxBoutLength, ncol=ncol(mat1))
    densityMatrix2 <- matrix(data=NA, nrow=maxBoutLength, ncol=ncol(mat2))
    avgBoutLength1 <- replicate(ncol(mat1), NA)
    avgBoutLength2 <- replicate(ncol(mat2), NA)
    nCol <- max(ncol(mat1), ncol(mat2))
    medianBoutLength1 <- replicate(nCol, NA)
    medianBoutLength2 <- replicate(nCol, NA)
    for (c in 1:ncol(mat1)){
       r <- rle((mat1>0)[,c]) #rle for first column (first animal) of genotype speed matrix where speed is >0 (so a value of TRUE means the animal was moving, FALSE = inactive)
       len <- unlist(r[1])[unlist(r[2])==TRUE]*10/60 #this takes the lengths of active bouts (when r[2] is true, the animal's speed is >0) and converts 10 second observations to minutes
       temp <- hist(len, breaks=bins)
       avgBoutLength1[c] <- mean(len)
       medianBoutLength1[c] <- median(len)
       densityMatrix1[,c] <- temp$density
    }
    for (c in 1:ncol(mat2)){
       r <- rle((mat2>0)[,c]) #rle for first column (first animal) of genotype speed matrix where speed is >0 (so a value of TRUE means the animal was moving, FALSE = inactive)
       len <- unlist(r[1])[unlist(r[2])==TRUE]*10/60 #this takes the lengths of active bouts (when r[2] is true, the animal's speed is >0) and converts 10 second observations to minutes
       temp <- hist(len, breaks=bins)
       avgBoutLength2[c] <- mean(len)
       medianBoutLength2[c] <- median(len)
       densityMatrix2[,c] <- temp$density
    }
    length(avgBoutLength1) <- length(avgBoutLength2)
    length(avgBoutLength2) <- length(avgBoutLength1)
    avgBoutLength <- data.frame(avgBoutLength1,avgBoutLength2)
    avgBoutLengthSEM <- apply(avgBoutLength, 2, sd)/c(nrow(avgBoutLength[1]), nrow(avgBoutLength[2]))
    names(avgBoutLengthSEM) <- c(geno1, geno2)
    names(avgBoutLength) <- c(geno1, geno2)
    avgBoutLength <- gather(avgBoutLength, genotype, avgBoutLength)
    #Plot avg bout length length by genotype
    p <- ggplot(avgBoutLength, aes(x = factor(genotype), y= avgBoutLength, fill=genotype)) + 
      stat_summary(fun = mean, geom = "bar") + 
      geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
      stat_summary(fun.data = mean_se, geom = "errorbar",width=.2,
                   position=position_dodge(.9)) +
      labs(x = "Genotype")+
      labs(y="Average Bout Length (min)")+
      theme_classic()+
      theme(legend.position = "none")
    
    ggsave(
      filename = paste0(fileName,"_average.pdf"),
      plot = p,
      device = "pdf",
      path = here("Figs"),
      scale = 1,
      width = 3,
      height = 5,
      units = "in",
      dpi = 600,
      limitsize = TRUE
    )
    medBoutLength <- data.frame(medianBoutLength1,medianBoutLength2)
    medBoutLengthSEM <- apply(medBoutLength, 2, sd)/c(nrow(medBoutLength[1]), nrow(medBoutLength[2]))
    names(medBoutLengthSEM) <- c(geno1, geno2)
    names(medBoutLength) <- c(geno1, geno2)
    medBoutLength <- gather(medBoutLength, genotype, medBoutLength)
    #Plot median bout length length by genotype
    p <- ggplot(medBoutLength, aes(x = factor(genotype), y= medBoutLength, fill=genotype)) + 
      stat_summary(fun = mean, geom = "bar") + 
      geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
      stat_summary(fun.data = mean_se, geom = "errorbar",width=.2,
                   position=position_dodge(.9)) +
      labs(x = "Genotype")+
      labs(y="Median Bout Length (min)")+
      theme_classic()+
      theme(legend.position = "none")
    
    ggsave(
      filename = paste0(fileName,"_median.pdf"),
      plot = p,
      device = "pdf",
      path = here("Figs"),
      scale = 1,
      width = 3,
      height = 5,
      units = "in",
      dpi = 600,
      limitsize = TRUE
    )
    #HISTOGRAM TOO!
    avgDensity1 <- rowMeans(densityMatrix1)
    semDensity1 <- rowSds(densityMatrix1)/sqrt(ncol(mat1))
    avgDensity2 <- rowMeans(densityMatrix2)
    semDensity2 <- rowSds(densityMatrix2)/sqrt(ncol(mat2))
    #if subtract 0.5 from bins, puts the two bars between their boundaries
    newDF <- data.frame(AvgFractionTime=c(avgDensity1, avgDensity2), SEM=c(semDensity1, semDensity2), binLabel=c(bins[-1]-0.5, bins[-1]-0.5), Genotype=rep(c(geno1, geno2),each=length(avgDensity1)))
    dfAll <- data.frame(avgDensity1, avgDensity2, binLabel = bins[-1])
    dfmeansem <- data.frame(avgDensity1, avgDensity2,semDensity1,semDensity2, binLabel = bins[-1])
    dfAll <- melt(dfAll, id.vars='binLabel')
    dfAll <- melt(dfmeansem, id.vars='binLabel')
    #WORKSSSS below
    p <- ggplot(newDF, aes(x=binLabel, y=AvgFractionTime, fill=Genotype)) +
    geom_bar(position=position_dodge(), stat="identity", colour='black') +
    geom_errorbar(aes(ymin=AvgFractionTime-SEM, ymax=AvgFractionTime+SEM), width=.2,position=position_dodge(.9))+
    labs(x = "Bout length Bin Center(min)")+
    labs(y="Fraction Time Spent")+
    scale_y_continuous(breaks=(seq(0,0.8,by=0.1)), limits=NULL)+
    scale_x_continuous(breaks=(seq(0,30,by=1.5)), limits=NULL)+
    theme(legend.position=c(0.85,0.9),
          legend.background = element_rect(fill=alpha('white', 1)))
    ggsave(
      filename = paste0(fileName,"_histogram_avg.pdf"),
      plot = p,
      device = "pdf",
      path = here("Figs"),
      scale = 1,
      width = 8,
      height = 4,
      units = "in",
      dpi = 600,
      limitsize = TRUE
    )
 }
 #GRAPH 6: Average day wrapped around clock
 #GRAPH 7: Coherence with Lightcycle (% Active time IN active cycle)
 PercentActiveCoherenceGroup <- function(geno1, geno2){
   df1 <- get(geno1)
   mouse1=get(paste0(micelist[1],"dt"))
   df2 <- get(geno2)
   mat1 <- get(paste0(geno1, "_speeds")) # speed matrix
   mat2 <- get(paste0(geno2, "_speeds"))
   fileName <- str_replace_all(paste0("Percent_activity_IN_activeCycle",geno1,"_vs_",geno2),"/","-")
   timefactor <- df1$timefactor #day
   activefactor <- mouse1$activeFactor #MUST ADD THIS IN to the geno DF
   #8640 is sum(timefactor==1) or the number of time points to a day
   #to take into account changing day length, we need a way to feed in the number of points available...or write a different function?
   meandf1 <- setNames(aggregate(mat1,by=list(activefactor, timefactor), function(x){sum(x>0 )/4320*100}), c("ActiveCycle","Day", colnames(mat1)))
   meandf2 <- setNames(aggregate(mat2,by=list(activefactor), function(x){sum(x>0 )/8640*100}), c("ActiveCycle", colnames(mat2)))
   meandf1 <- melt(meandf1 ,  id.vars = c('ActiveCycle', "Day"), variable.name = 'Mouse', value.name = 'PercentRunning')
   meandf2 <- melt(meandf2 ,  id.vars = 'ActiveCycle', variable.name = 'Mouse', value.name = 'PercentRunning')
   f1 <- function(x) c(Mean = mean(x), SEM = sd(x)/sqrt(length(x)))
   group1ALL <- do.call(data.frame, aggregate(PercentRunning~ActiveCycle, meandf1, f1))
   group2ALL <- do.call(data.frame, aggregate(PercentRunning~ActiveCycle, meandf2, f1))
   p <- ggplot(group1ALL, aes(ActiveCycle, PercentRunning.Mean)) + 
     geom_line(data=meandf1, aes(ActiveCycle, PercentRunning, color = geno1, group=Mouse), size=1, alpha=0.1)+
     geom_line(data=meandf2, aes(ActiveCycle, PercentRunning, color = geno2, group=Mouse), size=1, alpha=0.1)+
     geom_line(aes(color = geno1), group=1, size=1)+
     geom_errorbar(aes(ymin=PercentRunning.Mean-PercentRunning.SEM, ymax=PercentRunning.Mean+PercentRunning.SEM), width=.2,color="tomato4",
                   position=position_dodge(.9))+
     geom_line(data=group2ALL, aes(color = geno2), group=2, size=1)+
     geom_errorbar(data=group2ALL, aes(ymin=PercentRunning.Mean-PercentRunning.SEM, ymax=PercentRunning.Mean+PercentRunning.SEM), width=.2, color="turquoise4",
                   position=position_dodge(.9))+
     theme_classic()+
     theme(legend.position=c(0.2,0.95),
           legend.background = element_rect(fill=alpha('white', 0)),
           legend.title = element_blank())
   ggsave(
     filename = paste0(fileName,".pdf"),
     plot = p,
     device = "pdf",
     path = here("Figs"),
     scale = 1,
     width = 6,
     height = 4,
     units = "in",
     dpi = 600,
     limitsize = TRUE
   )
 }
 #GRAPH 8: Actogram
 actogram <- function(mouse_name, micelist, genotypes, cycleMappedOn, cycleMappedOff, binSize = 5){
   #takes mouse *name* for data frame object, bins data, and makes bar plot/actogram of that data across days
   #default bin size is 10 min
   fileName <- paste0(mouse_name,"_actogram")
   #1st # is 288 for 5 and 144 for 10, second number is 30 for 5 and 60 for 10
   mouse_df=get(paste0(mouse_name, "dt"))
   binsCol <- rep(1:288, each=30) #factor representing 5 min bins for 1 day
   epochBins = rep(binsCol, length.out = length(mouse_df$timefactor))
   Mode <- function(vector) {
     choices <- unique(vector)
     choices[which.max(tabulate(match(vector, choices)))]
   }
   data <- setNames(aggregate(mouse_df$speed, by=list(epochBins, mouse_df$timefactor), sum),c("Bin", "Day", "Sum Speed/5min")) #divided by 6 epochs per min #CHECK THIS, UNSURE
   hourInfo <- setNames(aggregate(as.numeric(mouse_df$hour), by=list(epochBins, mouse_df$timefactor), Mode),c("Bin", "Day", "Hour")) #divided by 6 epochs per min #CHECK THIS, UNSURE
   data$Hour <- hourInfo$Hour
   realDay=1
   realDayArray=rep(0,length(data$Hour))
   for (i in 1:length(data$Hour)){
     if (i==1){
       realDayArray[i]=1
     }
     else if ((data$Hour[i-1]==23)&(data$Hour[i]==0)){
       realDay <- realDay+1
     }
     realDayArray[i]<- realDay
   }
   data$realDay=realDayArray
   sortedData <- data[
     with(data, order(Day, Hour)),
     ]
   numDays<- length(unique(data$Day))
   beginHour<-0
   pdf(file = paste0(here("Figs"),"/", mouse_name,"_actogram.pdf"),   # The directory you want to save the file in
       width = 6, # The width of the plot in inches
       height = 4) # The height of the plot in inches
   par(mfrow = c(numDays,1))
   par(mar = c(0.1, 2, 0.1, 2), oma=c(4, 4, 4, 4)) # Set the margin on all sides to 2
   ymax = 150
   xleft <- cbind(floor(cycleMappedOn)+1, (cycleMappedOn-floor(cycleMappedOn))*max(binsCol))
   
   #xleft = (cycleMappedOn[1])*max(binsCol)-1 # fraction of the 288 entries where the change occurs, -1 because we begin at 1 on the x axis
   #xright = ((cycleMappedOff+1)%%dayMod)*max(binsCol)-1
   xright <- cbind(floor(cycleMappedOff)+1, (cycleMappedOff-floor(cycleMappedOff))*max(binsCol)-1)
   ytop=1
   ybottom=0
   for (day in 1:numDays){
     whereDay <- as.numeric(mouse_df$timefactor)==day
     activeStart=mouse_df$activeFactor[whereDay][1]==1
     activeEnd=mouse_df$activeFactor[whereDay][sum(whereDay)]==1
     if (activeStart){
       xleft <- rbind(xleft, c(day,0))
     }
     if (activeEnd){
       xright <- rbind(xright, c(day,288))     
       }
   }
   #print(xleft)
   #print(xright)
   barFill = "black"
   bgFill = "gray70"
   #plot(c(1, 288), c(0, 1), type= "n", xlab = "", ylab = "", xaxt='n', ann=FALSE, yaxt='n', axes=FALSE)
  # rect(1, ybottom, 288, ytop, col=NULL)
   #rect(xleft, ybottom, xright, ytop, col="black")
   thres = 10 #meters
   thres= thres*6 #convert to "speed"
   excludedDays=c()
   for (day in 1:numDays){
     if (sum(data$`Sum Speed/5min`[data$Day==day])>=thres){
       barplot(data$`Sum Speed/5min`[data$Day==day]/6, ylim=c(0,ymax), xlim=c(1,288), space=0, col=barFill)
       axis(side=2, labels=FALSE, tick=FALSE)
       left <- xleft[,2][(xleft[,1])==day]
       right <- xright[,2][(xright[,1])==day]
       rect(left, ybottom, right, ymax, col=bgFill, border=bgFill)
       barplot(data$`Sum Speed/5min`[data$Day==day]/6, ylim=c(0,ymax), xlim=c(1,288), space=0, col=barFill, add=TRUE)
       #scale_x_continuous(expand = c(0, 0))+
       text(x=max(binsCol)-15,y=ymax-50, labels = paste0("Day ",day), cex = 0.75)
       #Sys.sleep(2) 
     }
     else{
       excludedDays = append(excludedDays, day)
     }
   }
   idx = match(mouse_name, micelist)
   myGeno = genotypes[idx]
   mtext(text=paste("Distance (m)"), outer=TRUE, side=2, line=1, cex = 0.6)
   mtext(text=paste0(mouse_name, " (", myGeno, ")"), outer=TRUE, side=1, cex = 0.7, line=1)
   dev.off()
   print(paste0(mouse_name, " excluded day: ", excludedDays))
 }
 #GRAPH 9: Check activity with Actogram
 quickActogram <- function(binSize = 5, micelist, genotypes){
   #takes mouse *name* for data frame object, bins data, and makes bar plot/actogram of that data across days
   #default bin size is 10 min
   pdf(file = paste0(here("Figs"),"/", "check_actogram.pdf"),   # The directory you want to save the file in
       width = 8, # The width of the plot in inches
       height = 11) # The height of the plot in inches
   par(mfrow = c(length(micelist),1))
   par(mar = c(0.1, 2, 0.1, 2), oma=c(4, 4, 4, 4)) # Set the margin on all sides to 2
   #1st # is 288 for 5 and 144 for 10, second number is 30 for 5 and 60 for 10
   for (m in 1:length(micelist)){
     mouse_name=micelist[m]
     mouse_df=get(paste0(mouse_name, "dt"))
     binsCol <- rep(1:288, each=30) #factor representing 5 min bins for 1 day
     epochBins = rep(binsCol, length.out = length(mouse_df$timefactor))
     data <- setNames(aggregate(mouse_df$speed, by=list(epochBins, mouse_df$timefactor), sum),c("Bin", "Day", "Sum Speed/5min")) #divided by 6 epochs per min #CHECK THIS, UNSURE
     numDays<- length(unique(data$Day))
     ymax = 150
     xleft = cycleMappedOn[1]*max(binsCol)-1 # fraction of the 288 entries where the change occurs, -1 because we begin at 1 on the x axis
     xright = cycleMappedOff[1]*max(binsCol)-1
     ytop=1
     ybottom=0
     barFill = "black"
     bgFill = "gray70"
     barplot(data$`Sum Speed/5min`[data$Day==numDays]/6, ylim=c(0,ymax), xlim=c(1,288), space=0, col=barFill)
     axis(side=2, labels=FALSE, tick=FALSE)
     rect(xleft, ybottom, xright, ymax, col=bgFill, border=bgFill)
     barplot(data$`Sum Speed/5min`[data$Day==numDays]/6, ylim=c(0,ymax), xlim=c(1,288), space=0, col=barFill, add=TRUE)
     text(x=max(binsCol)-60,y=ymax/2, labels = paste0("Mouse ",m,": ",mouse_name), cex = 1, adj=0)
   }
   mtext(text=paste("Distance (m)"), outer=TRUE, side=2, line=1, cex = 1)
   dev.off()
 }
 #GRAPH 10: Onset of running times to onset of active period
 photoperiodSync <- function(micelist, genotypes){
   #get first maxima in dataset for a given active period
   # what about reaching 50% of max for that day?
   cumActiveFactor <- cumsum(df$activeFactor[df$timefactor==day])
   maxSpeed <- max(df$speed[df$timefactor==day])
   if (!is.na(maxSpeed)){
     pass <- cumActiveFactor[df$speed>(maxSpeed/2)]
   }
   #make new DF: name, genotype, day, passTime
   #compare the time for that maxima to the time for the onset of the active period
   #plot the difference between onset of active period and first maxima for all mice
 }

 #Graph 1: Quick actogram to check activity--------------
 quickActogram(5,micelist,genotypes)
 
 #Graph 2: Individual mice smoothed----------------
 for (i in 1:length(micelist)){
    mouse <- micelist[i]
    genotype <- genotypes[i]
    smoothActivity(i, micelist, genotypes, cycleMapped, downSample = 15)
 }
 
 #Graph 3: Individual mice actogram------------------
 for (i in 1:length(micelist)){
   mouse <- micelist[i]
   actogram(mouse, binSize = 5, micelist, genotypes, cycleMappedOn, cycleMappedOff)
 }
 
 # Setup for grouped graphs (USER CAN CHANGE)--------------
 genoOptions = unique(genotypes)
 print(genoOptions) # these are the options (1...n) for the genotype names in your dataset
 group1Name = genoOptions[1] # change this line to change the name of the first group
 group2Name = genoOptions[4] # change this line to change the name of the second group
 
 # Grouped graphs-----------
 smoothActivityGroups(group1Name,group2Name, cycleMapped)
 smoothScaledActivityGroups(group1Name,group2Name, cycleMapped)
 # Graph 3: Avg speed per group per day 
 avgSpeedPerDayGroup(group1Name,group2Name)
 # Graph 4: % Time running per group per day 
 PercentTimeRunGroup(group1Name,group2Name)
 # Graph 5: Running bout length graphs
 BoutLengthGroup(group1Name,group2Name)