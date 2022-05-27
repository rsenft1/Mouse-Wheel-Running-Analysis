### R scripts for loading and cleaning data----

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
  
  data.frame(time=myTime, speed=z, expDay = factor(ceiling(myTime))) #basic df
}
addTimeInfo <- function(df, myTime, binSize=10, lightsOnTime=7, lightsOffTime= 19){
  #now add in time information
  #Current time
  currentHour <- as.numeric(substring(myTime,1,2))
  currentExpHour <- 1
  currentMinute <- as.numeric(substring(myTime,4,5))
  nRows <- length(df$time)
  #initialize variables that will be hours and whether it's the animal's active cycle
  Hours <- replicate(nRows,0)
  expHours <- replicate(nRows,0)
  activeFactor <- replicate(nRows,0)
  minute <- 60/binSize #num of entries for each minute
  byMin <- seq(from=1,to=nRows, by=60/binSize)
  for (k in byMin){
    Hours[k:(k+minute-1)] <- currentHour
    expHours[k:(k+minute-1)] <- currentExpHour
    currentMinute <- currentMinute+1
    currentMinute <- currentMinute%%60
    if (currentMinute==0){
      currentHour <- currentHour+1
      currentExpHour <- currentExpHour+1
    }
    currentHour <- currentHour%%24
  }
  #add hours to df
  Hours <- head(Hours, nRows) #fix extra entries added at the end by going in minutes
  expHours <- head(expHours, nRows) #fix extra entries added at the end by going in minute
  df$hour <- Hours
  df$expHour <- expHours
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
  day <- df$expDay
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
  filelist = list.files(directory, pattern = ".*.log")
  myTime <- scan(here(directory, filelist[1]), '', skip = 2, nlines = 1, sep = '\n')
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
  mice = list.files(directory)
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
filter.data <- function(df, thres=20){
  temp <- df
  thres= thres*6 #convert meters to "speed"
  for (d in 1:max(as.numeric(temp$expDay))){
    if (sum(temp$speed[temp$expDay==d])<thres){
      temp$speed[temp$expDay==d]<- NA
      temp$smoothSpeed[temp$expDay==d]<- NA
    }
  }
  return(temp)
}