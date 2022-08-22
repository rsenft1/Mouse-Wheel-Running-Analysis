addSpeedTime <- function(df, binSize=10){ # binSize in seconds
  # conversion factor for meter/min, 0.0625 is radius in meters (=6.25 cm)
  df <- df %>% mutate(speed = raw*2*pi*0.0625 *60) # not bin size dependent 	
  
  # create column for elapsed time in days, beginning of experiment = 0
  df <- df %>% mutate(elapsedDays = (1:length(speed)) / 3600 * binSize / 24) 
  
  return(df)
}

addHoursActive <- function(mouseName, metadata, binSize=10, lightsOnTime=7, lightsOffTime= 19){
  # This function adds columns to the df for clock time and subjective/relative time. 
  # Clock time is relative to the local midnight (midnight = 0).
  # Subjective/relative time is relative to the beginning of the first light cycle (i.e., first
  # lights ON) on the first experimental day (this typically happens before the mice are put in
  # the cages). 
  
  # Variables:
  # mouseName:     str, name of mouse(e.g., "C676")
  # metadata:      dataframe, the metadata dataframe for the experiment
  # binSize:       int, how frequently measurements are taken in seconds
  # lightsOnTime:  int, hour corresponding to lights on
  # lightsOffTime: int, hour corresponding to lights off
  
  # Get the time the experiment started recording from the mouse metadata table
  df <- get(mouseName)
  myTime <- metadata %>% filter(MouseName==mouseName) %>% select(StartTime) %>% pull()
  currentHour <- as.numeric(substring(myTime,1,2))
  currentExpHour <- 1
  currentMinute <- as.numeric(substring(myTime,4,5))
  currentSeconds <- behavr::hours(currentHour)+behavr::mins(currentMinute)
  nRows <- length(df$raw)
  
  # set t time (in seconds) relative to lights on, set such that the last lights on before starting measurements is t=0
  t <- seq(currentSeconds-behavr::hours(lightsOnTime), nRows*10+currentSeconds-10-behavr::hours(lightsOnTime), by=binSize)
  t2 <- seq(currentSeconds, nRows*10+currentSeconds-10, by=binSize)
  
  df <- df %>% mutate(expDay = ceiling((t+behavr::hours(lightsOnTime)+10)/behavr::days(1)))
  df <- df %>% mutate(hour = floor((t2)/behavr::hours(1)) %% 24)
  df <- df %>% mutate(minutes = floor((t2)/behavr::mins(1)) %% 60)
  
  activeFactor <- rep(makeDayActivityWrapper(lightsOnTime, lightsOffTime, startSecs=t[1]), ceiling(max(df$expDay)))[seq(1,length(df$raw))]
  #add active factor to df
  df$activeFactor <- activeFactor
  
  #add t for behavr
  df$t <- t
  return (df)
}

addExpPhase <- function(df, phaseNames, phaseStart){
  # phaseNames = array of names for phases in experiment, set in settings
  # phaseStart = array of days each phase starts on
  numDays <- max(df$expDay)
  phases <- c(phaseNames, "end")
  phaseStart <- c(phaseStart, numDays+1)
  currentPhase <- rep(NA, numDays)
  for(i in seq(1, length(phases)-1)){
    start <- phaseStart[i]
    end <- phaseStart[i+1]
    currentPhase[seq(start, end-1)] <- rep(phases[i], abs(start-end))
  }
  phase_df <- data.frame(expDay = 1:numDays, expPhase = currentPhase)
  df <- merge(x = df, y = phase_df, by = "expDay") 
  return(df)
}

makeDayActivityWrapper <- function(lightsOn, lightsOff, startSecs, originalLightsOn=7, binSize=10){
  #lightsOn = hour lights come on, relative to midnight = 0
  #lightsOff = hour lights go off, relative to midnight = 0
  #binSize = length of time between measurements, in seconds
  #startSecs = starting time for recording data, in seconds
  print(startSecs)
  startHour <- startSecs/behavr::hours(1)+originalLightsOn # Start time in hrs
  print(startHour)
  if(min(lightsOn, lightsOff)==lightsOn){
    if(0 < lightsOn){
      start <- T
    }else{
      print("start is false")
      start <- F
    }
  }
  if(min(lightsOn, lightsOff)==lightsOff){
    if(0 < lightsOff){
      start <- F
    }else
    {start <- T}
  }
  print(start)
  first <- min(lightsOn, lightsOff)
  second <- abs(lightsOn-lightsOff)
  third <- 24 - sum(first, second)
  makeDayActivity(hours=c(first, second, third), start=start, startSecs=startSecs+behavr::hours(originalLightsOn), binSize=binSize)
}

makeDayActivity <- function(hours, startSecs, start=T, binSize=10){
  #hours = vector of hours alternating active and inactive (e.g., c(7, 12, 5))
  #start = whether to start the day (midnight) in the active or inactive state, T=start active
  #binSize = length of time between measurements, in seconds
  
  if (sum(hours) !=24){
    stop(message=paste("Hours must sum to 24. They currently sum to",sum(hours)))
  }
  activeState <- start
  myDay <- c()
  for(h in hours){
    myDay<- append(myDay, rep(activeState, behavr::hours(h)*1/binSize))
    activeState <- !activeState
  }
  myDay <- rep(myDay, 2) #double up the day
  start <- startSecs/binSize+1 # Start time in hrs
  
  #start <- (startSecs/binSize)+1 #get index
  myDay<- myDay[seq(start, start+(behavr::days(1)/binSize)-1)]
  return(myDay)
}


addTimeChange <- function(df, dayOfChange, dayChangeEnds = F, binSize=10, newlightsOnTime=newlightsOnTime, newlightsOffTime= newlightsOffTime){
  #binSize = length of time between measurements, in seconds
  currentSeconds <- df$t[1]
  newDay <- makeDayActivityWrapper(newlightsOnTime, newlightsOffTime, startSecs=currentSeconds)
  if(dayChangeEnds==F){
    dayChangeEnds <- max(as.numeric(df$expDay))
  }
  newDays <- rep(newDay, dayChangeEnds-dayOfChange)
  changeTheseTimepoints <- seq((behavr::days(dayOfChange))*(1/binSize)+1, (behavr::days(dayChangeEnds))*(1/binSize)-1)
  if(max(changeTheseTimepoints > nrow(df))){
    #if the timepoints to change extend beyond the number of rows in the dataset, chop off the data
    #this would happen if the last day isn't a full 24 hrs long
    subsetindex <- c(changeTheseTimepoints <= nrow(df))
    changeTheseTimepoints <- changeTheseTimepoints[subsetindex]
    newDays <- newDays[subsetindex]
  }
  df$activeFactor[changeTheseTimepoints] <- newDays
  return(df)
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

everyThird <- function(x){
  x[seq(2, length(x), 3)] <- ""
  x[seq(3, length(x), 3)] <- ""
  x
}