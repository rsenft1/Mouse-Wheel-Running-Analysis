actogram <- function(mouse_name, micelist, genotypes, cycleMappedOn, cycleMappedOff, binSize = 5){
  #takes mouse *name* for data frame object, bins data, and makes bar plot/actogram of that data across days
  #default bin size is 10 min
  fileName <- paste0(mouse_name,"_actogram")
  #1st # is 288 for 5 and 144 for 10, second number is 30 for 5 and 60 for 10
  mouse_df=get(paste0(mouse_name, "dt"))
  binsCol <- rep(1:288, each=30) #factor representing 5 min bins for 1 day
  epochBins = rep(binsCol, length.out = length(mouse_df$expDay))
  Mode <- function(vector) {
    choices <- unique(vector)
    choices[which.max(tabulate(match(vector, choices)))]
  }
  data <- setNames(aggregate(mouse_df$smoothSpeed, by=list(epochBins, mouse_df$expDay), sum),c("Bin", "Day", "Sum Speed/5min")) #sum of speeds for that epoch (sum)
  hourInfo <- setNames(aggregate(as.numeric(mouse_df$hour), by=list(epochBins, mouse_df$expDay), Mode),c("Bin", "Day", "Hour")) #most common hour per epoch (mode)
  timeInfo <- setNames(aggregate(mouse_df$time, by=list(epochBins, mouse_df$expDay), max),c("Bin", "Day", "Time")) #max time at end of each epoch (max)
  data$Hour <- hourInfo$Hour
  data$time <- timeInfo$Time
  # FINDING onset and offset of activity for each day...
  #r <- rle(data$`Sum Speed/5min`>0) # any activity
  activityThres <- .01*max(data$`Sum Speed/5min`, na.rm=TRUE)
  r <- rle(data$`Sum Speed/5min`>activityThres)
  moving <- unlist(r[1])[unlist(r[2])==TRUE] #this takes the lengths of active bouts (when r[2] is true, the animal's speed is >0) and converts 10 second observations to minutes
  both_move_and_not <-unlist(r[1])
  bout_borders <- cumsum(both_move_and_not)
  moving_ind<-unlist(r[2])
  start_bout <- c()
  stop_bout <- c()
  epochThres <- 10
  for(epoch in 1:(length(both_move_and_not)-1)){
    if(!moving_ind[epoch] & (both_move_and_not[epoch]>epochThres)){
      if(moving_ind[epoch+1]&(both_move_and_not[epoch+1]>epochThres)){
        start_bout <- append(start_bout, as.numeric(bout_borders[epoch]))
      }
    }
    if(moving_ind[epoch] & (both_move_and_not[epoch]>epochThres)){
      if(!moving_ind[epoch+1]&(both_move_and_not[epoch+1]>epochThres)){
        stop_bout <- append(stop_bout, as.numeric(bout_borders[epoch]))
      }
    }    
  }

  
  start_time <- data$time[start_bout]
  start_hour <- data$Hour[start_bout]
  
  stop_time <- data$time[stop_bout]
  stop_hour <- data$Hour[stop_bout]
  
  plot(data$time[data$Day==c(1:10)], data$`Sum Speed/5min`[data$Day==c(1:10)], type="l")
  points(stop_time, rep(activityThres, length(stop_time)), col="red")
  points(start_time, rep(activityThres, length(start_time)), col="green")
  #plots starts and stops on the activity graph!
  
  
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
    whereDay <- as.numeric(mouse_df$expDay)==day
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