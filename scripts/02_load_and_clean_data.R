# Step 3. Load data into data frames --------------------------------------
loadAndCleanData <- function(dataPath, resultsPath, micelist, myTime, toFilter=TRUE, dayofChange=0, newlightsOnTime=7, newlightsOffTime= 12, oldlightsOnTime = 7, oldlightsOffTime=19){
  for (i in 1:(length(micelist))){
    direct <- here(dataPath,micelist[i])
    name <- paste(micelist[i],"dt", sep="")
    percentDone <- (i-1)/length(micelist)*100
    print(paste0("Loading data for Mouse ",i, " ... ", round(percentDone), "% complete"))
    temp <- getData(direct,binSize=10)
    temp$smoothSpeed <- smo.dist(temp$speed, 200)
    temp <- addTimeInfo(temp, myTime)
    if(dayofChange != 0){
      temp <- addTimeChange(temp, dayofChange, binSize=10, newlightsOnTime=7, newlightsOffTime= 12, oldlightsOnTime = 7, oldlightsOffTime=19)
    }
    if (toFilter==TRUE){
        temp <- filter.data(temp)
    }
    assign(name,temp)
    write.csv(get(name),here(resultsPath, paste0(name, ".csv")), row.names = FALSE)
  }
}
