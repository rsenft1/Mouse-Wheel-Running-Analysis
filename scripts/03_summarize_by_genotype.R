# Create data tables based on genotype --------------------------------------
rowNum <- length(df$speed)
allGeno = unique(genotypes)
all_geno_df <- data.frame(genotype = character(0), time = double(0), expDay = integer(0), avg.raw.speed = double(0), sem.raw = double(0), avg.smooth.speed = double(0), sem.smooth = double(0), avg.scaled.raw.speed = double(0), sem.scaled.raw.speed = double(0), avg.scaled.smooth.speed = double(0), sem.scaled.smooth = double(0))
all_geno_avgSpeedPerDay <- data.frame(Day = integer(0), Mouse = character(0), average.raw.speed=double(0), genotype = character(0))
all_geno_PercentRunningPerDay <- data.frame(Day = integer(0), Mouse = character(0), PercentRunning=double(0), genotype = character(0))

for (geno in allGeno){
  myMice <- micelist[genotypes==geno]
  indSpeedMatrix <- matrix(data=NA, nrow=rowNum, ncol=length(myMice)) #initialize
  colnames(indSpeedMatrix) <- myMice
  indSmoothSpeedMatrix <- indSpeedMatrix
  #avgSpeedDf <- data.frame(genotype = character(0), time = double(0), expDay = integer(0), avg.raw.speed = double(0), sem.raw = double(0), avg.smooth.speed = double(0), sem.smooth = double(0), avg.scaled.raw.speed = double(0), sem.scaled.raw.speed = double(0), avg.scaled.smooth.speed = double(0), sem.scaled.smooth = double(0))
  mouseNum <- 0
  for (mouse in myMice){
    df <- get(paste0(mouse,"dt"))
    indSpeedMatrix[,mouse] <- df$speed
    indSmoothSpeedMatrix[,mouse] <- df$smoothSpeed
    mouseNum <- mouseNum +1
  }
  scaledRawSpeedMat <- apply(indSpeedMatrix, 2, scaleCol)
  scaledSmoothSpeedMat <- apply(indSmoothSpeedMatrix, 2, scaleCol)
  
  #smooth and scale speed for each mouse in matrix
  avgSpeedDf = data.frame(
    time = df$time,
    genotype = rep(geno, length(df$time)),
    expDay = df$expDay,
    avg.raw.speed = rowMeans(indSpeedMatrix, na.rm=TRUE),
    sem.raw = rowSems(indSpeedMatrix),
    avg.smooth.speed = rowMeans(indSmoothSpeedMatrix, na.rm=TRUE),
    sem.smooth = rowSems(indSmoothSpeedMatrix),
    avg.scaled.smooth.speed = rowMeans(scaledSmoothSpeedMat),
    sem.scaled.smooth.speed = rowSems(scaledSmoothSpeedMat),
    avg.scaled.raw.speed = rowMeans(scaledRawSpeedMat),
    sem.scaled.raw.speed = rowSems(scaledRawSpeedMat) 
  )
  #get avg speed per day for individuals of a given genotype
  avgPerDayInd <- setNames(aggregate(indSpeedMatrix,by=list(df$expDay), mean), c("Day", colnames(indSpeedMatrix)))
  avgPerDayInd_with_geno <- melt(avgPerDayInd ,  id.vars = 'Day', variable.name = 'Mouse', value.name = 'average.raw.speed')
  avgPerDayInd_with_geno$genotype <- rep(geno, length(avgPerDayInd_with_geno$Day))
  
  #get % of time running per day for individuals of a given genotype. 8640 is the number of 10 second epochs in a day. We look for % of those with registered activity
  percentRunningInd <- setNames(aggregate(indSpeedMatrix,by=list(timefactor), function(x){(sum(x>0)/8640)*100}), c("Day", colnames(indSpeedMatrix)))
  percentRunningInd_with_geno <- melt(percentRunningInd ,  id.vars = 'Day', variable.name = 'Mouse', value.name = 'PercentRunning')

  #save the tables
  write.csv(indSpeedMatrix, here(cleanDataDir, paste0(str_replace_all(geno,"/", "-"), "_ind_speeds.csv")), row.names = FALSE)
  write.csv(avgSpeedDf, here(cleanDataDir, paste0(str_replace_all(geno,"/", "-"), "_subtable.csv")), row.names = FALSE)
  write.csv(avgPerDayInd, here(speedEachDayDir, paste0(str_replace_all(geno,"/", "-"), "_avg_speed_per_day.csv")), row.names = FALSE)
  write.csv(percentRunningInd, here(cleanDataDir, paste0(str_replace_all(geno,"/", "-"), "_percent_time_running_per_day.csv")), row.names = FALSE)
  
  
  #get sd and from sd, the sem of the smoothed and scaled data
  all_geno_df <- rbind(all_geno_df, avgSpeedDf)
  all_geno_avgSpeedPerDay <- rbind(all_geno_avgSpeedPerDay, avgPerDayInd_with_geno)
  all_geno_PercentRunningPerDay <- rbind(all_geno_PercentRunningPerDay, percentRunningInd_with_geno)
  
}
#combine all genotype averages into one csv
mean_sem <- function(x) c(Mean = mean(x), SEM = sd(x)/sqrt(sum(!is.na(x))))

#groupAvgPerDay <- do.call(data.frame, aggregate(average.raw.speed~Day, meltAvgPerDayInd, semBasic))
#groupAvgPerDay <- do.call(data.frame, aggregate(avg.raw.speed~expDay, all_geno_df, semBasic))
all_genoAvgPerDay <- do.call(data.frame, aggregate(avg.raw.speed~expDay+genotype, all_geno_df, mean_sem)) #check this one. sem is wrong based on # of points in a day rather than # of animals
all_geno_avgSpeedPerDay <- do.call(data.frame, aggregate(average.raw.speed~Day+genotype, all_geno_avgSpeedPerDay, mean_sem))
all_geno_PercentRunningPerDay <- do.call(data.frame, aggregate(PercentRunning~Day+genotype, all_geno_PercentRunningPerDay, mean_sem))

#Save these group tables with all genotypes together
write.csv(all_geno_df,here(groupCleanDataDir, "all_geno_avg_speed_all_points.csv"), row.names = FALSE)
write.csv(all_genoAvgPerDay,here(speedEachDayDir, "all_geno_avg_speed_per_day.csv"), row.names = FALSE)
write.csv(all_geno_PercentRunningPerDay,here(groupCleanDataDir, "all_geno_percent_running_per_day.csv"), row.names = FALSE)


