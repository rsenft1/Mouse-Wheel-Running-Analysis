# R Script to organize your files!
#This package-checking chunk of code is from vikram: https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
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

library(here)
#dir_rawData <- here("data", "Raw_data") #name of folder within project directory containing the raw data files (.log) USER MUST UPDATE

#1. Get names of the .log files:
logfiles <- list.files(here(dir_rawData), full.names=FALSE)

#2. Get names of mice from log file names
mouseNames <- sapply(str_split(logfiles, pattern="_"), getElement,1) #names are before first underscore. Don't use underscores in the mouse IDs or this will not work

#3. Make folders based on unique names
folderNames <- unique(mouseNames)
#dir_orgData <- here("data", "Organized_data") #uncomment this line if not defined in top-level script
for (i in seq_along(folderNames)){
  if (!file.exists(here(dir_orgData,folderNames[i]))){
    folder<-dir.create(here(dir_orgData,folderNames[i]))
  }
}

#4. Now move the .log files into the appropriate folders!
for (i in seq_along(logfiles)){
  file.copy(here(dir_rawData, logfiles[i]), here(dir_orgData,mouseNames[i]))
}

