###############################################################################
### IMPORT AND CLEAN SCRIPT                                                   #
###   ADNI participants Genetic Scores                                        #
###   Version 1                                                               #
###   December 2022                                                           #  
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed: 
###    * Import genetic PRS data from independent .txt files.
###    * Data cleaning and merge all of them in a single dataset
###                                                                           #
###############################################################################

source("Rscripts/libraries.R")
source(here(RscriptsDir,"1_importAndClean","importAndClean_auxFunctions.R"))


#=======================================================================
# A) IMPORT SCORE DATA
#=======================================================================

### Files names and paths----------------------------------
path <- here(rawdataDir,"scores_genetics")
aux <- grep(list.files(path), pattern="^(.*)\\.txt") # select .txt files
files <- list.files(path,full.names = T)[aux]

### FUNCTION: import data-------------------------------
read_geneticScore <- function(file){
  #--------------------------------------
  # READ DATA
  #--------------------------------------
  score <- read.table(file, header=T)
  
  #-----------------------------------------------------------
  # Replace PTID current value by only its last 4-digit code 
  #-----------------------------------------------------------
  regex_PTID <- "^(\\d{3})_(\\w)_(\\d{4})$" # expected format
  
  #any(!grepl(regex_PTID,score$IID)) #is there any PTID with a different format?
  # we are good;all PTID are consistent format-wise.
  score$PTID <- as.factor(gsub(regex_PTID,"\\3",score$IID))
  
  #------------------------------------------------------
  # Keep Pt_5e.06 column and PTID
  #------------------------------------------------------
  aux <- names(score)[grep("(*.)5e.06$",names(score))]
  score <- score[,c("PTID",aux)]
  
  #------------------------------------------------------
  # Get target name
  #------------------------------------------------------
  aux <- tail(strsplit(file,split="_")[[1]],1)
  targetName <- sub("(*.)\\.txt","\\1",aux)
  colnames(score) <- c("PTID",paste(targetName,"Pt_5e.06",sep="_"))
  
  # RETURN
  return(score)
}

### Apply import function to all files, and merge them by PTID -------------
ptid_scores <- lapply(files[-1],read_geneticScore)
ptid_scores <- Reduce(function(x,y){merge(x,y,by="PTID")},ptid_scores)

### Apply import for AD files (differnet format) ---------------------
ptid_ADScores <- read_geneticScore(files[1])[1:2]
names(ptid_ADScores) <- c("PTID","ADnoAPOE_Pt_5e.06")
ptid_scores <- merge(ptid_scores,ptid_ADScores,by="PTID")


### Rename variables and reorder ----------------------------------------
prs_names <- sub("(*.)(_Pt_5e.06)$",replacement = "\\1",
                 names(ptid_scores)[-1])
names(ptid_scores)[-1] <- prs_names
ptid_scores <- ptid_scores[,c("PTID","AD","ADnoAPOE","FTD","PKSON",
                              "IEAA","EEAA","TELOM",
                              "FRAILTY","LEXP","LONGEVITY",
                              "LONGEVITYnoAPOE")]

names(ptid_scores)[-1] <- paste("PRS",names(ptid_scores)[-1],sep="_")

#=======================================================================
# B) Normalize score value: Z-score center+scale
#=======================================================================

ptid_scores[,-1] <- lapply(ptid_scores[,-1],scale)
ptid_scores[,-1] <- lapply(ptid_scores[,-1],as.numeric)

#=======================================================================
# C) Histograms
#=======================================================================

file <- here("results","Others","PRS_histograms.png")
png(file,width = 800,height = 480)
par(mfrow=c(3,5))
aux <- ptid_scores[,-1]

lapply(1:ncol(aux),function(i){
  h <- hist(aux[,i], xlab="Score value",ylab="",main=colnames(aux)[i],freq=F)
})

dev.off()



#=======================================================================
# C) SAVE clean dataset
#=======================================================================
save(ptid_scores,
     file=here(cleandataDir,"patient_PRS.RData"))

rm(list=ls())
