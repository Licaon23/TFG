###############################################################################
### IMPORT AND CLEAN SCRIPT                                                   #
###   ADNI participants data for analysis                                      #
###   Version 1                                                               #
###   November-December 2022                                                  #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed: 
###     * Select target dataset for analysis (3T, 1.5T or mixed)              #
###     * Select variable of interest, calculate some variables and create    #
###         a basal characteristics dataset for cross-sectional analysis.     #
###     * RESULT: list with two dataframes: longitudinal data and baseline    #
###############################################################################
source("Rscripts/libraries.R")
load(here(cleandataDir,"volumesAndPRS.RData")) # load data


### FUNCIO ====================================================================
ready2Analysis <- function(dd, codebook, saveName){
  
  # Elimianr els PRS que no farem servir ---------------------------------------
  prsVars <- names(dd)[grep("^(PRS_)(.*)$",names(dd))]
  work_prs <- paste("PRS",c("AD","ADnoAPOE","FTD","IEAA","EEAA","LONGEVITY",
                            "LONGEVITYnoAPOE"),
                    sep="_")
  deletePRS <- prsVars[!(prsVars %in% work_prs)]
  aux <- !(names(dd) %in% deletePRS)
  
  dd <- dd[,aux]
  
  # Calcular PRS variables com a categòriques. percentil 80 ---------------
  aux <- grep("^(PRS_)(.*)$",names(dd))
  prs_names <- sub("^(PRS_)(.*)$","\\2",names(dd)[aux])
  PRS_catVars <- paste(prs_names,"predisposition",
                       sep="_")
  dd[PRS_catVars] <- lapply(dd[aux], FUN=function(x){
    threshold <- quantile(x,probs = 0.8)
    PRS_cat <- cut(x,breaks = c(-Inf,threshold,+Inf),labels = c("Low-Intermediate","High"))
    return(PRS_cat)
  })
  
  
  # Dades de característqiues basals ----------------------------------------
  deleteVars <- c("enrolDate","tesla","measureDate","VISCODE2",
                  "timelapse","time_age","time_visit" )
  deleteVars <- which(names(dd) %in% deleteVars)
  basal <- lapply(split(dd,dd$PTID),FUN="[",1,-deleteVars)
  basal <- do.call(rbind,basal)
  

  
  # Add labels codebook for categorical PRS
  names <- names(dd)[grep("predisposition",names(dd))]
  labels <- sub("^(.*)(_.*)$","\\1",names)
  aux <- data.frame(Variable = names, Label = labels,
                    description = rep("PRS categorized as low or high predisposition, from percentil 80",
                                      length(names)))
  
  codebook <- rbind(codebook,aux)
  
  # Save clean results
  name <- paste0("dataForAnalysis_",saveName,".RData")
  save(dd,basal,codebook,
       file = here(cleandataDir,name))
  
}
#===============================================================================

ready2Analysis(dd = data_singleT_overall$data$T_1.5,
               codebook = data_singleT_overall$codebook,
               saveName = "1T")

ready2Analysis(dd =data_singleT_overall$data$T_3,
               codebook = data_singleT_overall$codebook,
               saveName = "3T")

ready2Analysis(dd = data_singleT_mixed$data,
               codebook = data_singleT_mixed$codebook,
               saveName = "mixed")




rm(list=ls())
