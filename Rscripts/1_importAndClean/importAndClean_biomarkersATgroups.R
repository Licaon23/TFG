###############################################################################
### IMPORT AND CLEAN SCRIPT                                                   #
###   ADNI patient biomarkers data                                            #
###   Version 1                                                               #
###   Desember 2022                                                           #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###     * Import biomarkers dataset:                                          #
###         - Assign correct variable format                                  #
###         - Quality control: identify missing values, correct range of      #
###             values and internal coherence.                                #
###         - Join information from other sources (biomarkers ATgroup)        #
###                                                                           #
###     
###############################################################################
source("Rscripts/libraries.R")


#==============================================================================
# 1) IMPORT BIOMARKERS DATASET
#==============================================================================

biomarkers <- read.table(here(rawdataDir,"amylode_status",
                              "biomarkers_csvfile_defNVT_MM.txt"),
                         header=T)

summary(biomarkers)
#------------------------------------------------------------------------------
# A) Quality control and validation exam
#------------------------------------------------------------------------------
  # Missing values have been detected to be coded as both NA and -99
  # PTAU measurments have an lower and upper precision limit. Values below
  #   and above these thresholds appear to be coded as "<" or ">".
  #   Since a numerical value is needed, minimum and maximum threshold are
  #   reported instead: 8 and 120.

vars <- c("ABETA","TAU","PTAU","ABETA_ex")
  # -99 as NA value
biomarkers[biomarkers==-99] <- NA

  # Correct lower and upper thresholds as numeric.
biomarkers[vars] <- lapply(biomarkers[vars], FUN=function(x){
                        as.numeric(gsub("^(>|<)([0-9]+)$","\\2",x))})

summary(biomarkers)
head(biomarkers)
str(biomarkers)

#==============================================================================
# 2) AMYLODE AND TAU STATUS. A/T groups
#==============================================================================

  #
  
  # Given the fact that some patients have two measurements, a mean value
  # will be calculated in order to asses what status they have in both
  # Amylode beta (ABETA_ex) and TAU status.
baselineBiomarkers  <- summaryBy(cbind(PTAU,ABETA_ex) ~ RID, data=biomarkers,
                                 FUN=mean,na.rm=T,keep.names = T)

  # Classify patient status in ABETA, TAU and A/T GROUP according to the 
  # following thresholds
  #   * ABETA_ex <= 1100 -> AB+ ; AB- otherwise
  #   * PTAU >= 24 -> T+ ; T- otherwise
  #   * A/T group: AB-/T-, AB-/T+, AB+/T-, AB+/T+

averagedBiomarkers <- within(averagedBiomarkers,{
  PTAU_cat <- cut(PTAU,
                  breaks=c(-Inf,24,+Inf),
                  right=F,
                  labels=c("T-","T+"))
  
  ABETA_ex_cat <- cut(ABETA_ex,
                      breaks=c(-Inf,1100,+Inf),
                      right=T,
                      labels=c("A+","A-"))
  
#  ATgroup <- as.factor(paste(ABETA_ex_cat,PTAU_cat,sep="/"))
  ATgroup <- as.factor(paste(ABETA_ex_cat,PTAU_cat,sep=""))
})

  # Result dataframe
patientATgroup <- averagedBiomarkers[c("RID","ABETA_ex_cat",
                                       "PTAU_cat","ATgroup")]

#save(patientATgroup,  file=here(cleandataDir,"patientATgroup.RData"))










#=======================================================================
# 3) COMPARE RESULTS WITH OFFICIAL CLASSIFICATION
#=======================================================================

# Grups oficials AT
at_oficial <- read.table(here(rawdataDir,"amylode_status",
                "ATgroup_df.txt"),header=T,na.strings = "")

at_oficial_cc <- at_oficial[complete.cases(at_oficial),]

#Valors ABETA i PTAU
biomarkersInfo <- lapply(split(biomarkers,biomarkers$RID),FUN="[",1,)
biomarkersInfo <- do.call(rbind,biomarkersInfo)

nrow(at_oficial)
nrow(biomarkersInfo)
setdiff(at_oficial_cc$PTID,biomarkersInfo$RID) # tots els ID amb grup tenen
                                               # info de biomarcadors.

# Tots els ptid amb atgroup estan al biomarcadors?
merge <- merge(biomarkersInfo,at_oficial_cc,by.x="RID",by.y="PTID",
               all=T)

ptid_withInfo_andNoATgroup <- subset(merge,is.na(ATgroup))
nrow(ptid_withInfo_andNoATgroup)
ptid_noInfo <- sprintf("%04d",ptid_withInfo_andNoATgroup$RID)
# 374 ptid dels quals hi ha info PTAU i ABETA_ex, pero no grup AT...
