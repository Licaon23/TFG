###############################################################################
### IMPORT AND CLEAN SCRIPT                                                   #
###   ADNI participants demographic data                                      #
###   Version 1                                                               #
###   October-November 2022                                                   #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###     * Import demographic dataset:                                         #
###         - Assign correct variable format                                  #
###         - Quality control: identify missing values, correct range of      #
###             values and internal coherence.                                #
###         - Join information from other sources (biomarkers ATgroup and PRS)#
###         - Reclassify diagnositc groups according                          #
###                                                                           #
###     * Analysis ADNI path for each patient.                                #
###         - Get in-out-follow up values                                     #
###                                                                           #
###     * Unique patient info dataset: create dataframe with unique           #
###         demographic information for each patient, excluding redundant     #
###         multiple visits.                                                  #
###############################################################################
source("Rscripts/libraries.R")
load(here(cleandataDir,"patient_PRS.RData"))

# =============================================================================
# 1) IMPORT DEMOGRAPHIC DATA FROM PACIENTS
# =============================================================================

# It's been noticed at first glimpse that missing values appear encoded both as
# NA and "". What is more, it seems that several cases are corrupted,
# with strange values in some fields and even less fields than other registers.

file <- list.files(here(rawdataDir))
file <- file[grep("demographic",file)]
demo <- read.table(here(rawdataDir, file),
                   sep=" ", header=T, fill=T, na.strings = c("",NA))
rm(file)

# -----------------------------------------------------------------------------  
# A) Assign correct format to variables
# -----------------------------------------------------------------------------  
  # i) RID as integer, to get rid of corrupted values and set them to NA
demo$RID <- as.integer(demo$RID) #RID com a enter

  # ii) EXAMDATE as Date format
demo$EXAMDATE <- as.Date(demo$EXAMDATE, format="%d/%m/%Y") # a format data

  # iii) String as factors
catvars <- sapply(demo,is.character) # caracters com a factor
demo[catvars] <- lapply(demo[catvars],as.factor)

# -----------------------------------------------------------------------------  
# B) Quality control and cleaning procedures
# -----------------------------------------------------------------------------  

### RID ### ---------------------------------

  # There appear to be 28 cases with missing RID value. These are corrupted
  # cases from raw data, as all the other fields don't have values either.
  # Therefore, they will be removed from the dataset
#summary(demo)
#demo[is.na(demo$RID),]
demo <- demo[!is.na(demo$RID),]


### PTID ### ---------------------------------

  # On the other hand, it seems that there are several cases with RID value
  # but not PTID. All of them have nothing but NA values. They are to be 
  # removed as well.
#summary(demo)
#summary(demo[is.na(demo$PTID),]) 
demo <- demo[!is.na(demo$PTID),]


  # Since many empty factor values have been removed, it is necessary
  # to drop their empty levels as well from the factor object.
demo[catvars] <- lapply(demo[catvars],droplevels)


  # With regard to patient identifier, PTID, check whether all of them
  # have the same format
ptid <- demo$PTID
regex_PTID <- "^(\\d{3})_(\\w)_(\\d{4})$" # expected format
any(!grepl(regex_PTID,ptid)) #is there any PTID with a different format?
                             # we are good;all PTID are consistent format-wise.

  # Replace PTID current value by only its last 4-digit code
demo$PTID <- as.factor(gsub(regex_PTID,"\\3",demo$PTID))
#head(demo)



### VISCODE2 ### ----------------------------------------

  # Factor variable that points out when a measure was taken, whether it is
  # a baseline measure or it is a follow-up measure after 6,12,18... months.
  # Since it is an ordinal factor, it must be encoded as such.
viscode_levels <- levels(demo$VISCODE2)
viscode_regex <- "(m)(\\d+)"
measure_levels <- grepl(viscode_regex,viscode_levels) #follow-up measures
order_measures <- sub(viscode_regex,"\\2",
                      viscode_levels[measure_levels]) #extract only months
order_measures <- order(as.numeric(order_measures)) #as numeric and order

  # Sort viscode levels acording to its natural time order
viscode_levels <- c(viscode_levels[!measure_levels],
                    viscode_levels[measure_levels][order_measures])

demo$VISCODE2 <- ordered(demo$VISCODE2,levels=viscode_levels)

### Cohort indicator: ORIGPROT, COLPROT ### --------------------------

  # ADNI project has several cohorts, with different patients
  # going in at different times, and in some cases leaving the project
  # at the end of a phase. 
  # These different phases are: ADNI1, ADNIGO, ADNI2 and ADNI3.
  # ORIGPORT indicates which cohort a given patient was enrolled at,
  # and COLPROT points out during what cohort a given measure was taken.
  # Levels from both variables should reflect this time ordering.
adni_levels <- c("ADNI1","ADNIGO","ADNI2","ADNI3")
demo[c("COLPROT","ORIGPROT")] <- lapply(demo[c("COLPROT","ORIGPROT")],
                                        FUN=ordered,levels=adni_levels)


### Add AT-group information ### -----------------------------------------

atgroup <- read.table(here(rawdataDir,"amylode_status",
                         "ATgroup_df.txt"),
                    header=T,na.strings = "")

# Descartem pacients sense info
#atgroup <- atgroup[complete.cases(atgroup),]

#Merge
demo <- merge(demo,atgroup, by.x = "RID", by.y = "PTID", all.x=T)
demo$ATgroup <- as.factor(demo$ATgroup)

### Reclassify SMC as Control patients CN and ### ------------------------- 
### merge EMCI and LMCI into MCI category     ### -----------------------------
  
  # NOTE: SMC to CN because proportion of AB+ is greater in CN than SMC
  # Output crosstable DX.bl groups and amyloid status
    #     |   AB+         AB-
    #  -----------------------------
    # CN  | 0.435(124)  0.565(161)
    # AD  | 0.927(216)  0.073(17)
    # SMC | 0.368(35)   0.632(60)
    # MCI | 0.677(436)  0.323(208)

  # Add new factor levels
  levels(demo$DX.bl) <- c(levels(demo$DX.bl), "MCI") 

  # Reclassify
  demo[demo$DX.bl %in% c("EMCI","LMCI"),"DX.bl"] <- "MCI"
  demo[demo$DX.bl %in% c("SMC"),"DX.bl"] <- "CN"
  
  # Drop empty levels
  demo$DX.bl <- droplevels(demo$DX.bl)
  
  # Revel CN control as reference level and MCI as intermediate status.
  demo$DX.bl <- factor(demo$DX.bl, levels = c("CN","MCI","AD"))
  demo$DX.bl <- relevel(demo$DX.bl, ref="CN")


### Adding patient genetic information available ### ------------------------- 
### Polygenic risk scores                        ### -------------------------
demo <- merge(demo,ptid_scores,by="PTID",all.x=T)

### Ordering cases by PTID and EXAMDATE ### ------------------------------

demo <- demo[order(demo$PTID,demo$EXAMDATE),]
rownames(demo) <- NULL
# head(demo)
# summary(demo)
# str(demo)



# =============================================================================
# 2) PACIENT TRAJECTORY THROUGH ADNI PROJECT
# =============================================================================

# GOAL: to know for each pacient in the dataset from which project they started
#       and what their trajectory is thorugh the different project phases.

#---------------------------------------------------------------------
# i) Data wrangling: for each patient, get one observation for each
#     project during which they have been measured.
#---------------------------------------------------------------------
aux <- duplicated(demo[,c("PTID","COLPROT")])
patientCohorts <- demo[!aux,]
#head(patientCohorts)


#---------------------------------------------------------------------
# ii) For each patient, get their path through ADNI phases
#---------------------------------------------------------------------

# FUNCTION: get_PatientPath
# Given all patient rows from 'patientCohorts' dataset, get a string
# with a sequence of ADNI phases they have gone through
get_PatientPath <-  function(x){
  cohorts <- x$COLPROT
  path <- paste0(cohorts,collapse="-")
  return(path)
}

# Get ADNI paths for each patient
aux <- split(patientCohorts, f = patientCohorts$PTID)
patientPaths <- lapply(aux,FUN=get_PatientPath)
patientPaths <- do.call("rbind", patientPaths)
patientPaths <- data.frame(PTID = rownames(patientPaths),
                           Path=as.factor(patientPaths[,1]),
                           row.names = NULL)
#head(patientPaths)

# Get number of patients with the same trajectory through ADNI project
unique.paths <- table(patientPaths$Path)
#unique.paths

save(patientPaths,patientCohorts,unique.paths,
     file=here(cleandataDir,"ADNIPathsInfo.RData"))

# =============================================================================
#  3) Patient baseline information: one row for each patient
# =============================================================================
aux <- split(demo,demo$PTID)
  # Number of measures for each patient
numObs_perPatient <- sapply(aux,nrow,USE.NAMES = T)
numObs_perPatient <- data.frame(PTID = names(numObs_perPatient),
                                numVisits = numObs_perPatient,
                                row.names = NULL)
  
  # Build dataframe with patient info (one row for each subject)
patientInfo <- lapply(aux,FUN="[",1,-c(2:5))
patientInfo <- do.call(rbind,patientInfo)
colnames(patientInfo)[colnames(patientInfo) == 'EXAMDATE'] <- "enrolDate"

  # Add ADNI cohorts trajectory
#setdiff(patientPaths$PTID,patientInfo$PTID)
patientInfo <- merge(patientInfo,patientPaths,by="PTID")

# Add number of longitudinal measures for each pacient
patientInfo <- merge(patientInfo,numObs_perPatient,by="PTID")



# =============================================================================
# Codebook
# =============================================================================

codebook <- data.frame(Variable = colnames(patientInfo))

codebook$Label <- c("Patient ID",
                     "Date of enrollment",
                     "Diagnostic group",
                     "Age at baseline (years)",
                     "Sex",
                     "Education level (years)",
                     "APOE4",
                     "AT status",
                     paste("PRS",c("AD","AD_noAPOE","FTD","PKSON",
                                   "IEAA","EEAA","TELOM","Frailty",
                                   "LEXP","LONGEVITY","LONGEVITYnoAPOE")),
                    "ADNI project Path",
                    "Number of visits")

codebook$description <- c("Patient identification number",
                          "Date of enrollment",
                          "Patient status: Control (CN), Mild cognitive imperment and Alzheimer (AD)",
                          "Age at baseline",
                          "Reported sex",
                          "Reported number of years of education",
                          "Number of alleles APOE4",
                          "Amyloide-Tau patient status",
                          "Polygenic Risk Score - Neurodegenerative: Alzheimer's disease, with APOE",
                          "Polygenic Risk Score - Neurodegenerative: Alzheimer's disease, without APOE",
                          "Polygenic Risk Score - Neurodegenerative: Frontotemporal dementia",
                          "Polygenic Risk Score - Neurodegenerative: Parkinson disease",
                          "Polygenic Risk Score - Aging (molecular): Intrinsic epigenetic age acceleration",
                          "Polygenic Risk Score - Aging (molecular): Extrinsic epigenetic age acceleration",
                          "Polygenic Risk Score - Aging (molecular): Telomere length",
                          "Polygenic Risk Score - Aging (social): Fraility",
                          "Polygenic Risk Score - Aging (social): Life expectancy",
                          "Polygenic Risk Score - Aging (social): Longevity",
                          "Polygenic Risk Score - Aging (social): Longevity without APOE variants",
                          "ADNI project trajectory",
                          "Total number of visits of any kind")




patientInfo[names(patientInfo)] <- lapply(1:nrow(codebook),function(x){
  
  var <- patientInfo[,codebook[x,"Variable"]]
  attributes(var)$Description <- codebook[x,"description"]
  return(var)
})


# =============================================================================
# SAVE CLEAN DATASETS IN RDATA OBJECT
# =============================================================================

# save clean data ------------------------------------------
save(demo,patientInfo,codebook,
     file=here(cleandataDir,"patientCharacteristics.RData"))
rm(list=ls())


