###############################################################################
### IMPORT AND CLEAN SCRIPT                                                   #
###   ADNI participants MRI measurements                                      #
###   Version 1                                                               #
###   November-December 2022                                                  #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:
###     * Import MRI volume data and cleaning tasks                           #
###     * Merge volume data with PRS and demographic                          #
###############################################################################
source("Rscripts/libraries.R")
source(here(RscriptsDir,"1_importAndClean","importAndClean_auxFunctions.R"))
load(here(cleandataDir,"patientCharacteristics.RData")) 
 

# aux <- levels(patientInfo$Path)[1:6]
# patientInfo <- subset(patientInfo,Path %in% aux)
# nrow(patientInfo)
# 
# patientInfo <- patientInfo[!is.na(patientInfo$PRS_AD),]
# nrow(patientInfo)
# 
# with(patientInfo,table(DX.bl,ATgroup))
# exclude <- with(patientInfo,
#                 (DX.bl %in% c("MCI","AD") & ATgroup %in% c("A-T-","A-T+")))
# sum(exclude)
# patientInfo <- patientInfo[!exclude,]
# nrow(patientInfo)


# =============================================================================
# A) IMPORT MRI VOLUMNE DATA and DATA WRANGLING
# =============================================================================

# It's been noticed at first glimpse that missing values appear encoded both as
# NA and "". What is more, it seems that several cases are corrupted,
# with strange values in some fields and even less fields than other cases.

volumes <- read.table(here(rawdataDir, "20211025_HSvolumes.txt"),
                   sep=" ", header=T, fill=T, na.strings = c("",NA))


#----------------------------------------------------------
# Save variable names for hippocampal volume measurements and PRS
#----------------------------------------------------------
aux <- grep("(right|left)",names(volumes))
volumeHippoVars <- names(volumes)[aux]

aux <- grep("PRS",names(patientInfo))
PRSVars <- names(patientInfo)[aux]

#---------------------------------------
# Extract PTID from Subject2
#---------------------------------------
ptid <- volumes$Subject2
regex_PTID <- "^(\\d{3})_(\\w)_(\\d{4})$" # expected format

  # any(!grepl(regex_PTID,ptid)) #is there any PTID with a different format?
  # #we are good;all PTID are consistent format-wise.

#-----------------------------------------------------------
# Replace PTID current value by only its last 4-digit code 
#-----------------------------------------------------------
volumes$PTID <- as.factor(gsub(regex_PTID,"\\3",volumes$Subject2))


#---------------------------------------
# Extract Measure Date  from Subject
#---------------------------------------
regex_date <- "^(.*)_(\\d{4}-\\d{2}-\\d{2}).*$"
volumes$measureDate <- gsub(regex_date,"\\2",volumes$Subject)
volumes$measureDate <- as.Date(volumes$measureDate)

#------------------------------------------------------------------------
# From patientInfo dataset, join demograhpic characteristics and PRS 
#------------------------------------------------------------------------
volumes <- merge(volumes,
                 patientInfo[,c("PTID","enrolDate","AGE","PTGENDER","PTEDUCAT",
                                "DX.bl","ATgroup","APOE4","Path",PRSVars)], 
                 by="PTID")

#-------------------------------
# Calculate timelapse variable 
#-------------------------------
  # For each patient, add their enrollment Date and calculate the time 
  # difference between the measure date and the enrollment date.
volumes$timelapse <- with(volumes, measureDate-enrolDate)

  # Timelapse in months...
volumes$timelapse <- as.numeric(volumes$timelapse)/30

#----------------------------------------------------
# Calculate VISCODE2 variable according to timelapse
#---------------------------------------------------
volumes$VISCODE2 <- get_viscode(volumes$timelapse)


#------------------------------------------------------------------------
# Reorder variables in a more sensible order, and get rid of those that
# now are redundant.
#------------------------------------------------------------------------

# Get rid of Subject and Subject2 variables
volumes[c("Subject","Subject2")] <- NULL

# Change names for EstimatedTotalIntraCranialVol
aux <- which(colnames(volumes)=="EstimatedTotalIntraCranialVol")
colnames(volumes)[aux] <- "ICV"

# Reorder variables
volumes <- volumes[,c("PTID","AGE","DX.bl","PTGENDER","PTEDUCAT","ATgroup","APOE4",
                      PRSVars,"measureDate","enrolDate","timelapse","VISCODE2",
                      "tesla","Path","ICV",volumeHippoVars)]

# Arrange observations in a ascending order of PTID and measure date --------
volumes <- volumes[order(volumes$PTID,volumes$measureDate),]
rownames(volumes) <- NULL

#volumes[1:6,1:8]
# Delete objects from Glob.Env with no further utility
rm(list=c("ptid","regex_PTID","regex_date"))

# =============================================================================
# B) MISSING VALUES: volumes measurement
# =============================================================================

#------------------------------------------------------------------------
# 1) Calculate percentage of complete cases, proportion of missing values
#    by variable and by patient
#------------------------------------------------------------------------

# i) PERCENTATGE OF COMPLETE CASES 
  (perc_cc <- mean(complete.cases(volumes[volumeHippoVars]))*100)

# ii) PROPORTION OF MISSINGS BY VOLUME VARIABLES
  (NAsByVar <- sapply(volumes[volumeHippoVars],function(x) mean(is.na(x))*100))
    # since all variables have the same percentage of missings,
    # it is likely that incomplete cases are complete empty cases.

# iii) PROPORTION OF MISSING VOLUMES BY PATIENT
  aux <-  sapply(split(volumes[volumeHippoVars],volumes$PTID),
               FUN=function(x) mean(is.na(x))*100)
  (ptid_withNA <- sort(aux[which(aux>0)]))
  
  
  # There are 20 patients total with some missing measurements.
  # 10 of them have no valid volume measurement whatsoever; therefore,
  # all their cases will be dropped from dataset.
  # A part from that, there are another 10 patients with some degree
  # of missingness. These patients will be considered in a separated 
  # dataframe, and removed from the main one.

#--------------------------------------------------------------------
# 2) Remove all patients with any measure at all. aka: all cases are NA
#    values for all visits
#--------------------------------------------------------------------
  aux <- names(which(ptid_withNA==100)) #10 ptid and 37 cases
  volumes <- volumes[!(volumes$PTID %in% aux),]
  
  # After removing them, there will be 1395 patients i 6954 cases left

#------------------------------------------------------------------------
# 3) Identify patients whose measurement record is incomplete,
#    save them in a separate dataframe and remove them from main
#    dataset.
#-------------------------------------------------------------------------
  aux <- names(which(ptid_withNA!=100)) #ptid with some missing values
  rareCases <- volumes[volumes$PTID %in% aux,] #10 ptid, 73 cases
  rareCases <- rareCases[with(rareCases,order(PTID,tesla,measureDate)),]
  volumes <- volumes[!(volumes$PTID %in% aux),]

# Remove objects from the global environment that are of no use anymore,
# and drop empty ptid levels
  volumes$PTID <- droplevels(volumes$PTID)
  rm(list=c("aux"))
  
  # length(unique(volumes$PTID))
  # nrow(volumes)
  # Final: 6881 cases for 1385 different patients

  
# =============================================================================
# C) MISSING VALUES: PRS scores
# =============================================================================
  
  #------------------------------------------------------------------------
  # 1) Calculate percentage of complete cases, proportion of missing values
  #    by variable and by patient
  #------------------------------------------------------------------------
  
  # i) PERCENTATGE OF COMPLETE CASES 
  (perc_cc <- mean(complete.cases(volumes[PRSVars]))*100)
  
  # ii) PROPORTION OF MISSINGS BY VOLUME VARIABLES
  (NAsByVar <- sapply(volumes[PRSVars],function(x) mean(is.na(x))*100))
  # since all variables have the same percentage of missings,
  # it is likely some patients does not have any PRS value
  
  # iii) PROPORTION OF MISSING VOLUMES BY PATIENT
  aux <-  sapply(split(volumes[PRSVars],volumes$PTID),
                 FUN=function(x) mean(is.na(x))*100)
  (ptid_withNA <- names(sort(aux[which(aux>0)])))
  
    length(ptid_withNA) # 125 patients does not have PRS information

  #------------------------------------------------------------------------
  # 2) Check whether missing distribution is independent of diagnostic or 
  #     sex group.
  #------------------------------------------------------------------------
  # Dades de característqiues basals ----------------------------------------
  basal <- lapply(split(volumes,volumes$PTID),FUN="[",1,)
  basal <- do.call(rbind,basal)
  #head(basal)
  
  # Missing values distribution by DX.bl ------------
  aux <- data.frame(Grup = basal$DX.bl, 
                    Missing = is.na(basal$PRS_AD))
  tab <- table(aux$Grup,aux$Missing)
  chisq.test(tab)
  
    # There is no significant difference in missing proportion between 
    # diagnostic group, regarding PRS.
  
  
  # Missing values distribution by PTGENDER ------------
  aux <- data.frame(Grup = basal$PTGENDER, 
                    Missing = is.na(basal$PRS_AD))
  tab <- table(aux$Grup,aux$Missing)
  chisq.test(tab)

    # !!Proportion of missing values of PRS is significantly higher in female
    # than in male...!!
  
  #------------------------------------------------------------------------
  # 3) Check whether complete cases analysis modify significantly mean values
  #     for age and volumes
  #------------------------------------------------------------------------
  
  testMissingEfect <- function(numVar,toExclude){
    aux <- data.frame(numVar = basal[[numVar]], 
                      Missing = toExclude)
    mean_var_prev <- mean(aux$numVar)
    median_var_prev <- median(aux$numVar)
    
    test <- t.test(aux$numVar[aux$Missing], mu = mean_var_prev, data=aux)
    #wilcox.test(aux$numVar[aux$Missing], mu = median_var_prev, data=aux)
    
    res <- data.frame(Variable = numVar,
                      pvalue = round(test$p.value,3))
    
    return(res)
  }
  
  aux <- lapply(names(basal)[c(2,25:51)], FUN=testMissingEfect,
                toExclude = is.na(basal$PRS_AD))
  missingEffect <- do.call(rbind,aux)
  missingEffect
  
  # No hi ha efecte sobre els valors mitjans de les numèriques al treballar
  # amb casos complets.
  # A més, els NA en els PRS estan distribuits at random per grup diagnòstic,
  # però no per sexe (més missing en dones). BIAIX??
  
### THEREFORE... we will analysis complete cases in terms of PRS  
volumes <- volumes[!is.na(volumes$PRS_AD),]
volumes$PTID <- droplevels(volumes$PTID)
  # 6288 observations for 1260 individuals.
  length(unique(volumes$PTID))



# =============================================================================
# D) EXCLUSION CRITERIA: exclude AD and MCI subjects with Amiloyde negative 
#     status
# =============================================================================
# Dades de característqiues basals ----------------------------------------
basal <- lapply(split(volumes,volumes$PTID),FUN="[",1,)
basal <- do.call(rbind,basal)
with(basal,table(DX.bl,ATgroup,useNA = "ifany"))  


# exclude1 <-with(basal,
#   (DX.bl %in% c("MCI","AD") & ATgroup %in% c("A-T-","A-T+")) |
#   (DX.bl == "CN" & basal$ATgroup == "A-T+") |
#   is.na(ATgroup))

exclude2 <-with(basal,
  (DX.bl %in% c("MCI","AD") & ATgroup %in% c("A-T-","A-T+")))

#  sum(exclude1) # 588 patients to be excluded
  sum(exclude2) # 157 patients to be excluded
  
  # Differences in PTGENDER? .........................................
  chisq.test(table(exclude2,basal$PTGENDER))
  # No

  # Differences in numeric variables? .................................
  aux <- lapply(names(basal)[c(2,25:51)], FUN=testMissingEfect,
                toExclude = exclude2)
  missingEffect <- do.call(rbind,aux)
  missingEffect
    # !Excluding AD cases with AD A- status changes significantly the mean
    # value for several volume measurements.
    # Nevertheless, it is necessary to exlude them.
  
# EXCLUDE CASES -------------------------------------
  # exclude <- with(volumes,
  #   (DX.bl %in% c("MCI","AD") & ATgroup %in% c("A-T-","A-T+")) |
  #   (DX.bl == "CN" & ATgroup == "A-T+") |
  #   is.na(ATgroup))

  exclude <-with(volumes,
    (DX.bl %in% c("AD","MCI") & ATgroup %in% c("A-T-","A-T+")))
  
volumes <- volumes[!exclude,]
volumes$PTID <- droplevels(volumes$PTID)

   length(unique(volumes$PTID));dim(volumes)
  # 1103 individuals, 5433 observations



# =============================================================================
# E) MEASUREMENTS BY TESLA PARAMETER
# =============================================================================
  
#-------------------------------------------------
# PTID with repeated measures, both in 1.5 and 
# 3 teslas for similar dates
#-------------------------------------------------
  aux <- lapply(split(volumes,volumes$tesla),function(x)unique(x$PTID))
  PTID_withBothTeslaMeasure <- Reduce(intersect,aux)
  
  length(PTID_withBothTeslaMeasure)
    # among 1103 PTID with complete data available, 149 of them have been 
    # measured in both resolutions,whereas 954 subjects left have only 1 
    #type-resolution measurements.
  
  ### Subset
  volumes_bothT <- volumes[volumes$PTID %in% PTID_withBothTeslaMeasure,] 
  volumes_bothT$PTID <- droplevels(volumes_bothT$PTID)
  
  length(unique(volumes_bothT$PTID));dim(volumes_bothT)
  #1171 observations, 149 patients
  
  ### Order by tesla and measureDate
  aux <- with(volumes_bothT, order(PTID,tesla,measureDate))
  volumes_bothT <- volumes_bothT[aux,]
  
  ### Identify Patients with two series of measurements in different resolutions 
  ### and patients with one sequence during which resolution was updated.
  
  aux <- split(volumes_bothT,volumes_bothT$PTID)
  dups <- sapply(aux,function(x)any(table(x$VISCODE2)>1))
  
  # One sequence. Add patients missclassified by visual inspection
  volumes_oneSeq <- c(aux[!dups],aux[c("0089","0996")]) #el "0996" ja no hi és
  volumes_oneSeq <- do.call(rbind,volumes_oneSeq)
  rownames(volumes_oneSeq) <- NULL
    
  length(unique(volumes_oneSeq$PTID))
  # 32 individuals
  
  # Two sequences
  volumes_twoSeq <- aux[dups]
  volumes_twoSeq[c("0089","0996")] <- NULL # el "0996" ja no hi es
  volumes_twoSeq <- do.call(rbind,volumes_twoSeq)
  rownames(volumes_twoSeq) <- NULL
  # 117 individuals
  ptid_twoSeq <- as.character(unique(volumes_twoSeq$PTID))
  length(ptid_twoSeq)
  
  volumes_twoSeq <- split(volumes_twoSeq,volumes_twoSeq$tesla)
  
  
  
  
#-------------------------------------------------
# PTID with a single tesla resolution type
#-------------------------------------------------
  aux <- unique(volumes_oneSeq$PTID)
  
  volumes <- volumes[!(volumes$PTID %in% aux),]
  #volumes_singleT$PTID <- droplevels(volumes_singleT$PTID)
  
  length(unique(volumes$PTID));dim(volumes)
  # 5168 cases de 1071 different patients. 

  volumes_singleT <- split(volumes,volumes$tesla)
  
  

  volumes_singleT <- lapply(volumes_singleT, function(x){
                            x$PTID <- droplevels(x$PTID)
                            x <- x[order(x$PTID,x$measureDate),]
                            rownames(x) <- NULL
                            return(x)})

lapply(volumes_singleT,function(x)c(length(unique(x$PTID)),nrow(x)))
# 1.5T: 580 pacients, amb 2567 obs.
# 3T: 608 pacients, am 2601 obs. 


#Comprovació dels merge
# a <- summaryBy(cbind(AGE,DX.bl,PTGENDER,ATgroup)~PTID,data=volumes_singleT$`1.5`,
#           FUN=function(x)length(unique(x)))
# summary(a)


# =============================================================================
# D) CLEANING TASKS
# =============================================================================

# In order to clean up some inconsistencies with data, several procedures
# have been implemented in the main function 'cleanTasks', which makes use
# of several auxiliary functions (see auxFunctions.R) that performs the
# following tasks:

# A) No baseline cases----------------------------------
#
# i) Identify patients with no bl measure, aka, first measure they have was
#     taken at least three months after baseline visit (enrolmentDate)
# ii) For each one of these cases, first measurement date is considered as
#     enrolment date. Consequently, baseline AGE is also modified in order to
#     match their actual age in the first measurement date. These tasks are
#     performed by 'fix_nobl' auxiliary function.
# iii) Timelapse is then recalculated, and VISCODE as well (by using
#     'get_viscode' auxiliary function).
# iv) Finally, previous values are replaced by these updates.

# B) Repeated measurements-------------------------
#   v) Some patients had two measurements taken very close in time, meaning
#      its VISCODE is the same. We assumed that there was some sort of error
#      during that scan session, and therefore patient was asked to repeat the
#      session. Thus, in case of double VISCODE (two measurements taken close one
#      to another), only the latest will be kept as the corrected.
#      This tasks is implemented and performed by 'get_repeatedViscode' auxiliary
#      function.

#   vi) Once these repeated measurements are identified, the first one in time
#       is dropped. This task is implemented and performed by the auxiliay
#       function 'dropCases'.

cleanTasks <- function(dd){
  
  # A) No baseline cases ###################################################
  
  # STEP i)
  nobl <- names(which(sapply(split(dd,dd$PTID),
                             FUN=function(x) table(x$VISCODE2)[["bl"]]==0)))
  
  position <- which(dd$PTID %in% nobl)
  data_nobl <- dd[position,]
  
  # STEP ii) 
  data_nobl <- fix_nobl(data_nobl)
  
  # STEP iii)
    # Timelapse in months...
  data_nobl$timelapse <- with(data_nobl, measureDate-enrolDate)
  data_nobl$timelapse <- as.numeric(data_nobl$timelapse)/30
  
    #VISCODE
  data_nobl$VISCODE2 <- get_viscode(data_nobl$timelapse)
  
  # STEP iv)
  dd[position,names(data_nobl)] <- data_nobl
  
  
  #B) Repeated viscodes ##################################################
  
  # STEP v)
  labels <- levels(dd$VISCODE2)
  aux <- lapply(labels,get_repeatedViscode,dd=dd)
  aux <- aux[sapply(aux,length)!=0]
  
  # STEP vi)
  drop <- lapply(aux,dropCases,dd=dd)  
  drop <- Reduce('c',drop)
  dd <- dd[-drop,]  
  
  ### RETURN #############################################################
  return(dd)
}


### Apply cleanTasks function over volumes_singleT
volumes_singleT <- lapply(volumes_singleT, cleanTasks)

lapply(volumes_singleT,function(x)c(length(unique(x$PTID)),nrow(x)))
# 1.5T: 580 patients, and 2556 obs. (38 nobl and 11 patients duplicated viscode)
# 3T: 608 patients, and 2577 obs. (44 nobl and 24 patients duplicated viscode)

#CHECKS
# a <- summaryBy(cbind(AGE,DX.bl,PTGENDER,ATgroup)~PTID,data=volumes_singleT$`1.5`,
#                FUN=function(x)length(unique(x)))
# summary(a)
# labels <- levels(demo$VISCODE2)
# dd <- volumes_singleT$`1.5`
# sapply(levels(demo$VISCODE2),f,dd=dd) #no hi ha etiquetes repetides
# names(which(sapply(split(dd,dd$PTID),
#                    FUN=function(x) table(x$VISCODE2)[["bl"]]==0))) #0 no-bl


# C) Drop patients with only baseline measurement------------------------------

volumes_singleT <- lapply(volumes_singleT, function(x){
  nmeasures <- sapply(split(x,x$PTID),nrow)
  drop_ptid <- names(which(nmeasures == 1))
  res <- x[!(x$PTID %in% drop_ptid),]
  res$PTID <- droplevels(res$PTID)
  
  return(res)
  })

lapply(volumes_singleT,function(x)c(length(unique(x$PTID)),nrow(x)))
#1.5T : 18 patients with only one mesurement, 562 patients, 2538 obs.
#3T : 3 patients with only one mesurement, 605 patients, 2574 obs.



# =============================================================================
# E) CREATING A CONTINOUS LONGITUDINAL TIME MEASUREMENT 
# =============================================================================

# A part from VISCODE variable, which is a discrete time longitudinal
# measurement, we can create a continuous version by using patient's baseline
# age and adding up timelapse variable (aka, time between measurement and
# baseline visit)

volumes_singleT <- lapply(volumes_singleT, function(x){
  # Patient age every time they underwnet a MRI session
  x$time_age <- round(x$AGE + x$timelapse/12,1)
  
  #Numeric variable corresponding to months after baseline visit
  x$time_visit <- sub("(m)(\\d+)","\\2",x$VISCODE2)
  x$time_visit <- ifelse(x$time_visit=="bl",0,x$time_visit)
  x$time_visit <- as.integer(x$time_visit)
  x <- x[,c("PTID","AGE","DX.bl","PTGENDER","PTEDUCAT","ATgroup","APOE4",
            "enrolDate","Path",PRSVars,
            "tesla","measureDate","VISCODE2","timelapse",
            "time_age","time_visit","ICV",
            volumeHippoVars)]
  rownames(x) <- NULL
  return(x)
})




# =============================================================================
# F) CREATING OVERALL HIPPOCAMPAL MEASUREMENTS
# =============================================================================

# For each hippocampal subfield, create a new variable adding up the volume
# from both left and right hemispheres

# i) Get hippocampal field names from variable names
hippocampal_subfields <- unique(gsub("^(left|right)_(*.)",
                                     "\\2",
                                     volumeHippoVars))



# ii) Apply get_overallVolume function to all pairs of measurements
#     and add each result in a new variable.
volumes_singleT_overall <- volumes_singleT
#1.5T
volumes_singleT_overall$`1.5`[hippocampal_subfields] <- 
  lapply(hippocampal_subfields, get_overallVolume, dd=volumes_singleT$`1.5`)

#3T
volumes_singleT_overall$`3`[hippocampal_subfields] <- 
  lapply(hippocampal_subfields, get_overallVolume,dd=volumes_singleT$`3`)

# iii) Keep only overall volumes for each hippocampal field.
#
volumes_singleT_overall <- lapply(volumes_singleT_overall,function(x){
  x[volumeHippoVars] <- NULL
  return(x)})

### CHECK
# dd1 <- volumes_singleT$`1.5`
# dd2 <- volumes_singleT_overall$`1.5`
# measure <- sample(hippocampal_subfields,1)
# ptid <- sample(unique(dd1$PTID),1)
# 
# check <- subset(dd1,PTID==ptid,
#               select=c(1:4,21:23,grep(pattern = paste0("(left|right)_",measure),
#                                         names(dd1))))
# check <- merge(check,dd2[,c("PTID","VISCODE2",measure)],by=c("PTID","VISCODE2"))
# check$check <- identical(check[,8]+check[,9],check[,10])
# check

# =============================================================================
# G) Changing list names and saving clean datasets and adding attributes
# =============================================================================

#Adding "T" to the tesla resolution ------------------------------
names(volumes_singleT) <- paste0("T_",names(volumes_singleT))
names(volumes_singleT_overall) <- paste0("T_",names(volumes_singleT_overall))


# =============================================================================
# H) ADD A JOINT DATASET WITH INDIVIDUALS WITH EITHER 1.5T AND 3T SCANS.
#    FOR THOSE WITH TWO SEQUENCES, ADD THE 3T ONE
# =============================================================================
volumes_mixed <- rbind(volumes_singleT_overall$T_3,
                       subset(volumes_singleT_overall$T_1.5,
                              !(PTID %in% ptid_twoSeq)))

volumes_mixed <- volumes_mixed[order(volumes_mixed$PTID),]
#length(unique(volumes_mixed$PTID)) # 1053 individuals 


#---------------------------------------------------------------
# i) CODEBOOK: give description to each variable in datasets
#---------------------------------------------------------------
### Codebook

codebookBoth <- data.frame(Variable = colnames(volumes_singleT$T_1.5))
codebookBoth$Label <- c("Patient ID",
                     "Age at baseline (years)",
                     "Diagnostic group",
                     "sex",
                     "Education level (years)",
                     "AT status",
                     "APOE4 carriership",
                     "Date of enrollment",
                     "ADNI phases path",
                     paste("PRS",c("AD","AD_noAPOE","FTD","PKSON",
                                   "IEAA","EEAA","TELOM","Frailty",
                                   "LEXP","LONGEVITY","LONGEVITYnoAPOE")),
                     "MRI scan resolution (teslas)",
                     "Date of measurement",
                     "Visit code",
                     "Time lapse",
                     "time_age",
                     "time_visit",
                     "Intracranial Volume (mm3)",
                     paste("Left",hippocampal_subfields,"(mm3)"),
                     paste("Right",hippocampal_subfields,"(mm3)"))


codebookBoth$description <- c("Patient identification number",
                           "Age at baseline",
                           "Patient status: Control (CN), Mild cognitive imperment and Alzheimer (AD)",
                           "Reported biological sex",
                           "Reported number of years of education",
                           "Amyloide-Tau patient status",
                           "Number of alleles APOE4 (0,1,2)",
                           "Date of enrollment: when the first baseline visit took place",
                           "Different phases for which subject has gone through",
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
                           "Polygenic Risk Score - Aging (social): Longevity withou APOE variants",
                           "MRI resolution in teslas: either 1.5T or 3T",
                           "Date of MRI measurement",
                           "Visit code, coding number of months after baseline visit",
                           "Difference between measureDate and enrollDate, in months",
                           "Longitudinal time measure: current patient age, when a given MRI measurement was taken",
                           "Longitudinal time measure: time in months after baseline visit",
                           "Estimated Total Intra-cranial volumes, in cubic millimeters",
                           rep("Hippocampal subfields volume measurements, in cubic millimeters",26))

#Codebook for overall volumes datasets
codebookOverall <- data.frame(Variable = colnames(volumes_singleT_overall$T_1.5))
codebookOverall$Label <-c("Patient ID",
                          "Age at baseline (years)",
                          "Diagnostic group",
                          "Sex",
                          "Education level (years)",
                          "AT status",
                          "APOE4 carriership",
                          "Date of enrollment",
                          "ADNI phases path",
                          paste("PRS",c("AD","AD_noAPOE","FTD","PKSON",
                                        "IEAA","EEAA","TELOM","Frailty",
                                        "LEXP","LONGEVITY","LONGEVITYnoAPOE")),
                          "MRI scan resolution (teslas)",
                          "Date of measurement",
                          "Visit code",
                          "Time lapse",
                          "time_age",
                          "time_visit",
                          "Intracranial volume (mm3)",
                          paste(hippocampal_subfields,"(mm3)"))


codebookOverall$description <- c("Patient identification number",
                           "Age at baseline",
                           "Patient status: Control (CN), Mild cognitive imperment and Alzheimer (AD)",
                           "Reported biological sex",
                           "Reported number of years of education",
                           "Amyloide-Tau patient status",
                           "Number of alleles APOE4 (0,1,2)",
                           "Date of enrollment: when the first baseline visit took place",
                           "Different phases for which subject has gone through",
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
                           "MRI resolution in teslas: either 1.5T or 3T",
                           "Date of MRI measurement",
                           "Visit code, coding number of months after baseline visit",
                           "Difference between measureDate and enrollDate, in months",
                           "Longitudinal time measure: current patient age, when a given MRI measurement was taken",
                           "Longitudinal time measure: time in months after baseline visit",
                           "Estimated Total Intra-cranial volumes, in cubic millimeters",
                           rep("Hippocampal subfields volume measurements, in cubic millimeters",13))

# Set variable description as attributes-------------------------------------

volumes_singleT$T_1.5[codebookBoth$Variable] <-lapply(1:nrow(codebookBoth),
                                                  FUN=set_attributes,
                                                  dataset=volumes_singleT$T_1.5,
                                                  codebook=codebookBoth)

volumes_singleT$T_3[codebookBoth$Variable] <-lapply(1:nrow(codebookBoth),
                                                FUN=set_attributes,
                                                dataset=volumes_singleT$T_3,
                                                codebook=codebookBoth)

volumes_singleT_overall$T_1.5[codebookOverall$Variable] <-lapply(1:nrow(codebookOverall),
                                                        FUN=set_attributes,
                                                        dataset=volumes_singleT_overall$T_1.5,
                                                        codebook=codebookOverall)

volumes_singleT_overall$T_3[codebookOverall$Variable] <-lapply(1:nrow(codebookOverall),
                                                       FUN=set_attributes,
                                                       dataset=volumes_singleT_overall$T_3,
                                                       codebook=codebookOverall)


volumes_mixed[codebookOverall$Variable] <-lapply(1:nrow(codebookOverall),
                                                        FUN=set_attributes,
                                                        dataset=volumes_mixed,
                                                        codebook=codebookOverall)




# =============================================================================
# G) Save clean datasets
# =============================================================================
# Objects in lists
data_singleT <- list(data = volumes_singleT, 
                     codebook = codebookBoth,
                     Info = "subjects who have a sequence of MRI measurements taken with the same tesla parameter. Volumes are available for both hemispheres")
                     
                     
data_singleT_overall <- list(data=volumes_singleT_overall,
                             codebook = codebookOverall,
                             Info = "subjects who have a sequence of MRI measurements taken with the same tesla parameter. Hippocampal subfields volumes are reported as the sum of both hemispheres")

data_bothT <- list(data = volumes_oneSeq, 
                   Info = "35 subjects whose MRI scans changed tesla parameter at one point of their sequence")

data_withNA <- list(data = rareCases,
                    Info = "These subjects were excluded because had missing values for volume measurements. Nevertheless, some of them have two sequences of measurements, one of which do have complete observations. Thus, they could be merged with the regular dataset")

data_singleT_mixed <- list(data = volumes_mixed,
                           codebook = codebookOverall,
                           Info = "This dataset contains volume data from subjectes with a sequence of scans either in 1.5 tesla or 3 tesla. Merged from data in data_singleT")
#Save
save(data_singleT , data_singleT_overall, data_bothT,
     data_withNA, data_singleT_mixed,
     file=here(cleandataDir,"volumesAndPRS.RData"))
rm(list=ls())
