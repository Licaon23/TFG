###############################################################################
### EXPLORATORY ANALYSIS                                                      #
###   Number of measures distribution by Tesla resolution                     #
###   Version 1                                                               #
###   December 2022                                                           #  
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###                                                                           #
###############################################################################
source("Rscripts/libraries.R")
load(here(cleandataDir,"volumesAndPRS.RData"))
load(here(cleandataDir,"patientCharacteristics.RData"))



# =============================================================================
# H) MEASUREMENT DISTRIBUTION BY TESLA TYPE
# =============================================================================

# FOR EACH PACIENT, GET THEIR PTID, NUMBER OF MRI MEASURES AND 
# WHETHER THESE ARE ONLY IN 1.5 TESLA, 3 TESLA OR A MIX.
dd <- do.call(rbind,data_singleT$data)

ptid_tesla <- lapply(split(dd,dd$PTID),function(x){
  nmeasures <- nrow(x)
  type <- paste(unique(x$tesla), collapse = "/")
  res <- data.frame(PTID=unique(x$PTID),
                    tesla=type,
                    numMRI_measures=nmeasures)
  return(res)
})
ptid_tesla <- do.call(rbind,ptid_tesla)
ptid_tesla$tesla <- as.factor(ptid_tesla$tesla)
levels(ptid_tesla$tesla)[levels(ptid_tesla$tesla) == "3/1.5"] <- "1.5/3"

ptid_tesla <- merge(ptid_tesla,
                    patientInfo[,c("PTID","Path")],
                    by="PTID")



#-----------------------------------------------------------------------------
# PLOT: number of patients with a given ADNI project path,
#       explained by type of measurement: 1.5T, 3T or 1.5/3 T
#-----------------------------------------------------------------------------
aux <- with(ptid_tesla, table(Path , tesla))
end_point <- 0.5 + nrow(aux) + nrow(aux) - 1

png(file = here("results","Others","teslaMeasuresByAdniPhase.png"),
    850,450)

barplot(t(aux),xaxt="n",space=1,ylim=c(0,500),
        col=c("#0066CC","#33CCCC","#CC9933"))
text(seq(1.5, end_point, by = 2), par("usr")[3]-0.25, 
     srt = 35, adj = 1.05, xpd = TRUE, 
     labels = paste(rownames(aux)), cex = 0.65,font=2)
legend(22,500,legend=colnames(aux),
       fill=c("#0066CC","#33CCCC","#CC9933"),
       title="Tesla",cex=1.2)
mtext(side=2,text="Num. patients with a given ADNI path",
      line=2.3,cex=1,font=2)
axis(1,c(0,24),labels = F)
dev.off()

#------------------------------------------------------------------------
# PLOT: Distribution of number of MRI measurements for each patient,
#       both in 1 tesla and 3 teslas.
#-----------------------------------------------------------------------


barplot_MRImeasures <- function(dd){
  # Num.Measures
  aux <- split(dd,dd$PTID)
  numMeasures <- sapply(aux,function(x) nrow(x))
  tab <- table(numMeasures)
  
  # Make BARPLOT
  file <- here("results","Others",paste0("measureDist_byTesla",
                                unique(dd$tesla),".png"))
  png(file,width = 800,height = 480)
  
  barplot(prop.table(tab), las= 1, xlab=" ", ylab="", col="grey",
          cex.lab=1.7, cex.main=1.5, axes= F, ylim=c(0,.27),
          names.arg="")
  
  axis(1,seq(0.7,by=1.2,length.out=length(tab)),labels = names(tab),
       cex.axis=1.5)
  axis(2,seq(0,0.25,by=0.05),cex.axis=1.5,las=1,line=-1.5)
  mtext("Num. MRI measures",side=1, line=2.5, cex=1.5, font=2)
  mtext("Prop. of patients",side=2, line=2.5, cex=1.5, font=2)
  mtext(paste("Measure distribution in ",
              unique(dd$tesla)," teslas"),
        side=3, line=1, cex=1.5, font=2)
  text(x=1.7,y=0.26,labels = paste0("Num.patients: ",sum(tab)),
       cex=1.5)
  
  dev.off()
  
}

lapply(data_singleT$data, barplot_MRImeasures)
rm(list=ls())
