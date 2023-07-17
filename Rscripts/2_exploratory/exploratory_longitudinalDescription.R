###############################################################################
### EXPLORATORY DATA ANALYSIS                                                 #
###   Longitudinal measurements                                               #
###   Version 1                                                               #
###   March 2023                                                            #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###                                                                           #
###############################################################################
source("Rscripts/libraries.R")
#source(here(RscriptsDir,"2_exploratory","descSummary_bis.R")) # aux function

load(here(cleandataDir,"dataForAnalysis_3T.RData"))
paramT <- "Dataset_3T"

# load(here(cleandataDir,"dataForAnalysis_1T.RData"))
# paramT <- "Dataset_1T"

# =============================================================================
# 1) DESCRIPTIVA MESURES LONGITUDINALS
# =============================================================================
volVar <- sort(names(dd)[21:33]) # variables volums

#----------------------------------------------------------------------
# i) Nombre observacions per cada temps, segons grup
#---------------------------------------------------------------------
aux <- summaryBy(Hippocampal_tail ~ DX.bl+VISCODE2,
                              FUN=length,data=dd)
names(aux)[3] <- "numMeasures"

numMeasures_wide <- reshape(aux, direction = "wide",
                            v.names = "numMeasures",
                            idvar = "VISCODE2",
                            timevar = "DX.bl")

names(numMeasures_wide) <- c("VISCODE2","CN","MCI","AD")
numMeasures_wide[is.na(numMeasures_wide)] <- 0

# To proportions
numMeasures_wide_prop <- numMeasures_wide

numMeasures_wide_prop[c("CN","MCI","AD")]<-
  lapply(numMeasures_wide_prop[c("CN","MCI","AD")],
         FUN=function(x)round(x/x[1],3))

# Long format
numMeasures_long_prop <- reshape(numMeasures_wide_prop, direction="long",
                                 varying=c("CN","MCI","AD"),
                                 v.names="propMeasures", idvar="VISCODE2",
                                 times=c("CN","MCI","AD"),timevar = "DX.bl")
numMeasures_long_prop$DX.bl <- factor(numMeasures_long_prop$DX.bl,
                                      levels = c("CN","MCI","AD"))

numMeasures_long_prop <- merge(numMeasures_long_prop,aux,
                               by=c("DX.bl","VISCODE2"),all.x=T)
#numMeasures_long_prop[is.na(numMeasures_long_prop)] <- 0

#----------------------------------------------------------------------
# ii) PLOT
#---------------------------------------------------------------------

ggplot(numMeasures_long_prop,aes(x=VISCODE2, y=propMeasures,group=DX.bl)) +
  geom_vline(xintercept = "m60",linetype="dashed",col="red4") +
  geom_hline(yintercept = 0.1,linetype="dotted",col="red4")+
  geom_point(aes(color=DX.bl),show.legend = F,size=2) + 
  geom_line(aes(color=DX.bl),show.legend = F,size=.5) +
  geom_text_repel(aes(label=numMeasures,color=DX.bl),
                  box.padding = unit(0.35,"lines"),
                  point.padding = unit(0.3,"lines"),
                  show.legend = F) +
  
  ggstatsplot::theme_ggstatsplot() +
  scale_color_manual(values = c("#1B9E77","#D95F02","#7570B3")) + 
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank()) +
  ggtitle(paste0("Proportion of measurements at different follow-up visits by diagnostic group")) +
  labs(subtitle = "<span style = 'color: #1B9E77;'>**Controls**</span>, <span style = 'color: #D95F02;'>**Mild cognitive impairment**</span> and <span style = 'color: #7570B3;'>**Alzheimer`s disease**</span>" ,
       caption="**Source**: ADNI project | **Units** Proportion of available measurements at each visit with respect to total number of pacients per group. Time after baseline visit in months <br/> Annotations are number of patients who underwent MRI session at a given visit after baseline.") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = -.05))

ggsave(filename = "propMeasuresbyTime.png",
       path = here("results",paramT,"1_DescriptiveAnalysis",
                   "2_longitudinalMeasurements"))


#----------------------------------------------------------------------
# iii) TABLE
#---------------------------------------------------------------------
tab <- t(numMeasures_wide[-1])
colnames(tab) <- numMeasures_wide$VISCODE2

tab_measures <- kbl(tab,full_width=F,font="Cambria",
                    caption="Number of patients with measurement at each time")
tab_measures <- kable_classic(tab_measures,full_width=F)

#tab_measures <- column_spec(tab_measures,width="2cm")
tab_measures <- row_spec(tab_measures,0,bold=T)


file <- "table_numMeasuresByTime_bis.png"
path <- here("results",paramT,"1_DescriptiveAnalysis",
             "2_longitudinalMeasurements")

if(!(any(grepl(file,list.files(path))))){
  save_kable(tab_measures,file=here(path,file))
}


# =============================================================================
# 2) Lineplots de trajectòries longintudinals per pacients
# =============================================================================

#----------------------------------------------------------------------
# i) Funció
#---------------------------------------------------------------------
volume_trajectory <- function(volume,data,plot=TRUE){
  
  dd <- data[,c("PTID","DX.bl","time_age",volume)]
  names(dd)[4] <- "volume"
  
  
  ptid_sel <- sample(unique(dd$PTID),size = 25)
  highlighted <- subset(dd,PTID %in% ptid_sel)
  
  ggplot(data=dd,aes(time_age,y=volume,group=PTID)) + 
    geom_line(aes(color=DX.bl),alpha=0.3,show.legend = F)+
    
    geom_point(data=highlighted,
               aes(color=DX.bl),alpha=0.8,size=1,show.legend = F) +
    geom_line(data=highlighted,
              aes(color=DX.bl),alpha=0.8,linewidth=.5,show.legend = F)+
    
    facet_wrap(~DX.bl)+
    
    ggstatsplot::theme_ggstatsplot() +
    scale_color_manual(values = c("#1B9E77","#D95F02","#7570B3"))+
    theme(#strip.background = element_blank(),
      plot.margin = margin(1,1,.5,1,"cm"),
      axis.title = element_blank()) +
    ggtitle(paste("Individual patient evolution for",volume,"volume")) +
    labs(subtitle = "<span style = 'color: #1B9E77;'>**Controls**</span>, <span style = 'color: #D95F02;'>**Mild cognitive impairment**</span> and <span style = 'color: #7570B3;'>**Alzheimer`s disease**</span>" ,
         caption="**Source**: ADNI project | **Diagnostic groups** CN: controls, AD: Alzheimer`s disease, MRI: mild cognitive impairment |\n**Units** Volume in mm3. Patient age at each visit") +
    theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
          plot.subtitle = element_markdown(size = 12, hjust = 0),
          plot.title = element_text(size = 16, hjust = -.05)) -> p
  
  if(plot==T){
    cat("Plotting",volume,"...")
    ggsave(filename = paste0("individualPaths_",volume,".png"),
           path = here("results",paramT,"1_DescriptiveAnalysis",
                       "2_longitudinalMeasurements",
                       "patient_trajectories"),
           width=30,height=15,units="cm")
  }
  else{cat("Plotting patient trajectories for",volume,"...");p}
}
  

#----------------------------------------------------------------------
# ii) Plots sobre els 13 volums
#---------------------------------------------------------------------
lapply(volVar,volume_trajectory,data=dd)
volume_trajectory(volVar[5],data=dd,plot=F)



# =============================================================================
# 3) Gràfics de mitjanes i intervals de confiança segons edat
# =============================================================================

#----------------------------------------------------------------------
# ii) Funció
#---------------------------------------------------------------------
meanPlot_long <- function(volume,data,plot=T,window=0.5,minN = 3){
  
  # Prepare data-----------------------------------------------
  thres <- seq(54,95,by=window) # thresholds
  b <- cut(data$time_age,breaks=thres) #age groups

    # Edats de referència pels grups d'edat
  referenceAge <- data.frame(timeGroup =levels(b),
                             referenceAge = diff(thres)/2 + thres[1:length(thres)-1])
  data$timeGroup <- b
  data <- merge(data,referenceAge,by="timeGroup")

    # Seleccionar volum d'interès i calcular mitjana i CI per cada grup d'edat
  data <- data[,c("PTID","DX.bl","timeGroup","referenceAge",volume)]
  names(data)[5] <- "volume"
  
  summary <- summaryBy(volume~DX.bl+timeGroup,data=data,
                       FUN=c(length,smean.cl.normal),
                       id = "referenceAge",
                       fun.names = c("N","Mean","Lower","Upper"))
  summary <- subset(summary,volume.N > minN)
  
  # Plot lineplots--------------------------------------------------------
  ggplot(summary,aes(x=referenceAge,y=volume.Mean,group=DX.bl,
                     color = DX.bl)) +
    geom_point(size=.5,show.legend = F) +
    geom_line(show.legend = F,linewidth=.5) +
    geom_line(aes(y=volume.Lower),linetype="dashed",
              linewidth=.2,show.legend = F) +
    geom_line(aes(y=volume.Upper),linetype="dashed",
              linewidth=.2,show.legend = F) +
    
    # geom_pointrange(aes(y=volume.Mean,ymin=volume.Lower,ymax=volume.Upper,
    #                     color=DX.bl),show.legend = F,size=.2) +
    facet_wrap(~DX.bl) +
    
    ggstatsplot::theme_ggstatsplot() +
    scale_color_manual(values = c("#1B9E77","#D95F02","#7570B3")) + 
    theme(strip.background = element_blank(),
          plot.margin = margin(1,1,.5,1,"cm"),
          axis.title = element_blank()) +
    ggtitle(paste0("Longitudinal evolution of volum in ",volume)) +
    labs(subtitle = "Mean value and 95% confidence bands along patient age.<br/><span style = 'color: #1B9E77;'>**Controls**</span>, <span style = 'color: #D95F02;'>**Mild cognitive impairment**</span> and <span style = 'color: #7570B3;'>**Alzheimer`s disease**</span>" ,
         caption="**Source**: ADNI project | **Diagnostic groups** CN: controls, AD: Alzheimer`s disease, MRI: mild cognitive impairment | **Units** Mean volume in mm3. Patient age in years <br/>
       **Age groups**: patient ages along follow up have been categorized in equally size intervals, for which mean volume has been calculated.") +
    theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
          plot.subtitle = element_markdown(size = 12, hjust = 0),
          plot.title = element_text(size = 16, hjust = -.05)) -> p
  
  # Device: save  image or display -------------------------------------------
  if(plot==T){
    ggsave(filename = paste0("meanVolPlot_",volume,".png"),
           path = here("results",paramT,"1_DescriptiveAnalysis",
                       "2_longitudinalMeasurements","meanVolumePlots_byAgeGroup"),
           width = 25,height = 12,units = "cm")
  }
  else{cat("Plotting mean plots for",volume,"...");p}
  
}

#----------------------------------------------------------------------
# ii) Plots sobre els 13 volums
#---------------------------------------------------------------------

meanPlot_long("CA1",data=dd,plot=F,window = 1, minN = 3)
lapply(volVar,FUN = meanPlot_long, data = dd, window = .8,
       minN = 3, plot=T)




# CHECKS =============================
## A cada timeGroup, hi ha més d'una mesura del mateix pacient?
plotData <- dd
volume <- "CA1"
window <- 0.25
# Prepare data-----------------------------------------------
thres <- seq(54,95,by=window) # thresholds
b <- cut(plotData$time_age,breaks=thres) #age groups
referenceAge <- data.frame(timeGroup =levels(b),
                           referenceAge = diff(thres)/2 + thres[1:length(thres)-1])
plotData$timeGroup <- b
plotData <- merge(plotData,referenceAge,by="timeGroup")
plotData <- plotData[,c("PTID","DX.bl","timeGroup","referenceAge",volume,
                        "VISCODE2","time_age")]
names(plotData)[5] <- "volume"

a <- summaryBy(volume~PTID+timeGroup,data=plotData,FUN=length)
subs <- a[a$volume.length>1,1:2]
nrow(subs)

b <- subset(plotData,PTID %in% subs$PTID & 
         timeGroup %in% subs$timeGroup ,c(PTID,VISCODE2))

merge(subs,plotData[,c("PTID","timeGroup","VISCODE2","time_age")],
      by=c("PTID","timeGroup"))
# 
















  
       


#-----------------------------------------------------------------------------
# ii) Gràfics de mitjanes i intervals de confiança en el temps de visita
#-----------------------------------------------------------------------------

### FUNCIO ------------------------------------------
meanPlot_timeVisit <- function(volumeVar,dd){
  
  # Summarise data
  sum <- summaryBy(list(c(volumeVar),c("time_visit","DX.bl")),
                   data=dd,
                   FUN=function(x){
                     n <- length(x)
                     if(n<3){return(c(n,mean(x),NA,NA))}
                     
                     m <- mean(x)
                     se <- sd(x)/sqrt(n)
                     t <- qt(0.975,df=n-1)
                     
                     return(c(n , m , m-t*se , m+t*se ))
                   })
  sum <- sum[order(sum$DX.bl,sum$time_visit),]
  colnames(sum) <- c("time","DX.bl","n","mean","lower","upper")
  
  #Plot: line plot. Mean volume value and 95% CI versus time after baseline
  ggplot(data=sum,aes(x=time,y=mean,color=DX.bl)) + 
    geom_line(aes(group=DX.bl),show.legend = F)+
    geom_point(show.legend = F) +
    geom_pointrange(aes(ymin=lower , ymax=upper),show.legend = F) +
    ggstatsplot::theme_ggstatsplot() +
    scale_color_manual(values = c("#1B9E77","#D95F02","#7570B3")) + 
    theme(strip.background = element_blank(),
          plot.margin = margin(1,1,.5,1,"cm"),
          axis.title = element_blank()) +
    ggtitle(paste0("Longitudinal evolution of volum in ",volumeVar,".\n  Mean value and 95% CI after basaeline visit")) +
    labs(subtitle = "<span style = 'color: #1B9E77;'>**Controls**</span>, <span style = 'color: #D95F02;'>**Mild cognitive impairment**</span> and <span style = 'color: #7570B3;'>**Alzheimer`s disease**</span>" ,
         caption="**Source**: ADNI project | **Diagnostic groups** CN: controls, AD: Alzheimer`s disease, MRI: mild cognitive impairment | **Units** Mean volume in mm3. Time after baseline visit in months") +
    theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
          plot.subtitle = element_markdown(size = 12, hjust = 0),
          plot.title = element_text(size = 16, hjust = -.05))
  
}

### Aplicar funció per totes les mesures volums

lapply(volVar,FUN=function(x){
  p <- meanPlot_timeVisit(volumeVar = x, dd = dd)
  
  filename <- paste0("meanPlotTimeVisit_",x,".png")
  ggsave(filename = filename,
         path = here("results",paramT,"1_DescriptiveAnalysis",
                     "2_longitudinalMeasurements","meanVolumePlots_byTimeVisit"),
         plot = p, device="png",width = 10,height = 5)
})


