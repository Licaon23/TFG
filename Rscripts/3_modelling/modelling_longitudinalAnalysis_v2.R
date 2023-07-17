###############################################################################
### MODELLING                                                                 #
###   Longitudinal Analysis analysis                                                      #
###   Version 1                                                               #
###   March   2023                                                            #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###                                                                           #
###############################################################################
source("Rscripts/libraries.R")
source(here(RscriptsDir,"3_modelling","modelling_auxFunctions.R"))
       
### SEL·LECCIONAR DATASET DE TREBALL: 1.5T O 3T ----------------------------

  # ## DATASET 3T ·········································
  # load(here(cleandataDir,"dataForAnalysis_3T.RData"))
  # paramT <- "Dataset_3T"
  #   #Treure única mesura individu AD al m60
  #   aux <- dd$DX.bl == "AD" & dd$VISCODE2=="m60"
  #   dd <- dd[!aux,]
  #   dd$PTID <- droplevels(dd$PTID)
  
  #·························································
  ### DATASET 1.5T
  load(here(cleandataDir,"dataForAnalysis_1T.RData"))
  paramT <- "Dataset_1T"
  #·························································
  
### VARIABLES DE TREBALL + TRANSFORMACIONS PREVIES --------------------------

  # variables: volums hipocamp, PRS continu, PRS categoric
  
  volVars <- sort(names(dd)[24:36])
  PRSVars <- names(dd)[10:16]
  riskVars <- names(dd)[37:43]

  # Data: seleccionem variables interès, i transformem temps de mesos a anys
  dd2 <- dd[c("PTID","AGE","DX.bl","PTGENDER","PTEDUCAT","ICV",
             "VISCODE2","timelapse","time_age","time_visit",
             riskVars, volVars)]
  dd2$time_visit <- dd2$time_visit/12
  dd2$timelapse <- dd2$timelapse/12
  head(dd2)


#########
# Seleccionem mesures longitudinals fins a 5 anys. Utilitzarem la variable
# time_visit com a mesura de temps longitudianal (com a temps transcorregut
# desde visita baseline, en anys). Seran 0, 0.5, 1, 2, 3, 4 i 5.
# És equivalent a prendre l'edat de cada individu a la mesura volum Yit
# i centrar-la a l'edat a l'inici de l'estudi (AGE, edat baseline)
# Mesures bl,m3,m6,m12,m24,m36,m48,m60
dd2 <- subset(dd2,VISCODE2 %in% c("bl","m06","m12","m24","m36","m48","m60"))
rownames(dd2) <- NULL



#===============================================================================
# 1) Nombre observacions per cada temps, segons grup
#===============================================================================

### i) CONSTRUIR DATAFRAME AMB NOMBRE OBSERVACIONS A CADA INSTANT (temps des de
#      mesura baseline en anys)
aux <- summaryBy(PTID ~ DX.bl+VISCODE2,
                 FUN=length,data=dd2)
names(aux)[3] <- "numMeasures"

numMeasures_wide <- reshape(aux, direction = "wide",
                            v.names = "numMeasures",
                            idvar = "VISCODE2",
                            timevar = "DX.bl")

names(numMeasures_wide) <- c("VISCODE2","CN","MCI","AD")
numMeasures_wide[is.na(numMeasures_wide)] <- 0

#numMeasures_wide

### ii) AFEGIR TEMPS DES DE BASELINE REAL ASSOCIAT A CADA INSTANT DE LA VARIABLE
#       TEMPORAL 'time_visit' : mitjana(sd) [min,max]
f <- function(x){
  stats <- c(m=mean(x),sd=sd(x),min(x),max(x))
  sprintf("%.2f (%.2f) [%.2f,%.2f]",stats[1],stats[2],stats[3],stats[4])
}

timelapseTab <- summaryBy(timelapse ~ VISCODE2,data=dd2,FUN=f)

taulaLongi <- merge(numMeasures_wide,timelapseTab,by="VISCODE2")
names(taulaLongi) <- c("Time Point","CN","MCI","AD","Time from baseline")

taulaLongi$`Time Point` <- c("Baseline", "Month 6",
                             paste("Year",1:5))


taulaLongi[1,"Time from baseline"] <- 0

#taulaLongi

### iii) CALCULAR NOMBRE DE SUBJECTES A INICI BASELINE
x <- basal$DX.bl
n <- sprintf("N = %.i",table(x))
n_header <- paste(paste(levels(dd$DX.bl),n,sep=", "),collapse="; ")


### IV) CONSTRUIR TAULA FORMATEJADA HTLM
tab <- kbl(taulaLongi,align = "lcccr",
           caption = paste0("Number and timing of scans per time point by diagnostic group (",n_header,")."))


tab <- kable_classic(tab,full_width=F,font="Cambria")
formattedTab <- footnote(tab,
                         general = c("Time from baseline (in years) is in mean (standard deviation); ranges are listed in square brackets.",
                                     "Diagnostic groups: Controls (CN), Mild Cognitive Impairment (MCI) and Alzheimer's Disease (AD) "),
                         title_format = c("bold"))

formattedTab_longiMeasures <- row_spec(formattedTab,0,bold=T)

# Guardar i visualitzar taula de resultats
file <- "longitudinalMeasures_descriptiveTable.png"
path <- here("results",paramT,"2_Modelling",
             "2_longitudinalTrajectories","taules")

if(!(any(grepl(file,list.files(path))))){
  save_kable(formattedTab_longiMeasures,file=here(path,file),zoom=2)
}

#formattedTab_longiMeasures

### Nombre de mesures per participant
aux <- sapply(split(dd2,dd2$PTID),nrow)
measures <- data.frame(PTID=names(aux),numObs = aux)
measures <- merge(measures,basal[c("PTID","DX.bl")],by="PTID")
obsTab <- summaryBy(numObs ~ DX.bl, data=measures, FUN=summary)
names(obsTab) <- c("Diagnostic group","Min","Q1","Q2","Mean","Q3","Max")
obsTab


#===============================================================================
# 2) GRÀFICS DE TRAJECTÒRIES INDIVIDUALS I MITJANES 
#===============================================================================


# Mean trend
pos <- position_dodge(width=0.2)

plotData <- dd2[,c("PTID","DX.bl","time_visit",volVars)]
plotData <- reshape(plotData,
                    varying = names(plotData[-c(1:3)]),
                    timevar = "subfield",
                    times = names(plotData[-c(1:3)]),
                    idvar = c("PTID","time_visit"),
                    v.names = "volume",
                    direction = "long")
rownames(plotData) <- NULL

ggplot(plotData,aes(x=time_visit,y=volume,group=DX.bl)) +
  geom_smooth(aes(linetype=DX.bl),method = "loess",se=F,color="black") +
  stat_summary(fun.data= "mean_cl_boot",size=.4)+
  scale_y_continuous(expand=c(0,0)) +
  scale_linetype_discrete("Diagostic group")+
  theme_classic()+
  facet_wrap(~as.factor(subfield),scales = "free_y") +
  ggtitle("Smoothed Mean measurement trajectories for different hippocampal subfield volumes")+
  labs(y = "Hippocampal volume (mm3)",
       caption="Hippocampal subfield mean volume at each measurement time after baseline (years), and 95% confidence interval.<br/>
                Smoothed trajectory was plotted using loess curve for each diagnositc group.<br/>
                It should be taken into account that subjects at baseline have different age.<br/>
       **Source**: ADNI project | **Units** Average hippocampal volume (mm3) of subjects at each measurement time after baseline (in years)") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = 0),
        axis.title.x = element_blank())






#===============================================================================
# 3) Anàlisi del canvi en el temps dels volums hipocampals
# --------------------------------------------------------------------------
#    OBJECTIUS:
#
#     0) Comprovar diagnòstic de residus i detectar outliers.
#
#     1) Comprovar que existeix una reducció significativa al llarg del temps 
#        dels volums de l'hipocamp (Model 1)
#
#     2) Veure si aquesta trajectòria és o no diferent segons el grup 
#        diagnòstic (Model 2)
#
#     3) Amb el model 2 com a efectes fixos i random slope i intercept com
#        a aleatoris, contrastar si el random slope millora el model.
#
#     4) Verificar si en algun dels volums el terme quadràtic pel temps és
#         significatiu, per considerar no-linealitat en la trajectòria.
#         (Model 3)
#
#     5) Comprovar si existeix un efecte específic del sexe sobre el ritme de
#        reducció dels volums hipocampals diferents per grup diagnòstic.
#
#
#------------------------------------------------------------------------------
#   PROCEDIMENT:
#
#     Obj.0: revisar els gràfics diagnòsitc i analitzar subjectes que tinguin
#            observacions amb residus elevats.
#
#     Obj.1: ajustar un model lineal mixte per cada volum amb time_visit
#            com a variable temporal i corregint per grup diagnòstic (DX.bl),
#            edat basal (AGE), sexe (PTGENDER), anys educació (PTEDUCAT) i volum
#            intracranial total en cm3 (ICV/1000).
#            Per considerar la correlació entre les mesures d'un mateix subjecte
#            es permet que cada individu tingui un punt de partida i pendent
#            temporal propi (random intercept i slope segons PTID).
#            (Model additiu, 1):
#
#         Vol_i ~ time_visit + DX.bl + Covariables + (1+time_visit|PTID)
#
#     Obj.2: considerar al model anterior que l'efecte de time_visit sobre
#            el volum resposta depèn del grup diagòsitc, afegint una interacció
#            time_visit * DX.bl. (Model diagnositc, 2)
#
#         Vol_i ~ time_visit * DX.bl + Covariables + (1+time_visit|PTID)
#
#     Obj.3: aplicar un test de raó de versemblança per testar 
#            els models aniuats amb i sense random slope, basat en un
#            estadístic mixtura de ji-quadrats (funció ranova() de lmerTest)
#            Corregir per comparacions múltiples segons FDR al 5%
#
#     Obj.4: afegir el terme quadràtic de time_visit tant com a efecte principal
#            com amb interacció amb DX.bl. Realitzar el contrast simultani
#            de que tots els coeficients que involucrin el terme quadràtic
#            del temps són 0. Se'n recupera el p-valor del contrast. Si és
#            inferior a 0.05 es considera que algun dels coeficients del terme
#            quadràtic té un efecte. (Model 3)
#            Corregir per comparacions múltiples amd FDR 5% pels 13 volums.
#
#         Vol_i ~ (time_visit + time_visit^2) * DX.bl + Covariables + 
#                  (1+time_visit|PTID)
#
#     Obj.5: afegir un terme interacció time_visit*PTGENDER al model 2 i 
#            contrastar la significació del coeficient.
#
#         Vol_i ~ time_visit * DX.bl + time_visit * PTGENDER + 
#           Covariables + (1+time_visit|PTID)
#
#
#   NOTA: covariables són PTGENDER, PTEDUCAT, AGE, ICV/1000
#         Efectes aleatoris: random intercept i pendent time_visit segons PTID
#         No es considera efecte quadràtic del time_visit, ja que no s'ha
#         trobat significatiu a l'apartat 1.
#===============================================================================

#--------------------------------------------------------------------------------
# FUNCIÓ timeEffect(): per cada volum, realitza els procediments pels quatre
# objectius anteriors, i retorna com a resultat una llista d'objectes per
# respondre'ls.
#--------------------------------------------------------------------------------

# Valor per debuggin funció
# volume <- "CA1"; data<-dd2

timeEffect <- function(volume,data){
  # OBJECTIU 1: Ajustar model cru································
  predictors <- "time_visit + DX.bl + AGE + I(ICV/1000) + PTGENDER + PTEDUCAT + (1+time_visit|PTID)"
  form <- paste(volume,predictors, sep = " ~ ")
  mod1 <- lme4::lmer(formula = form, data=data,
                     control = lmerControl(optimizer ="Nelder_Mead"))
  
  # OBJECTIU 2: Model amb interacció temps*DX.bl························
  predictors <- "time_visit * DX.bl + AGE + I(ICV/1000) + PTGENDER + PTEDUCAT + (1+time_visit|PTID)"
  form <- paste(volume,predictors, sep = " ~ ")
  mod2 <- lme4::lmer(formula = form, data=data,
                     control = lmerControl(optimizer ="Nelder_Mead"))
  
  # OBJECTIU 3: LRT test pels efectes aleatoris ··························
  ranovaTest <- ranova(mod2)$`Pr(>Chisq)`[2]
  
  # OBJECTIU 4: Terme quadràtic pel temps?······························
  mod3 <- update(mod2,~.+I(time_visit^2)+I(time_visit^2):DX.bl)
  posTerms <- grep("time_visit\\^2",names(fixef(mod3)))
  terms <- names(fixef(mod3))[posTerms] # coeficients time_visit^2 
  
  # Contrast H0: termes time_visit^2 són tots 0, aka, no hi ha efecte
  # quadràtic significatiu del terme temporal time_visit
  timeQuad <- car::linearHypothesis(mod3,terms,c(0,0,0))$`Pr(>Chisq)`[2]
  
  # OBJECTIU 5: efecte diferencial del sexe sobre el time rate of change? ····
  predictors <- "time_visit * DX.bl + time_visit * PTGENDER + AGE + I(ICV/1000) + PTEDUCAT + (1+time_visit|PTID)"
  form <- paste(volume,predictors, sep = " ~ ")
  mod4 <- lme4::lmer(formula = form, data=data,
                     control = lmerControl(optimizer ="Nelder_Mead"))
  
  # OBJECTIU 0: Residus··············································
  #resid <- plot(mod2,main = paste("Diagnositc",volume))
  resid <- resid_panel(mod2)
  
  
  # Resultats························································
  return(list(modelAdditiu = mod1,
              modelDiagEffect = mod2,
              modelSexEffect = mod4,
              ranovaTest = ranovaTest,
              efecteQuadratic = timeQuad,
              diagnostic = resid))
}
#--------------------------------------------------------------------------------

### APLICAR FUNCIO MODELS PER TOTS ELS VOLUMS HIPPOCAMPALS
timeEffectModels <- lapply(volVars,FUN=timeEffect, data=dd2)
names(timeEffectModels) <- volVars


#-------------------------------------------------------------------------------
# A) OBJECTIU 0
#    Revisar els gràfics diagnòstic dels models
#-------------------------------------------------------------------------------

lapply(timeEffectModels,"[[","diagnostic")

# RESPOSTA:
# En general, els residus 'within-group' semblen normals, centrats al
# 0 i amb variància constant, i distribuits aleatòriament.
# Tanmateix, sistemàticament sembla que tres observacions presenten
# un error molt elevat. També en la distribució dels efectes aleatoris.

# De quins individu son les mesures outliers?
outliers <- lapply(timeEffectModels,function(x){
  car::outlierTest(x$modelDiagEffect)
})

aux <- lapply(outliers,function(x)as.integer(names(x$rstudent))) #poisicio mesura
aux <- lapply(aux,function(x)as.character(unique(dd2[x,"PTID"]))) #ptid corresponent
(tab <- table(do.call(c,aux))) # nombre de vegades cada ptid presenta mesura outlier
(ptid_outlier <- names(which.max(tab))) # ptid amb mes observacions outliers

# Incoherències en les mesures d'aquest individu... 
dd2[dd2$PTID == ptid_outlier,c("PTID","DX.bl","time_visit",volVars)]

# DECISIÓ: excloure les mesures d'aquest participant
# (Caldrà tornar a correr l'anàlisi anterior amb aquest pacient exclòs)
dd2 <- subset(dd2, PTID != ptid_outlier)
dd2$PTID <- droplevels(dd2$PTID)
rownames(dd2) <- NULL
#length(unique(dd2$PTID))


### APLICAR FUNCIO MODELS PER TOTS ELS VOLUMS HIPPOCAMPALS
timeEffectModels <- lapply(volVars,FUN=timeEffect, data=dd2)
names(timeEffectModels) <- volVars
lapply(timeEffectModels,"[[","diagnostic")

## Save diagnostics
lapply(names(timeEffectModels),FUN=function(x){
  path <- here("results",paramT,"2_Modelling","2_longitudinalTrajectories",
               "modelDiagnostics")
  file <- paste0(path,"/modelResiduals_",x,".png")  
  ggsave(filename = file, plot=timeEffectModels[[x]]$diagnostic,device = "png")
})

#-------------------------------------------------------------------------------
# A) OBJECTIU 1
#    Hi ha un canvi en el temps del volums hippocampals?
#--------------------------------------------------------------------------------
rateChangeVolumes <- 
  lapply(timeEffectModels,function(x){
    fixef <- extract_fixed_effects(x$modelAdditiu)[2,]
    effect <- sprintf("%.2f (%.2f,%.2f); %.3f", 
                      fixef$value,fixef$lower_2.5,fixef$upper_97.5,fixef$p_value)
    
    return(effect)
  })

rateChangeVolumes <- do.call(c,rateChangeVolumes)
rateChangeVolumes <- data.frame(timeEffect = rateChangeVolumes)

rateChangeVolumes

# RESPOSTA:
#  Per a tots els subcamps de l'hipocamp s'observa una reducció del volum
#  significativa.


#-------------------------------------------------------------------------------
# B) OBJECTIU 2
#    Hi ha alguna diferècia en el canvi de volum segons el grup diagnòstic?
#-------------------------------------------------------------------------------
effectDXbl <- 
  sapply(timeEffectModels, function(x){
    pval_contrast <- car::linearHypothesis(x$modelDiagEffect,
                                           c("time_visit:DX.blMCI","time_visit:DX.blAD"),
                                           c(0,0))[2,"Pr(>Chisq)"]
    return(round(pval_contrast,3))
  })

effectDXbl 
round(p.adjust(effectDXbl,method = "fdr"),3)

# RESPOSTA:
# Sí, la trajectòria del volum en el temps és diferent segons grup
# diagnòsitc, sent el declivi més accentuat en els individus MCI i AD respecte
# els controls sans.

timeTrendDx <- lapply(timeEffectModels,function(volume){
  mod <- volume$modelDiagEffect
  
  timeTrend <- emtrends(mod,specs = "DX.bl",var="time_visit",
                        at=list(PTGENDER="Female"))
  contrast <- pairs(timeTrend)
  
  # Time trend
  aux <- t(summary(timeTrend)[,-1])
  slopes <- apply(aux,2,function(x){sprintf("%.2f (%.2f,%.2f)",x[1],x[4],x[5])})
  names(slopes) <- summary(timeTrend)[,1]
  
  #Contrast
  pvals <- summary(contrast)$p.value
  pvals <- ifelse(pvals<0.001,"<.001",sprintf("%.3f",pvals))
  names(pvals) <- summary(contrast)$contrast

  # Result
  res <-as.data.frame(t(matrix(c(slopes,pvals))))
  names(res) <- Reduce("c",names(c(slopes,pvals)))
  return(res)
})

timeTrendDx <- do.call(rbind,timeTrendDx)
timeTrendDx

#Taula Formatejada dels trends
tab <- kbl(timeTrendDx,align = "cccccc",
           caption = "Estimated mean annual volume changes by group diagnostic")


tab <- kable_classic(tab,full_width=F,font="Cambria")
tab <- add_header_above(tab,c(" " = 1,
                              "Estimated linear time trend"=3,
                              "Pairwise comparisons"=3), 
                        bold = T,line = T)
formattedTab <- footnote(tab,
                         general = c("Linear time trend (95% CI) as linear combination of time-related coefficients. It indicates mean annual volume change in mm3 estimated for a given diagnostic
                                     group and hippocampal region. Differences in time trend between diagnostic groups were assessed by Wald tests for linear contrasts. Pvalues are reported",
                                     "Diagnostic groups: Controls (CN), Mild Cognitive Impairment (MCI) and Alzheimer's Disease (AD) "),
                         title_format = c("bold"))

formattedTab_timeTrendDX <- row_spec(formattedTab,0,bold=T)
formattedTab_timeTrendDX

# # Save paràmetres dels models
#   params <- lapply(timeEffectModels,function(x){
#     extract_fixed_effects(x$modelDiagEffect)
#   })
#   path <- here("results",paramT,"2_Modelling","2_longitudinalTrajectories",
#                "taules","modelL0_Parameters.xlsx")
#   write.xlsx(x=params,file=path)
  
#-------------------------------------------------------------------------------
# C) OBJECTIU 3
#    El pendent aleatori en el temps millora els models?
#-------------------------------------------------------------------------------

ranovaTest <- do.call(c,lapply(timeEffectModels,"[[","ranovaTest"))
p.adjust(ranovaTest,method = "fdr")

#-------------------------------------------------------------------------------
# D) OBJECTIU 4
#    Hi ha un efecte quadràtic del temps?
#-------------------------------------------------------------------------------

quadraticTimeEffect <- do.call(c,lapply(timeEffectModels,"[[","efecteQuadratic"))
round(p.adjust(quadraticTimeEffect,method = "fdr"),3)

# RESPOSTA:
# Després de corregir per 13 comparacions, cap dels volums mostra un efecte
# no lineal quadràtic significatiu en el canvi temporal.

# PER TANT: els models representen evolució lineal en el temps dels volums
#   hipocampals, amb trajectòries diferents segons grup diagnòsitc i
#   intercept i pendent temporal aleatories per cada individu.

# Com que interessa veure diferències en les trajectòries segons grup diagnòstic
# i altres factors, els propera anàlisis centren la seva atenció als termes
# de interacció amb el temps: time_visit*DX.bl  i  time_visit * PRS

#-------------------------------------------------------------------------------
# E) OBJECTIU 5
#    Sex effect sobre el rate of change: beta_timeVisit:PTGENDER
#-------------------------------------------------------------------------------
# Per debuggin: x <- timeEffectModels[[1]]



timeTrendSx <- lapply(timeEffectModels,function(volum){
  mod <- volum$modelSexEffect
  
  timeTrend <- emtrends(mod,specs = "PTGENDER",var="time_visit")
  contrast <- pairs(timeTrend)
  
  # Time trend
  aux <- t(summary(timeTrend)[,-1])
  slopes <- apply(aux,2,function(x){sprintf("%.2f (%.2f,%.2f)",x[1],x[4],x[5])})
  names(slopes) <- summary(timeTrend)[,1]
  
  #Contrast
  pvals <- summary(contrast)$p.value
  pvals <- ifelse(pvals<0.001,"<.001",sprintf("%.3f",pvals))
  names(pvals) <- summary(contrast)$contrast
  
  # Result
  res <-as.data.frame(t(matrix(c(slopes,pvals))))
  names(res) <- Reduce("c",names(c(slopes,pvals)))
  return(res)
})

timeTrendSx <- do.call(rbind,timeTrendSx)
timeTrendSx

# RESPOSTA:
# Retorna coeficient time_visit segons sexe, amb interval de confiança pel
# mètode de Wald i pvalor de significació del contrast. 
# Com a mètode per calcular els graus de llibertat aproximats, 
# s'utilitza el mètode de Satterthwaite. Correció de tukey

# No sembla que hi hagi un efecte específic del sexe en el rate de pèrduda 
# de volum dels diferents subcamps de l'hipocamp.





#===============================================================================
# 4) Anàlisi de l'efecte de la predisposició genètica sobre la trajectòria
#    de reducció dels volums hipocampals.
#    --------------------------------------------------------------------------
#    OBJECTIUS:
#     1) Comprovar si la predisposició genètica a trets neurodegeneratius i
#        al envelliment influeixen en el ritme de reducció dels volums hipocampals.
#
#     2) Si existeix aquest efecte genètic time_visit*PRS, és diferent segons 
#        l'estat de l'individu en el continu de la malaltia, DX.bl?
#
#     3) Aquest efecte genètic és diferent segon el sexe?
#
#------------------------------------------------------------------------------
#   PROCEDIMENT:
#     Assumpció: es considera que la reducció de volum al llarg de les cinc
#               mesures longitudinals segueix una tendència lineal.
#
#     Obj.1: ajustar un model lineal mixte per cada volum amb:
#              - time_visit com a variable temporal, 
#              - un determinat PRS com a variable explicativa principal, amb
#                interacció amb el temps, (time_visit * PRS),
#              - Ajustant per covariables per grup diagnòstic (DX.bl),
#                edat basal (AGE), sexe (PTGENDER), anys educació (PTEDUCAT) i 
#                volum intracranial total en cm3 (ICV/1000).
#            La trajectòria temporal depèn de l'estat diagnòstic, com s'ha vist,
#            time_visit * DX.bl.
#            Per considerar la correlació entre les mesures d'un mateix subjecte
#            es permet que cada individu tingui un punt de partida i pendent
#            temporal propi (random intercept i slope segons PTID).
#            El terme d'interès és el coeficient interacció time_visit*PRS.
#
#         Vol_i ~ time_visit * DX.bl + time_visit * PRS + 
#                   Covariables + (1+time_visit|PTID)
#
#     Obj.2: afegir interacció triple time_visit*PRS*DX.bl? 
#            O estratificar per DX.bl i estudiar efecte time_visit*PRS?
#
#         Vol_i ~ time_visit * DX.bl * PRS + Covariables + (1+time_visit|PTID)
#
#         * Interessa el coeficient associat al temps (time trend), segons
#           quin sigui el diagnòstic i per a un valor donat del PRS (0 o 1)
#           PER EXEMPLE:
#           i) TimeTrend|DX.bl=CN & PRS=Low-Intermediate
#               TimeTrend = beta_timeVisit
#
#           ii) TimeTrend|DX.bl=CN & PRS=High
#                 TimeTrend = beta_timeVisit + beta_timeVisit:PRS_higj
#
#           iii) TimeTrend|DX.bl=MCI & PRS=Low-Intermediate
#                 TimeTrend = beta_timeVisit + beta_timeVisit:DX.blMCI
#
#           iv)  TimeTrend|DX.bl=MCI & PRS=High
#                 TimeTrend = beta_timeVisit + beta_timeVisit:DX.blMCI + 
#                               beta_timeVisit:Dx.blMCI:PRS_High
#
#         * Aquestes combinacions lineals de coeficients es calculen amb la
#           funció emtrends del paquet emmeans (que es basa amb la funció
#           glht del paquet multcomp per). Es retora el valor d'aquestes 
#           combinacions betaTrend + errors estàndard + significació per test
#           de Wald.
#
#     Obj.2a: anàlisi estratificat. time_visit*PRS estratificant per DX.bl
#
#     Obj.3: ??
#===============================================================================

#--------------------------------------------------------------------------------
# FUNCIÓ prsTimeEffect(): per cada volum i PRS continu es realitza els 
# el procediment per objectiu 1.
# A més, també s'ajusta el model amb la interacció triple i es se'n retorna
# la significació del test F: model amb interacció triple (complet) vs model
# restringit sense interacció triple time_visit*PRS*DX.bl
#--------------------------------------------------------------------------------

# Si es vol provar linea per linea la funció...
# volume <- "Hippocampal_tail";PRS <- "AD_predisposition";data <- dd2
#-------------------------------------------------------------------------------
prsTimeEffect <- function(volume,PRS,data){
  
  # Formula ··································
  predictors <- paste("time_visit",c("DX.bl",PRS),sep=" * ")
  predictors_3intDX <- paste("time_visit","DX.bl",PRS,sep=" * ") #triple interact DX
  predictors_3intSx <- paste("time_visit","PTGENDER",PRS,sep=" * ") #triple interact Sex
  #covars <- c("AGE","ICV","DX.bl","PTGENDER","PTEDUCAT")
  covars <- c("AGE","I(ICV/1000)","DX.bl","PTGENDER","PTEDUCAT")
  random <- "(1 + time_visit|PTID)"
  
  form <- 
    sapply(list(predictors,predictors_3intDX,predictors_3intSx),function(x){
      res <- paste(c(x,covars,random),collapse= " + ")
      res <- paste(volume,res, sep = " ~ ")
      return(res)
    })
  
  ## OBJECTIU 1: Model amb interacció time_visit*PRS ····················
  mod <- lme4::lmer(formula = form[1], data=data,
                    control = lmerControl(optimizer ="Nelder_Mead"))
  
  attributes(mod)$PRS <- PRS
  attributes(mod)$volume <- volume
  
  
  
  # Extraure efecte PRS sobre variacó volum en el temps 
  fixef <- extract_fixed_effects(mod)
  effect_pos <- which(fixef$term == paste0("time_visit:",PRS,"High"))
  
  # Resultat en un dataframe: valor, 95%CI i pvalor test wald
  # amb graus de llibertat per Satterthwaite
  res <- as.data.frame(fixef[effect_pos,c("value","se","lower_2.5",
                                          "upper_97.5","p_value")])
  
  rownames(res) <- paste(volume,PRS,sep=":")
  #······································································
  
  ## OBJECTIU 2: Model amb interacció triple time_visit*PRS*DX.bl·················
  mod_3IntDX <- lme4::lmer(formula = form[2], data=data,
                     control = lmerControl(optimizer ="Nelder_Mead"))
  
  testDX <- Anova(mod_3IntDX,type="III")[12,]
  attributes(mod_3IntDX)$PRS <- PRS
  attributes(mod_3IntDX)$volume <- volume
  
  #······································································
  
  ## OBJECTIU 3: Model amb interacció triple time_visit*PRS*PTGENDER·················
  form[3] <- sub("DX.bl","time_visit * DX.bl",form[3])
  mod_3IntSx <- lme4::lmer(formula = form[3], data=data,
                     control = lmerControl(optimizer ="Nelder_Mead"))
  
  testSx <- Anova(mod_3IntSx,type="III")[13,]
  attributes(mod_3IntSx)$PRS <- PRS
  attributes(mod_3IntSx)$volume <- volume
  
  ## Retornar resultats: llista amb efecte PRS*time_visit + test anova + 
  ## model amb interacció time_visit * PRS * DX.bl
  return(list("modGlobal"= mod,
              "PRS_timeEffect" = res,
              "AnovaTest_3wayIntDX" = testDX,
              "mod_3wayIntDX" = mod_3IntDX,
              "AnovaTest_3wayIntSx" = testSx,
              "mod_3wayIntSx" = mod_3IntSx))
}
#--------------------------------------------------------------------------
prsTimeEffect("subiculum","AD_predisposition",data=dd2)



#-------------------------------------------------------------------------------
# A) OBJECTIU 1
#    Efecte del PRS sobre la reducció dels volums?
#--------------------------------------------------------------------------------

### i) Matriu de combinacions PRS i volums hipocamp
modGrid <- expand.grid(PRS = riskVars,
                       Volume = sort(volVars),
                       stringsAsFactors = F)
# head(modGrid)

### ii) Apliquem la funció effectPRS_strat() sobre la GRID 
result <- mapply(FUN = prsTimeEffect,
                 volume = modGrid$Volume,
                 PRS = modGrid$PRS,
                 MoreArgs = list(data=dd2),
                 SIMPLIFY = F,USE.NAMES = F)

names(result) <- paste(modGrid$Volume,modGrid$PRS,sep=":")


# i arreglem el resultat...
res_prsTimeEffect <- lapply(result,"[[","PRS_timeEffect") # seleccionar efecte PRS de la list
res_prsTimeEffect <- do.call(rbind,res_prsTimeEffect) # consolidar en dataFrame
res_prsTimeEffect <- cbind(modGrid,res_prsTimeEffect) # Afegir PRS i volume 
rownames(res_prsTimeEffect) <- NULL # no rownames
res_prsTimeEffect$PRS <- factor(res_prsTimeEffect$PRS,levels=riskVars) 

head(res_prsTimeEffect,12)

### iii) Corregim pvalors segons FDR, considerant cada volum independent,
#       i corregint per 6 comparacions (6 PRS) o models per cada un dels volums

aux <- lapply(split(res_prsTimeEffect,res_prsTimeEffect$Volume),
              FUN=function(x)p.adjust(x$p_value,method="fdr"))
res_prsTimeEffect$pvalue_corrected <- do.call(c,aux)

### iv) Afegim simbols de significació
res_prsTimeEffect[paste(names(res_prsTimeEffect)[7:8],"symbol",sep="_")] <-
  lapply(res_prsTimeEffect[names(res_prsTimeEffect[7:8])], FUN=cut,
         breaks = c(-Inf,0.001,0.01,0.05,0.1,1),
         labels=c("***","**","*","·",""))

# RESULTAT
head(res_prsTimeEffect)
(significant_effects <- subset(res_prsTimeEffect, pvalue_corrected_symbol!=""))

# Taula
tab <- significant_effects[c(2,1,3:6,8)]
names(tab) <- c("Volume","PRS","Beta coef","SE","lower_CI","upper_CI","pvalue")
tab$pvalue <- ifelse(tab$pvalue<0.001,"<0.001",round(tab$pvalue,3))

tab <- kbl(tab,align = "llcccccc",row.names = F,
           caption = "Significant effects of PRS on time trend. Models L.1")


tab <- kable_classic(tab,full_width=F,font="Cambria")
formattedTab <- footnote(tab,
                         general = c("Significant linear coefficients for interaction between PRS and time. Pvalues were corrected for multiple comparisons using FDR method."),
                         title_format = c("bold"))

formattedTab_prsTimeEffect <- row_spec(formattedTab,0,bold=T)
formattedTab_prsTimeEffect
save_kable(x = formattedTab_prsTimeEffect,
           file = here("results",paramT,"2_Modelling",
                       "2_longitudinalTrajectories","taules",
                       "modelL1_betasPRSTime.png"),zoom=1.5)



### v) HEAT MAP amb els resultats
plotData <- subset(res_prsTimeEffect,Volume!="Whole_hippocampus")

heatMap <- 
  ggplot(data=plotData,
         aes(x=Volume, y=PRS, fill=value)) + 
  geom_tile() +
  geom_text(aes(label=pvalue_corrected_symbol),size=6)+
  scale_fill_gradient2(low="red4",mid = "white", high="#0066FF", midpoint = 0)+
  # facet_wrap(~DX.bl)+
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank(),
        axis.text = element_text(face="bold",size = 8),
        axis.text.x = element_text(angle=90,vjust = 0.6,hjust = 1)) +
  
  ggtitle("Corrected significant effect of investigated PRS on hippocampal\nsubfield volume's rate of change over time") +
  labs(caption="Linear coefficient for interaction between time and PRS. Effects adjusted for Age at baseline (years), intracranial volume(cm3), sex, diagnostic<br/>
  status and years of education. Model allows a diferent volume rate of change depending on diagnositc status. Also, intercept and rate of<br/>
  change in time are included as random effects. Significant interactions are reported as Wald-test pvalues at levels 0.1(·), 0.05(\\*), 0.01(\\**)<br/>
  and 0.001(\\***) and corrected for multiple comparisons by FDR method. Non-significant differences are not shown.<br/>
  **Source**: ADNI project | **Units** Average increment in volume rate of change (mm3) associated with having high predisposition for the PRS trait.") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = 0))

heatMap
# plotname <- here("results",paramT,"2_Modelling","2_longitudinalTrajectories",
#              "modelsL1_PRS_time_interaction.png")
# ggsave(filename = plotname, width=894,height=469,units="px",zoom=0.5)




#-------------------------------------------------------------------------------
# B) OBJECTIU 2
#    Efecte del PRS sobre la reducció dels volums, depèn de l'estat diagnòsitc?
#-------------------------------------------------------------------------------

# Interacció time_visit*PRS*DX.bl significativa?
testDX <- lapply(result,"[[","AnovaTest_3wayIntDX")
testDX <- lapply(testDX,"[[",'Pr(>Chisq)')
testDX <- do.call(c,testDX)
res_prsTimeEffect$pval_3wayIntDX <- round(testDX,3)

testDX[testDX<0.05]
subset(res_prsTimeEffect,pval_3wayIntDX<0.05,c(1:3,11,8:10))



### CALCULAR ELS TRENDS LINEALS SEGONS predisposition low or high ···················
# Debbuging:x<-result$`Hippocampal_tail:AD_predisposition` # per provar


trends <- lapply(result,function(x){
  mod <- x$mod_3wayIntDX
  PRS <- attributes(mod)$PRS
  specs <- as.formula(paste0("~",PRS,"*DX.bl"))
  at <- "list(PTGENDER='Female')"
  
  trend <- emtrends(mod, specs = specs, "time_visit",
                    at = eval(parse(text=at)))
  
  return(trend)
})

# A dataframe i canviar nom columna PRS a genèric per despres fer el merge
trends.df <- lapply(trends,summary)

trends.df <- lapply(trends.df,function(x){
  names(x)[1] <- "PRS_value"
  return(x)
})

trends.df <- do.call(rbind,trends.df)


### CALCULAR CONTRASTOS High-Low ···················
contrastPRS <- lapply(trends,function(x){
  n <- names(x@levels)[1]
  summary(contrast(x,"consec",simple=n))
})

contrastPRS <- as.data.frame(do.call(rbind,contrastPRS))

res_contrastPRS <- data.frame(Volume = rep(volVars,each=3*7),
                              PRS = rep(rep(riskVars,each=3),times=13))

res_contrastPRS <- cbind(res_contrastPRS,contrastPRS)

# Corregir per comparacions múltiples FDR
aux <- lapply(split(res_contrastPRS,res_contrastPRS$Volume),
              FUN=function(x)p.adjust(x$p.value,method="fdr"))

res_contrastPRS$pvalue_corrected <- do.call(c,aux)

# Afegir simbols de significació
res_contrastPRS$pvalue_corrected_symbol <-
  cut(res_contrastPRS$pvalue_corrected,
      breaks = c(-Inf,0.001,0.01,0.05,0.1,1),
      labels=c("***","**","*","·",""))

head(res_contrastPRS)

contrast1_signif <- subset(res_contrastPRS,pvalue_corrected<0.05)
prsDiff <- with(contrast1_signif,paste(Volume,PRS,sep=":"))

# Exportar taula significatius formatejada
contrast1_signif[5:8] <- lapply(contrast1_signif[5:8],round,2)
contrast1_signif[9:10] <- lapply(contrast1_signif[9:10],function(x){
  ifelse(x<0.001,"<0.001",sprintf("%.3f",x))
})

contrast1_signif$pvalue_corrected <- 
  with(contrast1_signif,paste0(pvalue_corrected,"(",pvalue_corrected_symbol,")"))

rownames(contrast1_signif) <- NULL
contrast1_signif$pvalue_corrected_symbol <- NULL

tab <- kbl(contrast1_signif,align = "lllccccccc",row.names = F,
           caption = "Significant effects of PRS on time trend. Models L.2, linear contrasts 1")


tab <- kable_classic(tab,full_width=F,font="Cambria")
formattedTab <- footnote(tab,
                         general = c("Wald test for $H_0: \\beta5_{time:High} + \\beta7_{time:High:MCI}=0$.", 
                                     "Pvalues were corrected for multiple comparisons using FDR method. Significance level: (*)<.05, (**)<.01, (***)<.001"),
                         title_format = c("bold"))

formattedTab_contrast1 <- row_spec(formattedTab,0,bold=T)
formattedTab_contrast1

# save_kable(x = formattedTab_contrast1,
#            file = here("results",paramT,"2_Modelling",
#                        "2_longitudinalTrajectories","taules",
#                        "modelL2_contrast1Signif.png"),zoom=1.5,)

# Exportar excel contrastos 1
  tab <- res_contrastPRS
  rownames(tab) <- NULL
  tab <- split(tab,tab$Volume)

  # Build Excel worksheet
  wb <- createWorkbook()
  invisible(lapply(names(tab),function(x)addWorksheet(wb,x)))
  invisible(
    lapply(1:length(tab),
           function(i) writeDataTable(wb = wb,
                                      x = tab[[i]],
                                      sheet = names(tab)[i],
                                      withFilter = T,
                                      tableStyle = "TableStyleLight1")
    )
  )

path <- here("results",paramT,"2_Modelling","2_longitudinalTrajectories",
             "taules","modelL2_contrast1.xlsx")
saveWorkbook(wb,file=path,overwrite = T)

# Diferències entre grups quant a efecte genètic?

aux <- lapply(prsDiff,function(x){
  
  linfct <- contrast(trends[[x]],interaction=c("consec","pairwise"))@linfct
  mod <- result[[x]]$mod_3wayIntDX
  joinTest <- summary(glht(mod,linfct),test=Chisqtest())
  
  res <- data.frame(Volumes = attributes(mod)$volume,
                    PRS = attributes(mod)$PRS,
                    W2_chisq = as.numeric(joinTest$test$SSH),
                    df = as.numeric(joinTest$test$df[[1]]),
                    Pval = joinTest$test$pvalue)
  res$symbolJoinTest <- ifelse(res$Pval<0.1,"*","")
  return(res)
})
contrast2 <- do.call(rbind,aux)

# Exportar taula contrast 2 formatejada

tab <- contrast2
tab$W2_chisq <- round(tab$W2_chisq,2)
tab$Pval <- ifelse(tab$Pval<0.001,"<0.001",sprintf("%.3f",tab$Pval))
tab$Pval <- with(tab,ifelse(symbolJoinTest=="*",paste0(Pval,"(*)"),Pval))
tab$symbolJoinTest <- NULL

tab <- kbl(tab,align = "llccc",row.names = F,
           caption = "Differences between diagnostic groups regarding PRS effects on volume rate of change. Models L2, linear contrasts 2")


tab <- kable_classic(tab,full_width=F,font="Cambria")
formattedTab <- footnote(tab,
                         general = c("Simultaneous Wald Test for 3 linear contrast. $H_0:\\beta7_{time:MCI:High}=0 ;\\beta7_{time:AD:High}=0; \\beta7_{time:AD:High}-\\beta7_{time:MCI:High}=0$",
                                     "Significance level: (*)<.05"),
                         title_format = c("bold"))

formattedTab_contrast2 <- row_spec(formattedTab,0,bold=T)
formattedTab_contrast2




### DATAFRAME DE RESULTATS PER PLOTEJAR ···································

finalResult <- data.frame(Volumes = rep(volVars,each = 6*length(riskVars)),
                          PRS = rep(riskVars,each=6))
finalResult <- cbind(finalResult,trends.df)

# Afegim significació corregida (el mateix pvalor pel low i high en cada
# volum-PRS-DX.bl)
finalResult$pvalue <- rep(res_contrastPRS$p.value,each=2)
finalResult$pvalue_corrected <- rep(res_contrastPRS$pvalue_corrected,each=2)


# Transparency per significació del contrast PRS
finalResult$alpha <- with(finalResult,ifelse(pvalue_corrected<0.05,1,.9))


# JoinTest
finalResult <- merge(finalResult,
                     contrast2[c("Volumes","PRS","symbolJoinTest")],
                     by=c("Volumes","PRS"),
                     all.x=T)

finalResult <- within(finalResult,{
  symbolJoinTest[is.na(symbolJoinTest) | alpha==0.9 |
                 (symbolJoinTest=="*"&PRS_value!="High")] <- ""
})

# PRS to factor
finalResult$PRS <- factor(finalResult$PRS,levels=riskVars)
names(finalResult)[names(finalResult)=="DX.bl"] <- "Group"

# PANEL PLOT ······························································
# Auxiliar
  # separacio pointrange predisposition
  pos_dodge <- position_dodge(width=0.4) 

  # Etiquetes panell
  volume.labs <- volVars
  names(volume.labs) <- volVars
  volume.labs[c(7:11,13)] <- c("Hippocampal\nfissure",
                              "Hippocampal\ntail",
                              "Molecular\nlayer hp",
                              "Para-\nsubiculum",
                              "Pre-\nsubiculum",
                              "Whole\nhippocampus")

  PRS.labs <- sapply(strsplit(riskVars,split="_"),"[",1)
  names(PRS.labs) <- riskVars
  PRS.labs[PRS.labs=="LONGEVITYnoAPOE"] <- "LONGEVITY\nno APOE"
  PRS.labs[PRS.labs=="ADnoAPOE"] <- "AD\nno APOE"
  
  finalResult$Group <- factor(finalResult$Group,levels = c("AD","MCI","CN"))

#PLOT
panel <- 
  ggplot(finalResult,aes(x=Group)) +
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  geom_pointrange(aes(y=time_visit.trend , ymin=lower.CL , ymax=upper.CL, 
                      color=as.factor(PRS_value), alpha=alpha),
                  position = pos_dodge, show.legend = F,size=.3,
                  linewidth=.5) +
  geom_text(aes(y=upper.CL,label=symbolJoinTest))+
  facet_grid(rows=vars(PRS),cols=vars(Volumes),scales = "free",
             labeller=labeller(Volumes=volume.labs, PRS=PRS.labs))+
  scale_color_manual(values = c("#D95F02","#7570B3")) + 
  scale_alpha(range = c(0.3,1))+
  coord_flip()+
  
  theme_bw() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90,hjust=1)) +
  ggtitle("PRS effect over time rate of change of hippocampal volumes by diagnostic group.") +
  labs(subtitle = "Adjusted for Age at baseline, sex, intracranial volume and years of education. <br/>
       Trends for <span style = 'color: #D95F02;'>**Low-Intermediate**</span> and <span style = 'color: #7570B3;'>**High**</span> predisposition to PRS trait or disease." ,
       caption="Linear coefficient and 95% confidence intervals for volumetric rate of change over time (mm3 change per year),by diagnostic group and level of predisposition for the PRS trait/disease. Highlighted<br/>
       estimates indicate significance differences between High and Low-Intermediate predispositions for a particular diagnostic group. Panels with an \\* present significant difference in PRS effect on<br/>
       volume change between groups. The red dashed line marks no time trend. Estimates correspond to models L.2") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_markdown(size = 16, hjust = 0))

panel

plotData <- subset(finalResult, 
                   PRS %in% paste(c("AD","ADnoAPOE","LONGEVITY","LONGEVITYnoAPOE"),
                                  "predisposition",sep="_")) 
panel %+% plotData





## Visualització dels pendents
library(ggeffects)
plotModels <- function(model){
  
  fit <- model$mod_3wayIntDX
  prs <- attributes(fit)$PRS
  volume <- attributes(fit)$volume
  
  # Plots
  m <- ggpredict(fit,terms=c("time_visit",prs,"DX.bl"),
                 condition = c(ICV=mean(dd2$ICV)))
  
  p <- plot(m,add.data = T,ci.style = "dash",dot.alpha = 0.2,dot.size = 1.5,
            colors = c("#D95F02","#7570B3"),line.size = 1,use.theme = T)
  
  # Save
  path <-here("results",paramT,"2_Modelling","2_longitudinalTrajectories",
              "trajectories","time_visit")
  filename <- paste0("modelTimeVisit_",volume,"_",prs,".png")
  ggsave(plot=p,filename = filename,path=path,device = "png")

  # Return
  return(plot=p)
}

plots<- lapply(result,FUN=plotModels)


#-------------------------------------------------------------------------------
# B) OBJECTIU 2.bis
#    Efecte del PRS sobre la reducció dels volums, depèn de l'estat diagnòsitc?
#    Efecte PRS sobre pendent beta_time (coeficient beta_timeVisit:PRS) 
#    estratificant per grup diagnòsitc.
#-------------------------------------------------------------------------------
# volume <- "Hippocampal_tail" ; PRS <- "AD_predisposition" ; group <- "All"
# data <- dd2

#--------------------------------------------------------------------------------
# FUNCIÓ prsTimeEffect_strata(): per cada volum, PRS continu s'ajusta el model:
#   
#   vol ~ time_visit*PRS + AGE + ICV/1000 + PTGENDER + PTEDUCAT + randomEffects.
#   
#  globalment ("all") o estratificant per grup diagnòstic.
#  Extraure el terme fix interacció time_visit*PRS (increment del rate of change
#  associat a increment unitari del PRS), amb interval de confiança i pvalor.
#  Els graus de llibertat aproximats es calculen amb el mètode de Satterthwaite
#  a fi de calcular pvalors.
#  Els intervalos de confiança es calculen
#--------------------------------------------------------------------------------


prsTimeEffect_strat <- function(volume,PRS,group = "All",data){
  data_strat <- data
  
  # i) Formula ··················································
  predictors <- paste("time_visit",c("DX.bl",PRS),sep=" * ")
  covars <- c("AGE","I(ICV/1000)","PTGENDER","PTEDUCAT")
  random <- "(1+time_visit|PTID)"
  
  # Estratificar per grup diagnòstic
  if(group %in% c("CN","MCI","AD")){
    data_strat <- subset(data,DX.bl == group)
    predictors <- paste("time_visit",PRS,sep=" * ")
  }
  
  # Estratificar per sexe
  if(group %in% c("Female","Male")){
    data_strat <- subset(data,PTGENDER == group)
    covars <-covars[!(covars=="PTGENDER")]
  }
  
  
  form <- paste(c(predictors,covars,random),collapse= " + ")
  form <- paste(volume,form, sep = " ~ ")
  
  # ii) Model amb interacció time_visit*PRS ························
  mod <- suppressMessages(
    lme4::lmer(formula = form, data=data_strat,
               control = lmerControl(optimizer ="Nelder_Mead")))
  
  # Si la matriu de covariancies dels factors aleatoris es singular,
  # reajustar el model sense random slope pel temps segons PTID...
  if(isSingular(mod)){
    cat("Singular random effects covariance matrix for ",volume," and ",PRS," in ",group,"strata.
Removing random slope to simplify model...\n\n")
    random <- "(1|PTID)"
    form <- paste(c(predictors,covars,random),collapse= " + ")
    form <- paste(volume,form, sep = " ~ ")
    mod <- lme4::lmer(formula = form, data=data_strat,
                      control = lmerControl(optimizer ="Nelder_Mead"))
  }
  
  # iii) Efectes fixos del model resultant ·······························
  fixef <- extract_fixed_effects(mod)
  
  
  # iv) Extraure efecte PRS sobre variacó volum en el temps ···············
  effect_pos <- which(fixef$term == paste0("time_visit:",PRS,"High"))
  
  # v) Resultat en un dataframe: valor, 95%CI i pvalor test wald
  # amb graus de llibertat per Satterthwaite ······························
  res <- as.data.frame(fixef[effect_pos,c("value","lower_2.5",
                                          "upper_97.5","p_value")])
  
  rownames(res) <- paste(volume,PRS,sep="_")
  
  # Retornar resultats: efecte PRS*time_visit
  return(list("PRS_time"=res,"model"=mod))
}
#--------------------------------------------------------------------------


### i) Matriu de combinacions PRS, volums hipocamp i estrat segons DX.bl
modGridDx <- expand.grid(PRS = riskVars,
                         Volume = sort(volVars),
                         DX.bl = levels(dd2$DX.bl),
                         stringsAsFactors = F)

### ii) Apliquem la funció effectPRS_strat() sobre la GRID 
result_stratDx <- mapply(FUN = prsTimeEffect_strat,
                         volume = modGridDx$Volume,
                         PRS = modGridDx$PRS,
                         group = modGridDx$DX.bl,
                         MoreArgs = list(data=dd2), 
                         SIMPLIFY = F,USE.NAMES = F)

# i arreglem el resultat...
res_prsTimeEffect_stratDx <- lapply(result_stratDx,"[[","PRS_time")
res_prsTimeEffect_stratDx <- do.call(rbind,res_prsTimeEffect_stratDx) # consolidar en dataFrame
res_prsTimeEffect_stratDx <- cbind(modGridDx,res_prsTimeEffect_stratDx) # Afegir PRS i volume 
rownames(res_prsTimeEffect_stratDx) <- NULL # no rownames
res_prsTimeEffect_stratDx$PRS <- factor(res_prsTimeEffect_stratDx$PRS,levels=riskVars) 
res_prsTimeEffect_stratDx$DX.bl <- factor(res_prsTimeEffect_stratDx$DX.bl,
                                          c("CN","MCI","AD"))
head(res_prsTimeEffect_stratDx)

### iii) Corregim pvalors segons FDR, considerant cada volum independent,
#       i corregint per 6 comparacions (6 PRS) o models per cada un dels volums
#       i estrat segons grup diagnòsitc
aux <- lapply(split(res_prsTimeEffect_stratDx,
                    list(res_prsTimeEffect_stratDx$Volume,
                         res_prsTimeEffect_stratDx$DX.bl)),
              FUN=function(x)p.adjust(x$p_value,method="fdr"))
res_prsTimeEffect_stratDx$pvalue_corrected <- do.call(c,aux)


### iv) Afegim simbols de significació
pvalueVars <- names(res_prsTimeEffect_stratDx)[7:8]
pvalueSymbol <- paste(pvalueVars,"symbol",sep="_")

res_prsTimeEffect_stratDx[pvalueSymbol] <-
  lapply(res_prsTimeEffect_stratDx[pvalueVars], 
         FUN=cut,
         breaks = c(-Inf,0.001,0.01,0.05,0.1,1),
         labels=c("***","**","*","·",""))


# RESULTAT
head(res_prsTimeEffect_stratDx)
(significant_effectsDx <- subset(res_prsTimeEffect_stratDx, 
                                 pvalue_corrected_symbol!=""))


### PLOT EFECTES time_visit*PRS per cada estrat ······························
plotDataDx <- subset(res_prsTimeEffect_stratDx,
                     Volume!="Whole_hippocampus")
names(plotDataDx)[names(plotDataDx)=="DX.bl"] <- "Group"

heatMapDx <- 
  ggplot(data=plotDataDx, aes(x=Volume, y=PRS, fill=value)) + 
  geom_tile() +
  geom_text(aes(label=pvalue_corrected_symbol),size=4)+
  scale_fill_gradient2(low="red4",mid = "white", high="#0066FF", midpoint = 0)+
  facet_wrap(~Group,ncol=2)+
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank(),
        axis.text = element_text(face="bold",size = 8),
        axis.text.x = element_text(angle=90,vjust = 0.6,hjust = 1)) +
  
  ggtitle("Corrected significant effect of investigated PRS on hippocampal\nsubfield volume's rate of change") +
  labs(subtitle="Stratified by diagnostic status",
       caption="Linear coefficient for interaction between time and PRS. Effects adjusted for Age at baseline (years), intracranial volume(cm3), sex and years<br/>
       of education. Model allows a diferent volume rate of change depending on diagnositc status. Also, intercept and rate of change in time<br/>
  are included as random effects. Significant interactions are reported as Wald-test pvalues at levels 0.1(·), 0.05(\\*), 0.01(\\**) and 0.001(\\***),<br/>
  corrected for multiple comparisons by FDR method. Non-significant differences are not shown.<br/>
  **Source**: ADNI project | **Units** Average increment in volume rate of change (mm3) associated with having high predisposition for the PRS trait.") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = 0))

heatMapDx






#-------------------------------------------------------------------------------
# B) OBJECTIU 3
#    Efecte del PRS sobre la reducció dels volums, depèn del sexe?
#-------------------------------------------------------------------------------
# volumes <- "Hippocampal_tail"
# PRS <- "AD_predisposition"
# data <- dd2

# Interacció time_visit*PRS*DX.bl significativa?
testSx <- lapply(result,"[[","AnovaTest_3wayIntSx")
testSx <- lapply(testSx,"[[",'Pr(>Chisq)')
testSx <- do.call("c",testSx)
res_prsTimeEffect$pval_3wayIntSx <- round(testSx,3)

testSx[testSx<0.05]
subset(res_prsTimeEffect,pval_3wayIntSx<0.05,c(1:3,8:10))



### CALCULAR ELS TRENDS LINEALS SEGONS predisposition low or high ···················
# Debbuging:x<-result$`Hippocampal_tail:AD_predisposition` # per provar


trendsSx <- lapply(result,function(x){
  mod <- x$mod_3wayIntSx
  PRS <- attributes(mod)$PRS
  specs <- as.formula(paste0("~",PRS,"*PTGENDER"))
  weights <- "prop"
  
  trend <- emtrends(mod, specs = specs, "time_visit",
                    weights = weights)
  
  return(trend)
})

# A dataframe i canviar nom columna PRS a genèric per despres fer el merge
trendsSx.df <- lapply(trendsSx,summary)

trendsSx.df <- lapply(trendsSx.df,function(x){
  names(x)[1] <- "PRS_value"
  return(x)
})

trendsSx.df <- do.call(rbind,trendsSx.df)


### CALCULAR CONTRASTOS High-Low ···················
contrastPRS_Sx <- lapply(trendsSx,function(x){
  n <- names(x@levels)[1]
  summary(contrast(x,"consec",simple=n))
})

contrastPRS_Sx <- as.data.frame(do.call(rbind,contrastPRS_Sx))

res_contrastPRS_Sx <- data.frame(Volume = rep(volVars,each=2*7),
                                PRS = rep(rep(riskVars,each=2),times=13))

res_contrastPRS_Sx <- cbind(res_contrastPRS_Sx,contrastPRS_Sx)

# Corregir per comparacions múltiples FDR
aux <- lapply(split(res_contrastPRS_Sx,res_contrastPRS_Sx$Volume),
              FUN=function(x)p.adjust(x$p.value,method="fdr"))

res_contrastPRS_Sx$pvalue_corrected <- do.call("c",aux)

# Afegir simbols de significació
res_contrastPRS_Sx$pvalue_corrected_symbol <-
  cut(res_contrastPRS_Sx$pvalue_corrected,
      breaks = c(-Inf,0.001,0.01,0.05,0.1,1),
      labels=c("***","**","*","·",""))


(res_contrastPRS_Sx_signif <- subset(res_contrastPRS_Sx,pvalue_corrected<0.05))

prsDiff_Sx <- res_contrastPRS_Sx_signif[c("Volume","PRS")]
prsDiff_Sx <- with(prsDiff_Sx,unique(paste(Volume,PRS,sep=":")))


# Exportar taula significatius formatejada
res_contrastPRS_Sx_signif[5:8] <- lapply(res_contrastPRS_Sx_signif[5:8],round,2)
res_contrastPRS_Sx_signif[9:10] <- 
  lapply(res_contrastPRS_Sx_signif[9:10],function(x){
    ifelse(x<0.001,"<0.001",sprintf("%.3f",x))
})

res_contrastPRS_Sx_signif$pvalue_corrected <- 
  with(res_contrastPRS_Sx_signif,paste0(pvalue_corrected,"(",pvalue_corrected_symbol,")"))

rownames(res_contrastPRS_Sx_signif) <- NULL
res_contrastPRS_Sx_signif$pvalue_corrected_symbol <- NULL

tab <- kbl(res_contrastPRS_Sx_signif,align = "lllccccccc",row.names = F,
           caption = "Significant effects of PRS on time trend. Models L.3, linear contrasts 1")


tab <- kable_classic(tab,full_width=F,font="Cambria")
formattedTab <- footnote(tab,
                         general = c("Wald test for $H_0: \\beta5_{time:High} =0$ when female, and $H_0: \\beta5_{time:High} + \\beta8_{time:High:Male}=0$ when male.",
                                     "Pvalues were corrected for multiple comparisons using FDR method. Significance level: (*)<.05, (**)<.01, (***)<.001"),
                         title_format = c("bold"))

formattedTab_contrast1_Sx <- row_spec(formattedTab,0,bold=T)
formattedTab_contrast1_Sx


# Exportar excel contrastos 1
  tab <- res_contrastPRS_Sx
  rownames(tab) <- NULL
  tab <- split(tab,tab$Volume)

  # Build Excel worksheet
  wb <- createWorkbook()
  invisible(lapply(names(tab),function(x)addWorksheet(wb,x)))
  invisible(
    lapply(1:length(tab),
           function(i) writeDataTable(wb = wb,
                                      x = tab[[i]],
                                      sheet = names(tab)[i],
                                      withFilter = T,
                                      tableStyle = "TableStyleLight1")
    )
  )

  path <- here("results",paramT,"2_Modelling","2_longitudinalTrajectories",
             "taules","modelL3_contrast1.xlsx")
  saveWorkbook(wb,file=path,overwrite = T)


# Diferències entre grups quant a efecte genètic?

aux <- lapply(prsDiff_Sx,function(x){
  
  test <- summary(contrast(trendsSx[[x]],interaction=c("consec","consec")))
  mod <- result[[x]]$mod_3wayIntSx
  
  res <- data.frame(Volumes = attributes(mod)$volume,
                    PRS = attributes(mod)$PRS,
                    t = test$t.ratio,
                    df = test$df,
                    Pval = test$p.value)
  res$symbolJoinTest <- ifelse(res$Pval<0.05,"*","")
  return(res)
})
rSx <- do.call(rbind,aux)


# Exportar taula contrast 2 formatejada

tab <- rSx
tab[c("t","df")] <- lapply(tab[c("t","df")], round, 2)
tab$Pval <- ifelse(tab$Pval<0.001,"<0.001",sprintf("%.3f",tab$Pval))
tab$Pval <- with(tab,ifelse(symbolJoinTest=="*",paste0(Pval,"(*)"),Pval))
tab$symbolJoinTest <- NULL

tab <- kbl(tab,align = "llccc",row.names = F,
           caption = "Differences between sex groups regarding PRS effects on volume rate of change. Models L3, linear contrasts 2")


tab <- kable_classic(tab,full_width=T,font="Cambria")
formattedTab <- footnote(tab,
                         general = c("Wald Test for linear contrast. $H_0:\\beta8_{time:Male:High}=0$",
                                     "Significance level: (*)<.05"),
                         title_format = c("bold"))

formattedTab_contrast2_Sx <- row_spec(formattedTab,0,bold=T)
formattedTab_contrast2_Sx



### DATAFRAME DE RESULTATS PER PLOTEJAR ···································

finalResultSx <- data.frame(Volumes = rep(volVars,each = 4*length(riskVars)),
                          PRS = rep(riskVars,each=4))
finalResultSx <- cbind(finalResultSx,trendsSx.df)

# Afegim significació corregida (el mateix pvalor pel low i high en cada
# volum-PRS-DX.bl)
finalResultSx$pvalue <- rep(res_contrastPRS_Sx$p.value,each=2)
finalResultSx$pvalue_corrected <- rep(res_contrastPRS_Sx$pvalue_corrected,each=2)


# Transparency per significació del contrast PRS
finalResultSx$alpha <- with(finalResultSx,ifelse(pvalue_corrected<0.05,1,.9))


# JoinTest
finalResultSx <- merge(finalResultSx,
                       rSx[c("Volumes","PRS","symbolJoinTest")],
                       by=c("Volumes","PRS"),
                     all.x=T)

finalResultSx <- within(finalResultSx,{
  symbolJoinTest[is.na(symbolJoinTest) | alpha==0.9 |
                   (symbolJoinTest=="*"&PRS_value!="High")] <- ""
})

# PRS to factor
finalResultSx$PRS <- factor(finalResultSx$PRS,levels=riskVars)
names(finalResultSx)[names(finalResultSx)=="PTGENDER"] <- "Group"
finalResultSx$Group <- factor(finalResultSx$Group,levels = c("Male","Female"))

panel %+% finalResultSx +
  ggtitle("PRS effect over time rate of change of hippocampal volumes by sex group.") +
  labs(subtitle = "Adjusted for Age at baseline, diagnostic group, intracranial volume and years of education. <br/>
       Trends for <span style = 'color: #D95F02;'>**Low-Intermediate**</span> and <span style = 'color: #7570B3;'>**High**</span> predisposition to PRS trait or disease." ,
       caption="Linear coefficient and 95% confidence intervals for volumetric rate of change over time (mm3 change per year),by sex group and level of predisposition for the PRS trait/disease. Highlighted<br/>
       estimates indicate significance differences between High and Low-Intermediate predispositions for a particular sex group. Panels with an \\* present significant difference in PRS effect on<br/>
       volume change male and female. The red dashed line marks no time trend. Estimates correspond to models L.3")
  

# EXPORT EXCEL TIME TRENDS MODELS L2 i L3

# Get table for parameter estimates

timeTrendList <- list(finalResult,finalResultSx)
names(timeTrendList) <- paste0("Models L",2:3)
timeTrendList <- lapply(timeTrendList,function(x){
  res <- x
  res$model <- with(x,paste(Volumes,PRS,sep=":"))
  res <- res[c("model",names(x)[1:11])]
  return(res)
})

# Build Excel worksheet
wb <- createWorkbook()
invisible(lapply(names(timeTrendList),function(x)addWorksheet(wb,x)))
invisible(
  lapply(1:length(timeTrendList),
         function(i) writeDataTable(wb = wb,
                                    x = timeTrendList[[i]],
                                    sheet = names(timeTrendList)[i],
                                    withFilter = T,
                                    tableStyle = "TableStyleLight1")
  )
)




path <- here("results",paramT,"2_Modelling","2_longitudinalTrajectories",
             "taules","longitudinalModels_timeTrends.xlsx")

saveWorkbook(wb,file=path,overwrite = T)





# EXPROT LONGITUDINAL MODEL PARAMETERS

  # Get table for parameter estimates

  modelParams <- lapply(c("modGlobal","mod_3wayIntDX","mod_3wayIntSx"),
                        FUN = ensambleParamsTable, data=result)
  
  names(modelParams) <- paste0("Models L",1:3)
  
  
  # Build Excel worksheet
  wb <- createWorkbook()
  invisible(lapply(names(modelParams),function(x)addWorksheet(wb,x)))
  invisible(
    lapply(1:length(modelParams),
           function(i) writeDataTable(wb = wb,
                                 x = modelParams[[i]],
                                 sheet = names(modelParams)[i],
                                 withFilter = T,
                                 tableStyle = "TableStyleLight1")
    )
  )

  
  
  
path <- here("results",paramT,"2_Modelling","2_longitudinalTrajectories",
             "taules","longitudinalModels_parameters.xlsx")

saveWorkbook(wb,file=path,overwrite = T)

