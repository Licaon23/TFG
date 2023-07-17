###############################################################################
### MODELLING                                                                 #
###   Cross sectional modeling: hippocampal volumes                           #
###   Version 1                                                               #
###   March   2023                                                            #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###                                                                           #
###############################################################################
source("Rscripts/libraries.R")
source(here(RscriptsDir,"3_modelling","modelling_auxFunctions.R"))

# load(here(cleandataDir,"dataForAnalysis_3T.RData"))
# paramT <- "Dataset_3T"

load(here(cleandataDir,"dataForAnalysis_1T.RData"))
paramT <- "Dataset_1T"
#===============================================================================
# 1) Establir diferències entre grups de diagnòstic o sexe biològic en els
#     volums hipocampals a mesura basal, ajustant per edat, sexe/grup,
#     i volum intracraneal total
#===============================================================================
volVars <- sort(names(basal)[17:29])
PRSVars <- names(basal)[9:15]
RiskVars <- names(basal)[30:36]

#------------------------------------------------------------------------------
# FUNCIO adjustedGroupDifferences 
#   INPUT:
#     - volume: variable volum hippocamp
#
#   TASQUES
#     - A) Ajustar model lineal per explicar volum hipocampal a instant basal
#           segons AGE, ICV, PTEDUCAT i amb interacció de DX.bl i PTGENDER.
#     - B) Extraure significació diferències en volum segons grup diagnòstic,
#           resultat a objecte 'res_signif'.
#     - C) Calcular mitjanes marginals, segons DX.bl i PTGENDER, fixant
#         variables AGE,ICV i PTEDUCAT a la mitjana
#     - D) Pairwise contrasts: diferències entre PTGENDER dins de cada grup
#           diagnòsitc. Resultat a objecte 'contrast'.
#     - E) Presentar resultats mitjanes marginals volum, a objecte 'res_emmeans'
#     - F) Fer gràfic diagnòstic de residus i guardar-los a directori
#
#   OUTPUT
#     - Llista amb els 3 objectes resultat: res_signif, contrast, res_emmeans
#------------------------------------------------------------------------------
adjustedGroupDifferences <- function(volume){
  
  # Model lineal amb interacció sexe-diagnostic status
  expr <- paste(volume,"AGE + ICV + PTEDUCAT + PTGENDER * DX.bl",
                sep = " ~ ")
  mod <- lm(formula = expr, data = basal)
  smod <- summary(mod)
  anovaTab <- anova(mod)

  # Significació diferències en volum segons grup diagnòstic
  pvals_DXbl <- c(anovaTab["DX.bl","Pr(>F)"])
  pvals_DXbl_formatted <- ifelse(pvals_DXbl<0.001,"<0.001",
                                 sprintf("%.3f",pvals_DXbl))
  
  res_signif <- data.frame(Volume = volume,
                           Group = "DX.bl",
                           pval = pvals_DXbl_formatted,
                           test = "FTest")
  
  
  # Marginal means
  RG <- ref_grid(mod)
  
  emm_conjunt <- emmeans(RG,spec=~PTGENDER|DX.bl) # with interaction
  emm_conjuntDF <- as.data.frame(emm_conjunt)
  
  # Pairwaise contrasts: for each DX.bl level, are the sex difference 
  # statisticaly significant?
  contrast <- as.data.frame(pairs(emm_conjunt))
  contrast <- data.frame(Volume=rep(volume,3),
                         DX.bl = contrast$DX.bl,
                         group1 = rep("Female",3),
                         group2 = rep("Male",3),
                         statistic = contrast$t.ratio,
                         df = contrast$df,
                         p.adj = contrast$p.value)
  
  
  
  # Result dataframe
  res_emmeans <- rbind(expand.grid(levels(basal$PTGENDER),levels(basal$DX.bl)))
  res_emmeans <- cbind(rep(volume,nrow(res_emmeans)),
                       res_emmeans)
  
  res_emmeans[c("emmean","lwr","upr")] <- 
    emm_conjuntDF[,c("emmean","lower.CL","upper.CL")]

  res_emmeans$value <- with(res_emmeans,
                            sprintf("%.f (%.f,%.f)",emmean,lwr,upr))
  
  colnames(res_emmeans)[1:3] <- c("Volume","PTGENDER","DX.bl")
  
  
  
  
  # Diagnostic residual plots
  file <- paste(volume,"diagnosticPlots.png",sep="_")
  path <- here("results",paramT,"2_Modelling",
               "1_baselineVolumes_groupAdjustedDifferences",
               "residualDiagnosticPlots")
  
  if(!(any(grepl(file,list.files(path))))){
    resid_panel(mod,qqbands = T,plots = "all")
    ggsave(file=file, path=path ,width = 8, height = 6)
  }
  
  # Return result
  return(list(significance = res_signif , 
              marginal_means = res_emmeans,
              contrast = contrast))
}
#------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# a) Aplicar la funció adjustedGroupDifferences()  a tots els volums
#-----------------------------------------------------------------------------
res <- lapply(volVars, FUN=adjustedGroupDifferences)

# Extraure objectes: mitjanes marginals, significació entre grups diagnòsic
#   i diferències entre sexes per cada grup diagnòsitc
res_significance <- do.call(rbind,lapply(res,"[[",1))
res_emmeans <- do.call(rbind,lapply(res,"[[",2))
res_contrast <- do.call(rbind,lapply(res,"[[",3))


#-----------------------------------------------------------------------------
# b) Generar taula de resultats: mitjanes marginals dels volums, per DX.bl
#     i PTGENDER, fixats els valors per AGE,ICV i PTEDUCAT
#-----------------------------------------------------------------------------

# Generar taula de resultats
emmeansData <- subset(res_emmeans,select=c(1,2,3,7))
emmeansData$class <- with(emmeansData, paste(DX.bl,PTGENDER,sep="_"))
emmeansData <- emmeansData[,c("Volume","class","value")]

# Contrastos diferències entre sexes per cada grup diagnòstic--------------
cont <- split(res_contrast,res_contrast$DX.bl)
cont <- lapply(cont, "[[","p.adj")
cont <- lapply(cont,function(x)sprintf("%.3f",x))

# Contruir taula------------------------------------
taula <- reshape(emmeansData, direction = "wide", idvar = "Volume",
                 timevar = "class")
taula[paste(c("CN","MCI","AD"),"sexContrast",sep="_")] <- cont

taulaDef <- taula[c(1:3,8,4,5,9,6,7,10)]
pvalData <- subset(res_significance,Group=="DX.bl",select=c(Volume,pval))
taulaDef <- merge(taulaDef,pvalData,by="Volume")

#Sample size: n per cada sexe i gup diagnòsitc
x <- with(basal,table(PTGENDER,DX.bl))
n <- sprintf("N = %.i",x)

# Format taula HTLM
tab <- kbl(taulaDef,align = "lcccccccccr",
           col.names = c("",n[1:2],"",n[3:4],"",n[5:6],"",""),
           caption = "Adjusted marginal means for hippocampal volumes at baseline
           by diagnostic group") 

tab <- add_header_above(tab,c("",
                              rep(c(levels(basal$PTGENDER),"P"),3),"Pval_groups"),
                        line = F,bold=T)

tab <- add_header_above(tab,c(" " = 1,"CN"=3,"MCI"=3,"AD"=3,""),
                        bold = T,line = T,)

tab <- column_spec(tab,2:10,width = "2in")
tab <- kable_classic(tab,full_width=F,font="Cambria")

formattedTab <- footnote(tab,
                         general = c("Results are reported as estimated marginal mean volume and 95% confidence interval,adjusted by Age, Intracraneal total volume and sex.",
                                     "Differences between diagnostic groups were tested using F-test for global significance. P-values are reported as Pval_groups.",
                                     "Differences between sex within each diagnostic group were tested pairwise t-test. P-values are reported as P in every diagnostic group",
                                     "Diagnostic groups: Controls (CN), Mild Cognitive Impairment (MCI) and Alzheimer's Disease (AD) "),
                         title_format = c("bold")) 


# Guardar i visualitzar taula de resultats
file <- "result_emmeansVolume.png"
path <- here("results",paramT,"2_Modelling",
             "1_baselineVolumes_groupAdjustedDifferences")

if(!(any(grepl(file,list.files(path))))){
  save_kable(formattedTab,file=here(path,file),zoom=2)
}

# RESULTAT
formattedTab

#-----------------------------------------------------------------------------
# c) Gràfic: volums hippocampals segons DX.bl i PTGENDER, ajustats
#       per AGE, ICV i PTEDUCAT. 
#       Mitjanes marginals calculades a valor mitjà de AGE, valor mitjà
#       de ICV i mediana PTEDUCAT
#-----------------------------------------------------------------------------

# Dades necessàires
emmeansData <- subset(res_emmeans,Volume != "Whole_hippocampus",
                      select=-value)

# Simbols significació diferències entre sexes
res_contrast$p.adj.signif <- cut(res_contrast$p.adj,
                                 breaks=c(-Inf,0.001,0.01,0.05,0.1,1),
                                 labels=c("***","**","*","·",""))

# Afegir info contrastos a la taula emmeansData
emmeansData <- merge(emmeansData,
                     res_contrast[,c("Volume","DX.bl","group1","p.adj.signif")],
                     by.x=c("Volume","DX.bl","PTGENDER"),
                     by.y=c("Volume","DX.bl","group1"),all.x=T)

emmeansData$p.adj.signif[is.na(emmeansData$p.adj.signif)] <- ""


# Gràfic de mitjanes i interval de confiança pels volums
pos_dodge <- position_dodge(width=0.2)

ggplot(emmeansData, aes(x=DX.bl)) + 
  geom_pointrange(aes(y=emmean,ymin=lwr,ymax=upr,color=PTGENDER),
                  position = pos_dodge,show.legend = FALSE,size=0.3) +
  geom_text(aes(y=upr,label=p.adj.signif),vjust=-1) +
  scale_y_continuous(expand = expansion(mult = c(0,0.2))) + 
  scale_color_manual(values = c("#D95F02","#7570B3")) +
  facet_wrap(~Volume,scales = "free_y") + 
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank()) +
  
  ggtitle("Estimated volume means at baseline by biological sex and diagnostic group.") +
  labs(subtitle = "Adjusted for Age at baseline, sex, intracranial volume, diagnostic status and years of education. <br/>
       <span style = 'color: #D95F02;'>**Female**</span> and <span style = 'color: #7570B3;'>**Male**</span>" ,
       caption="**Source**: ADNI project | **Units** Volume measurements in mm3 | Differences between sex are reported as significance level of t-test sliced by diagnostic group, <br/>at levels 0.1, 0.05, 0.01 and 0.001. Non-significant differences are not shown.") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = -.05))

ggsave(file="emmeansAdjustedVolumes.png",
       path = here("results",paramT,"2_Modelling",
                   "1_baselineVolumes_groupAdjustedDifferences"),
       width = 25,height = 15,units = "cm")







#==============================================================================
# 2) ASSOCIACIONS ENTRE PRS I VOLUMS HIPOCAMP
#==============================================================================

#------------------------------------------------------------------------------
# FUNCIO effectPRS_strat 
#   INPUT:
#     - volume: variable volum hippocamp
#     - PRS
#     - group: "all", sense estratificar; o "CN","MCI","AD si estratifiquem
#         segons grup diagnòstic
#
#   TASQUES
#     - A) Subset de les dades segons grup, si escau
#     - B) Ajusta el model lineal additiu amb covariables i PRS
#     - C) Recupera efecte PRS sobre volum, interval confiança i pvalor
#
#   OUTPUT
#     - Dataframe amb 4 variables: Effect, lower.ci, upper.ci i Significance
#       corresponent a l'efecte del PRS sobre volume per group "group".
#------------------------------------------------------------------------------

# Per debbuing: volume <- "CA1";PRS<-"PRS_AD";group <- "All"
effectPRS_strat <- function(volume,PRS,group="All"){
  data <- basal
  covars <- c("AGE","ICV","PTEDUCAT","PTGENDER","DX.bl")

  # Treballar amb tots els grups diagnòstic o estratificar
    if(group!="All"){
    data <- subset(basal,DX.bl == group)
    covars <- covars[covars!="DX.bl"]
    }
  
  # Formula per model
  covars <- paste(covars,collapse = " + ")
  expr <- paste(volume,
                paste(c(covars,PRS), collapse=" + "),
                sep = " ~ ")
  
  # Model: fit, summary i intervals confiança
  mod <- lm(formula = expr, data = data)        
  smod <- summary(mod)$coefficients
  CIs <- confint.default(mod)
  df <- summary(mod)$df[2]
  
  # Recuperar coeficient associat al PRS i el seu interval de confiança
  pos_beta <- grep(PRS,rownames(smod))
  beta_PRS <- smod[pos_beta,]
  CI_PRS <- CIs[pos_beta,,drop=T]
  res <- data.frame(Effect = beta_PRS[1],
                    SE = beta_PRS[2],
                    df = df,
                    t.ratio = beta_PRS[3],
                    lower.ci = CI_PRS[1],
                    upper.ci = CI_PRS[2],
                    Significance = beta_PRS[4])

  
  
  rownames(res) <- expr
  attributes(mod)$PRS <- PRS
  attributes(mod)$volume <- volume
  return(list("PRSEffect" = res,
              "ModelGlobal" = mod))
}
#------------------------------------------------------------------------------


#---------------------------------------------------------------------------
# a) Aplicar la funció effectPRS_strat() sobre la matriu modGrid, combinacions
#     de volume i pRS
#---------------------------------------------------------------------------

# Inicialitzem la matriu de combinacions
modGrid <- expand.grid(PRS = PRSVars,
                       Volume = sort(volVars),
                       DX.bl = c("All",levels(basal$DX.bl)),
                       stringsAsFactors = F)

# Apliquem la funció effectPRS_strat() sobre la GRID i arreglem el resultat
auxQ1 <- mapply(FUN = effectPRS_strat,
                volume = modGrid$Volume,
                PRS = modGrid$PRS,
                group = modGrid$DX.bl,
                SIMPLIFY = F,USE.NAMES = F)

result_strat <- lapply(auxQ1,"[[","PRSEffect")
result_strat <- do.call(rbind,result_strat)

modGrid <- cbind(modGrid,result_strat)
rownames(modGrid) <- NULL

modGrid$PRS <- factor(modGrid$PRS,levels=PRSVars)
modGrid$DX.bl <- factor(modGrid$DX.bl, levels=c("All","CN","MCI","AD"))

  # Resultat: modGrid

# Corregim pvalors segons FDR, considerant cada volum independent,
#  i corregint per 6 comparacions (6 PRS) o models per cada un dels volums
aux <- lapply(split(modGrid,list(modGrid$Volume,modGrid$DX.bl)),
              FUN=function(x)p.adjust(x$Significance,method="fdr"))
modGrid$Significance_corrected <- do.call(c,aux)

# Afegim simbols de significació
modGrid[paste(names(modGrid)[10:11],"symbol",sep="_")] <-
  lapply(modGrid[names(modGrid[10:11])], FUN=cut,
         breaks = c(-Inf,0.001,0.01,0.05,0.1,1),
         labels=c("***","**","*","·",""))

significant_effects <- subset(modGrid, Significance_corrected_symbol!="")



#-----------------------------------------------------------------------------
# b) Q1: existeix una associació entre predisposició genètica a neurodegeneració
#     i envelliment amb diferències en volums hipocampals?
#-----------------------------------------------------------------------------

# Mostrem els efetes de cada PRS i volum en un heatmap, sense estratificar
# i sense mostrar el volum hipocampal total (diferent escala de l'efecte)

plotData <- subset(modGrid,Volume!="Whole_hippocampus" & DX.bl=="All")
names(plotData)[3] <- "Group"

heatmapT1 <- ggplot(data=plotData,aes(x=Volume,y=PRS,fill=Effect)) + 
  geom_tile() +
  geom_text(aes(label=Significance_symbol),size=6)+
  scale_fill_gradient2(low="red4",mid = "white", high="#0066FF", midpoint = 0)+
  facet_wrap(~Group)+
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank(),
        axis.text = element_text(face="bold",size = 8),
        axis.text.x = element_text(angle=90,vjust = 0.6,hjust = 1)) +
  
  ggtitle("Nominal significant associations between investigated\nPRS and hippocampal subfield volumes") +
  labs(caption="PRS linear effects on hippocampal volumes, which corresponds to beta1 coefficient of T.1 models. Nominal pvalues are reported as \\** <.01, \\* <.05 and<br/>· <.1.
  After correcting for multiple comparisons, none of them remained significant. Effects adjusted for Age at baseline, intracranial volume, biologiacal<br/>
  sex, diagnostic status and years of education.<br/>
  **Source**: ADNI project") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = 0))


heatmapT1
subset(modGrid, Significance_symbol !="" & DX.bl=="All",
       c(Volume,PRS,Effect,Significance_symbol,Significance_corrected_symbol))


#-----------------------------------------------------------------------------
# b) Q2: l'efecte del PRS sobre les diferències en els volums hipocampals
#     al instant basal depèn de l'estat diagnòstic dels individus?
#     Contrast: tots els termes de la interacció PRS:DXbl son 0
#-----------------------------------------------------------------------------
#volume <- "Hippocampal_tail"; PRS <- "PRS_AD"; interact <- "DX.bl"

models_interact <- function(volume,PRS,interact){
  data <- basal
  covars <- c("AGE","ICV","PTEDUCAT","PTGENDER","DX.bl")
  
  # Formula per model
  covars <- paste(covars <- covars[covars!=interact],collapse = " + ")
  prs_interact <- paste(PRS,interact,sep = " * ")
  expr <- paste(volume,
                paste(c(covars,prs_interact), collapse=" + "),
                sep = " ~ ")
  
  # Model fit
  mod <- lm(formula = expr, data = data)
  attributes(mod)$PRS <- PRS
  attributes(mod)$volume <- volume
  
    
  # Beta effects
  specs <- as.formula(paste0("~",PRS,"|",interact))
  trend <- emtrends(mod, specs = specs, var=PRS, infer=c(F,T))
  
  attributes(trend)$PRS <- PRS
  attributes(trend)$volume <- volume
  
  return(list("PRSEffect" = trend,
              "modDXInteract" = mod))
}

modGridQ2 <- subset(modGrid,DX.bl=="All",select=c(PRS,Volume))

#x <- models_interact("Whole_hippocampus","PRS_AD",interact="DX.bl")

auxQ2 <- mapply(FUN = models_interact,
                      volume = as.character(modGridQ2$Volume),
                      PRS = as.character(modGridQ2$PRS),
                      MoreArgs = list(interact="DX.bl"),
                      SIMPLIFY = F)

betaPRS_DXInt <- lapply(auxQ2,"[[","PRSEffect")
names(betaPRS_DXInt) <- paste(modGridQ2$Volume,modGridQ2$PRS,
                              sep=" ~ ") 

# A dataframe i canviar nom columna PRS a genèric per despres fer el merge

aux <- lapply(betaPRS_DXInt, function(x){
  df <- summary(x)
  #names(df)[3] <- "PRS_effect"
  betas <- df[,c(2:7)]
  names(betas)[c(2,6)] <- c("Effect","Significance")
  
  # Result
  n <- nrow(betas)
  res <- data.frame(Volume = rep(attributes(x)$volume, times = n),
                    PRS = rep(attributes(x)$PRS, time = n))

  res <- cbind(res,betas)
  return(res)
})

resultQ2 <- do.call(rbind,aux)
resultQ2$PRS <- factor(resultQ2$PRS,levels=PRSVars)
resultQ2$DX.bl <- factor(resultQ2$DX.bl, levels=c("CN","MCI","AD"))

# Corregim pvalors segons FDR, considerant cada volum independent,
#  i corregint per 6 comparacions (6 PRS) o models per cada un dels volums
aux <- lapply(split(resultQ2,resultQ2$Volume),
              FUN=function(x)p.adjust(x$Significance,method="fdr"))
resultQ2$Significance_corrected <- do.call(c,aux)

# Afegim simbols de significació
resultQ2[paste(names(resultQ2)[8:9],"symbol",sep="_")] <-
  lapply(resultQ2[names(resultQ2[8:9])], FUN=cut,
         breaks = c(-Inf,0.001,0.01,0.05,0.1,1),
         labels=c("***","**","*","·",""))



plotData <- subset(resultQ2,Volume!="Whole_hippocampus")
names(plotData)[3] <- "Group"
heatMapT2 <- heatmapT1 %+% plotData +
  labs(caption="PRS linear effects on hippocampal volumes, which corresponds to beta1 or beta1+beta3 coefficients of T.2 models. Nominal pvalues are reported as \\* <.05 <br/>
  and · <.1. After correcting for multiple comparisons, none of them remained significant. Effects adjusted for Age at baseline, intracranial volume, biologiacal sex,<br/>
  and diagnostic status and years of education.<br/>
  **Source**: ADNI project") 
heatMapT2  
subset(resultQ2, Significance_symbol !="",
       c(Volume,PRS,DX.bl,Effect,Significance_symbol,Significance_corrected_symbol))



#-----------------------------------------------------------------------------
# c) Q3: quantificar la associació entre la predisposició genètica i el
#     volum hipocampal segons quin sigui l'estat diagnòsitc:
#     ANALISI ESTRATIFICAT
#-----------------------------------------------------------------------------

## HeatMap nominal stratified by diagnosis
plotData <- subset(modGrid,Volume!="Whole_hippocampus")
names(plotData)[3] <- "Group"
heatmapT1 %+%plotData +
  labs(caption="PRS linear effects on hippocampal volumes, which corresponds to beta1 coefficient of T.1 models. Nominal pvalues are reported as \\** <.01,<br/>
  \\* <.05 and · <.1. After correcting for multiple comparisons, none of them remained significant. Effects adjusted for Age at baseline,<br/>
  intracranial volume, biologiacal sex and years of education.<br/>
  **Source**: ADNI project")
  

# resultats iguals que amb interacció que estratificant


#-----------------------------------------------------------------------------
# d) Q4.1: l'associació entre PRS i volums, globalment per grup diagnòstic,
#     és diferent segons sexe?
#-----------------------------------------------------------------------------

# ## GLOBALMENT ·······················································
# modGridQ4 <- subset(modGrid,DX.bl=="All",select=c(PRS,Volume))
# 
# #x <- models_interact("Whole_hippocampus","PRS_AD",interact="DX.bl")
# 
# betaPRS_Sex <- mapply(FUN = models_interact,
#                       volume = as.character(modGridQ4$Volume),
#                       PRS = as.character(modGridQ4$PRS),
#                       MoreArgs = list(interact="PTGENDER"),
#                       SIMPLIFY = F)
# betaPRS_Sex <- lapply(betaPRS_Sex,"[[","PRSEffect")
# names(betaPRS_Sex) <- paste(modGridQ4$Volume,modGridQ4$PRS,
#                             sep=" ~ ") 
# 
# # A dataframe i canviar nom columna PRS a genèric per despres fer el merge
# 
# aux <- lapply(betaPRS_Sex, function(x){
#   df <- summary(x)
#   betas <- df[,c(2,3,7)]
#   names(betas)[c(2,3)] <- c("Effect","Significance")
#   
#   # Result
#   n <- nrow(betas)
#   res <- data.frame(Volume = rep(attributes(x)$volume, times = n),
#                     PRS = rep(attributes(x)$PRS, time = n))
#   
#   res <- cbind(res,betas)
#   return(res)
# })
# 
# resultQ4 <- do.call(rbind,aux)
# resultQ4$PRS <- factor(resultQ4$PRS,levels=PRSVars)
# resultQ4$PTGENDER <- factor(resultQ4$PTGENDER, levels=c("Female","Male"))
# 
# # Corregim pvalors segons FDR, considerant cada volum independent,
# #  i corregint per 6 comparacions (6 PRS) o models per cada un dels volums
# aux <- lapply(split(resultQ4,resultQ4$Volume),
#               FUN=function(x)p.adjust(x$Significance,method="fdr"))
# resultQ4$Significance_corrected <- do.call(c,aux)
# 
# # Afegim simbols de significació
# resultQ4[paste(names(resultQ4)[5:6],"symbol",sep="_")] <-
#   lapply(resultQ4[names(resultQ4[5:6])], FUN=cut,
#          breaks = c(-Inf,0.001,0.01,0.05,0.1,1),
#          labels=c("***","**","*","·",""))
# 
# 
# plotData <- subset(resultQ4,Volume!="Whole_hippocampus")
# names(plotData)[3] <- "Group"
# heatmapT1 %+%plotData
# subset(resultQ4, Significance_symbol !="",
#        c(Volume,PRS,PTGENDER,Effect,Significance_symbol,Significance_corrected_symbol))


#-----------------------------------------------------------------------------
# d) Q4.2: l'associació entre PRS i volums, estratificant per grup diagnòstic,
#     és diferent segons sexe?
#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# FUNCIO effectPRS_strat_inter 
#   INPUT:
#     - volume: variable volum hippocamp
#     - PRS
#     - group: "CN","MCI","AD indicant grup a estratificar
#
#   TASQUES
#     - A) Subset de les dades segons grup, si escau
#     - B) Ajusta el model lineal amb PRS interacció amb sexe, PRS*PTGENDER,
#          ajustant per edad basal AGE, educació PTEDUCAT i volum intracranial
#          total (ICV)
#     - C) Recupera efecte PRS sobre volum i interval de confiança, segons
#          sexe (PTGENDER) i per les dades estratificades segons grup DX.bl.
#          També retorna significació test wald pel terme interacció.
#
#   OUTPUT
#     - Dataframe amb 9 variables: 
#         * DX.bl: grup diagnòsitc del model estratificat.
#
#         * Volume: subcamp de l'hipocamp amb volum basal com a resposta.
#
#         * PRS: polygenic risk score usat com a variable explicativa.
#
#         * PTGENDER: male o female, segons quin sexe es refereixi l'efecte.
#
#         * Effect: efecte PRS sobre el volum resposta, considerant la
#             interacció amb sexe PTGENDER. Es calcula com:
#             i)  beta_PRS_male = beta_PRS
#             ii) beta_PRS_female = beta_PRS + beta_prs:female
#         * ci_lwr, ci_upr: limits del interval confiança per Effect al 95%
#
#         * Significance: pvalor contrast H0: beta_effect=0. O bé:
#             Male:   H0: beta_effect = beta_PRS = 0
#             Female: H0: beta_effect = beta_PRS + beta_prs:female = 0
#
#         * pvalInteraction: pvalor significacio test wald sobre coeficient
#             interacció beta_prs:female.
#------------------------------------------------------------------------------


effectPRS_strat_inter <- function(volume,PRS,group){
  
  # i) Subset dades segons estrat escollir ·················
  data_strata <- subset(basal,DX.bl == group)
  
  # ii) Construcció formula model ························
  covars <- paste(c("AGE","ICV","PTEDUCAT","PTGENDER"),
                  collapse=" + ")
  
  
  expr_additiu <- paste(volume,
                        paste(c(covars,PRS), collapse=" + "),
                        sep = " ~ ")
  expr_inter <- paste(volume,
                        paste(c(covars,PRS), collapse=" * "),
                        sep = " ~ ")
  
  # iii) Ajust model lineal estratificar per DX.bl: 
  #      main effect PRS + PRS effect segons PTGENDER
  mod_add <- lm(formula = expr_additiu, data = data_strata)        
  mod_int <- lm(formula = expr_inter, data = data_strata)
  
  attributes(mod_int)$PRS <- PRS
  attributes(mod_int)$volume <- volume
  
    
  # iv) Obtenció efecte principal PRS sense interacció amb PTGENDER
  prs_pos <- which(names(coef(mod_add)) == PRS)
  mainEffectPRS <- summary(mod_add)$coefficient[prs_pos,c(1:4)]
  ci <- confint.default(mod_add)[prs_pos,]
  df <- summary(mod_add)$df[2]

  # v) Obtenció efectes PRS segons PTGENDER i significació interacció ····
  emm <- emtrends(mod_int,~PTGENDER,var=PRS)
  emm.df <- as.data.frame(summary(emm,infer=c(T,T)))
  contrast <- as.data.frame(pairs(emm))
  signif <- contrast$p.value <0.05
  
  # v) Dataframe resultats
  res <- data.frame(DX.bl=group,
                    Volume = rep(volume,3),
                    PRS = rep(PRS,3),
                    PTGENDER = c("Global","Female","Male"),
                    Effect = c(mainEffectPRS[1],emm.df[,2]),
                    SE = c(mainEffectPRS[2],emm.df$SE),
                    df = c(df,emm.df$df),
                    t.ratio = c(mainEffectPRS[3],emm.df$t.ratio),
                    ci_lwr = c(ci[1],emm.df$lower.CL),
                    ci_upr = c(ci[2],emm.df$upper.CL),
                    Significance = c(mainEffectPRS[4],emm.df$p.value),
                    #pvalInteraction = c(NA,contrast$p.value,NA),
                    subset = as.logical(c(1,0,0) + signif*c(-1,1,1)))
  
  #rownames(res) <- expr
  return(list("PRSEffect" = res,
              "modSxInteract" = mod_int))
}
#-----------------------------------------------------------------------------

# Provar funció...
effectPRS_strat_inter("Hippocampal_tail","PRS_IEAA",group = "AD")

### i) Construcció matriu combinacions PRS, volums i estrat DX.bl
modGrid_strata <- expand.grid(PRS = PRSVars,
                              Volume = sort(volVars),
                              Dx.bl = levels(basal$DX.bl),
                              stringsAsFactors = F)

### ii) Aplicar funció effectPRS_Strat_inter() a cada combinació
result_strat <- mapply(FUN = effectPRS_strat_inter,
                       volume = modGrid_strata$Volume,
                       PRS = modGrid_strata$PRS,
                       group = modGrid_strata$Dx.bl,
                       SIMPLIFY = F,USE.NAMES = F)

auxQ4 <- lapply(result_strat,"[[","PRSEffect")

### ii.bis) Seleccionar aquells per mostrar
auxQ4 <- lapply(auxQ4,FUN=subset, subset=subset)

### iii) Consolidar llista resultats en un sol data frame
effectPRS_Sx <- do.call(rbind,auxQ4)
effectPRS_Sx$subset <- NULL

# PRS,DX.bl com a factor
effectPRS_Sx$PRS <- factor(effectPRS_Sx$PRS,
                        levels=PRSVars)
effectPRS_Sx$DX.bl <- factor(effectPRS_Sx$DX.bl,
                          levels=c("AD","MCI","CN"))
effectPRS_Sx$PTGENDER <- factor(effectPRS_Sx$PTGENDER,
                             levels = c("Global","Female","Male"))

effectPRS_Sx$alpha <- with(effectPRS_Sx,ifelse(Significance<0.1,1,.9))

# iv) Afegir simbols de significació
# effectPRS[paste(c("Significance","Interact"),"symbol",sep="_")] <-
#   lapply(effectPRS[c("Significance","pvalInteraction")],
#          FUN=cut,breaks = c(-Inf,0.001,0.01,0.05,0.1,1),
#          labels=c("***","**","*","^",""))
# effectPRS[is.na(effectPRS$Interact_symbol),"Interact_symbol"] <- ""

# RESULTAT
head(effectPRS_Sx)


# v) PANEL DE RESULTATS

# Auxiliar
  # separacio pointrange female/male
  pos_dodge <- position_dodge(width=0.8) 

  # Etiquetes panell
  volume.labs <- volVars
  names(volume.labs) <- volVars
  volume.labs[c(7:11,13)] <- c("Hippocampal\nfissure",
                           "Hippocampal\ntail",
                           "Molecular\nlayer hp",
                           "Para-\nsubiculum",
                           "Pre-\nsubiculum",
                           "Whole\nhippocampus")
  PRS.labs <- sapply(strsplit(PRSVars,split="_"),"[",2)
  names(PRS.labs) <- PRSVars
  PRS.labs[PRS.labs=="ADnoAPOE"] <- "AD\nno APOE"
  PRS.labs[PRS.labs=="LONGEVITYnoAPOE"] <- "LONGEVITY\nno APOE"
  
# PANEL
  
  
ggplot(effectPRS_Sx,aes(x=DX.bl)) +
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  geom_pointrange(aes(y=Effect , ymin=ci_lwr , ymax=ci_upr, 
                      color=PTGENDER, alpha=alpha),
                  position = pos_dodge,show.legend = F) +
#  geom_text(aes(y=Effect,label=Significance_symbol, color =PTGENDER),size=5,
#           position = position_dodge(width=0.1), show.legend = F)+
  facet_grid(rows=vars(PRS),cols=vars(Volume),scales = "free",
             labeller=labeller(Volume=volume.labs,PRS=PRS.labs))+
  scale_color_manual(values = c("#000000","#D95F02","#7570B3")) + 
  scale_alpha(range = c(0.3,1))+
  coord_flip()+
  
  theme_bw() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90,hjust=1)) +
  ggtitle("Association results for direct effects of polygenic risk scores on hippocampal volumes<br/>for female and male and stratified by diagnostic group.") +
  labs(subtitle = "Adjusted for Age at baseline, sex, intracranial volume and years of education. <br/>
       <span style = 'color: #000000;'>**Global**</span> , <span style = 'color: #D95F02;'>**Female**</span> and <span style = 'color: #7570B3;'>**Male**</span>" ,
       caption="**Source**: ADNI project | **Units** Linear coefficient and 95% confidence intervals of PRS effect, for female and male when sex interaction with PRS is significant; or main effect otherwise in models T.3.<br/>Non-significant coeficients are shown with transparency, according to Wald-test p-value at level 0.1.
       After correcting for multiple comparisons, none of these effects remained significant.") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_markdown(size = 16, hjust = 0))




### EXPORTAR PRS EFFECT models T1,T2,T3
auxT3 <- lapply(result_strat,"[[","PRSEffect")
auxT3 <- do.call(rbind,auxT3)
head(auxT3)

prsEffectList <- list("Models T1" = subset(modGrid,DX.bl=="All",c(2,1,4:7,10,11)),
                      "Models T2" = resultQ2[1:9],
                      "Models T3" = auxT3[c(2,3,1,4,5:12)]
                      )

prsEffectList <- lapply(prsEffectList,function(x){
  res <- x
  res$model <- with(res,paste(Volume,PRS,sep=":"))
  res <- res[c("model",names(x))]
  rownames(res) <- NULL
  return(res)
})


# Build Excel worksheet
wb <- createWorkbook()
invisible(lapply(names(prsEffectList),function(x)addWorksheet(wb,x)))
invisible(
  lapply(1:length(prsEffectList),
         function(i) writeDataTable(wb = wb,
                                    x = prsEffectList[[i]],
                                    sheet = names(prsEffectList)[i],
                                    withFilter = T,
                                    tableStyle = "TableStyleLight1")
  )
)

path <- here("results",paramT,"2_Modelling","1_baselineVolumes_groupAdjustedDifferences",
             "taules","crossSecModels_prsEffect.xlsx")

saveWorkbook(wb,file=path,overwrite = T)


### EXPORT parameter estimates models T1, T2, T3

# ModelsList
modelsT1 <- lapply(auxQ1[modGrid$DX.bl=="All"],"[[","ModelGlobal")
modelsT2 <- lapply(auxQ2,"[[","modDXInteract")
modelsT3 <- lapply(result_strat,"[[","modSxInteract")

# Parameters
modelParameters <- lapply(list(modelsT1,modelsT2,modelsT3),
                          FUN = enambleTable_xsecModels)
names(modelParameters) <- paste0("Models T",1:3)

# Build Excel worksheet
wb <- createWorkbook()
invisible(lapply(names(modelParameters),function(x)addWorksheet(wb,x)))
invisible(
  lapply(1:length(modelParameters),
         function(i) writeDataTable(wb = wb,
                                    x = modelParameters[[i]],
                                    sheet = names(modelParameters)[i],
                                    withFilter = T,
                                    tableStyle = "TableStyleLight1")
  )
)

path <- here("results",paramT,"2_Modelling","1_baselineVolumes_groupAdjustedDifferences",
             "taules","crossSecModels_parameters.xlsx")

saveWorkbook(wb,file=path,overwrite = T)
