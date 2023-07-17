###############################################################################
### EXPLORATORY DATA ANALYSIS                                                 #
###   Baseline values description                                             #
###   Version 1                                                               #
###   January 2023                                                            #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###                                                                           #
###############################################################################
source("Rscripts/libraries.R")
source(here(RscriptsDir,"2_exploratory","descSummary_bis.R")) # aux function

load(here(cleandataDir,"dataForAnalysis_3T.RData"))
paramT <- "Dataset_3T"

# load(here(cleandataDir,"dataForAnalysis_1T.RData"))
# paramT <- "Dataset_1T"

#============================================================================
# B) DESCRIPTIVA CARACTERISTIQUES BASALS DELS PACIENTS
#============================================================================

# Variables modificades
basal$ICV_cm3 <- basal$ICV/1000
codebook <- rbind(codebook,c("ICV_cm3","Intracranial volume (cm3)",NA))

basal$APOE4 <- with(basal,factor(APOE4,levels = 0:2, 
                                 labels = paste(0:2,"allele")))

# Posicio Variables
aux <- colnames(basal)
riskVars <- aux[grep("predisposition",aux)]
prsVars <- aux[grep("PRS",aux)]
volVars <- sort(aux[17:29])

#-----------------------------------------------------------------------------
# i) Taula estadístics descriptius segons grup diagnòstic
#-----------------------------------------------------------------------------
vars <- data.frame(Variable = c("AGE","PTGENDER","PTEDUCAT","APOE4","ICV_cm3",
                                riskVars,volVars))
vars <- merge(vars,codebook[,1:2],by="Variable",sort = F)
vars$nonNormal <- c(F,NA,T,NA,F,rep(NA,7),rep(F,13))

aux <-lapply(c(list(Total=basal) , split(basal,basal$DX.bl)), 
             descSummary,
             vars = vars$Variable, printNames=vars$Label,
             nonNormal=vars$nonNormal)


aux <- lapply(names(aux),function(x){
  d <- aux[[x]]
  names(d)[names(d) == "summary"] <- paste("summary",x,sep=".")
  return(x=d)})

taula_dx <- Reduce(function(x,y){merge(x,y,by="variable",sort=F)},aux)


pvals <- mapply(groupContrast,x=vars$Variable,nonNormal=vars$nonNormal,
                MoreArgs = list(by="DX.bl",data=basal),SIMPLIFY = F)

pvals <- do.call(c,pvals)
pvals <- c(pvals[1:2],rep("",2),pvals[3:4],rep("",3),pvals[-c(1:4)])
taula_dx$pvalue <- pvals
colnames(taula_dx) <- c("","Total",levels(basal$DX.bl),"Pvalue")

taula_dx <- rbind(c("Demographic characteristics",rep("",5)),
                  taula_dx[1:10,],
                  c("Genetic Scores",rep("",5)),
                  taula_dx[11:17,],
                  c("Hippocampal volumes",rep("",5)),
                  taula_dx[18:30,])

#Prevalence groups
counts <- table(basal$DX.bl)
n <- sprintf("N = %.i",c(sum(counts),counts))
p <- c("",paste0("(",sprintf("%.0f",prop.table(counts)*100),"\\%)"))
prev_dx <- format(c(n,p),justify = "centre")



tab <- kbl(taula_dx,align = "lccccc",
           col.names = c("",prev_dx[1:4],""),
           caption = "Baseline characteristics of subjects by diagnostic group",
           row.names = F) 
tab <- add_header_above(tab,c("","Total",levels(basal$DX.bl),"Pvalue"),
                        line = F,bold=T)
                        
tab <- kable_classic(tab,full_width=F,font="Cambria",)
formattedTab <- footnote(tab,
  general = c("Results are reported as mean(standrad deviation) or median [Q1,Q3] quartile for normal and non-normal numeric variables, respectively;
  and as number of cases and proportion n(p%) for categorical variables.",
  "Differences between diagnostic groups where tested using ANOVA or Kruskal-Wallis test for numeric Normal or non-Normal variables, respectively;
   and using a Chi-squared test for categorical variables. P-values are reported.",
  "Diagnostic groups: Controls (CN), Alzheimer's Disease (AD), Mild Cognitive Impairment (MCI)"),
   title_format = c("bold")) 

formattedTab <- row_spec(formattedTab,c(0,1,12,20),bold=T)
formattedTab_dxbl <- add_indent(formattedTab, position =c(2:3,6:7,11,13:19,21:33),
                                level_of_indent = .5)
formattedTab_dxbl <- add_indent(formattedTab_dxbl, position =c(4:5,8:10),
                                level_of_indent = 1)

# Guardar i visualitzar taula de resultats
file <- "result_baselineCharactByDxbl.png"
path <- here("results",paramT,"1_DescriptiveAnalysis",
             "1_baselineCharacteristics")

if(!(any(grepl(file,list.files(path))))){
  save_kable(formattedTab_dxbl,file=here(path,file),zoom=5)
}
formattedTab_dxbl

#-----------------------------------------------------------------------------
# ii) Taula estadístics descriptius segons sexe biològic
#-----------------------------------------------------------------------------
vars <- data.frame(Variable = c("AGE","DX.bl","PTEDUCAT","APOE4","ICV_cm3",
                                riskVars,volVars))
vars <- merge(vars,codebook[,1:2],by="Variable",sort = F)
vars$nonNormal <- c(F,NA,T,NA,F,rep(NA,7),rep(F,13))

aux <-lapply(c(list(Total=basal) , split(basal,basal$PTGENDER)), 
             descSummary,
             vars = vars$Variable, printNames=vars$Label,
             nonNormal=vars$nonNormal)

aux <- lapply(names(aux),function(x){
  d <- aux[[x]]
  names(d)[names(d) == "summary"] <- paste("summary",x,sep=".")
  return(x=d)})

taula_sex <- Reduce(function(x,y){merge(x,y,by="variable",sort=F)},aux)


pvals <- mapply(groupContrast,x=vars$Variable,nonNormal=vars$nonNormal,
                MoreArgs = list(by="PTGENDER",data=basal),SIMPLIFY = F)

pvals <- do.call(c,pvals)
pvals <- c(pvals[1:2],rep("",3),pvals[3:4],rep("",3),pvals[-c(1:4)])
taula_sex$pvalue <- pvals
colnames(taula_sex) <- c("","Total",levels(basal$PTGENDER),"Pvalue")

taula_sex <- rbind(c("Demographic characteristics",rep("",4)),
                   taula_sex[1:11,],
                  c("Genetic Scores",rep("",4)),
                  taula_sex[12:18,],
                  c("Hippocampal Volumes",rep("",4)),
                  taula_sex[19:31,])


#Prevalence groups
count <- table(basal$PTGENDER)
n <- sprintf("N = %.i",c(sum(count),count))
p <- c("",paste0("(",sprintf("%.0f",prop.table(count)*100),"\\%)"))
prev_sex <- format(c(n,p),justify = "centre")


tab <- kbl(taula_sex,align = "lcccc",
           col.names = c("",prev_sex[1:3],""),
           caption = "Baseline characteristics of subjects by sex",
           row.names = F) 
tab <- add_header_above(tab,c("","Total",levels(basal$PTGENDER),"Pvalue"),
                        line = F,bold=T)

tab <- kable_classic(tab,full_width=F,font="Cambria")
formattedTab <- footnote(tab,
                         general = c("Results are reported as mean(standrad deviation) or median [Q1,Q3] quartile for normal and non-normal numeric variables, respectively;
  and as number of cases and proportion n(p%) for categorical variables.",
                                     "Differences between sex groups where tested using ANOVA or Kruskal-Wallis test for numeric Normal or non-Normal variables, respectively;
   and using a Chi-squared test for categorical variables. P-values are reported.",
                                     "Diagnostic groups: Controls (CN), Alzheimer's Disease (AD), Mild Cognitive Impairment (MCI)"),
                         title_format = c("bold")) 

formattedTab <- row_spec(formattedTab,c(0,1,13,21),bold=T)
formattedTab_sex <- add_indent(formattedTab, position =c(2:3,7:8,12,14:20,22:34),
                               level_of_indent = .5)
formattedTab_sex <- add_indent(formattedTab_sex, position =c(4:6,9:11),
                               level_of_indent = 1)



file <- "result_baselineCharactBySex.png"
path <- here("results",paramT,"1_DescriptiveAnalysis",
             "1_baselineCharacteristics")

if(!(any(grepl(file,list.files(path))))){
  save_kable(formattedTab_sex,file=here(path,file),zoom=5)
}
formattedTab_sex


#-----------------------------------------------------------------------------
# v) Boxplots PRS segons grup diagnòstic
#-----------------------------------------------------------------------------
prs <- basal[,c("PTID","DX.bl",prsVars)]
prs_long <- reshape(prs, direction="long",
                    idvar="PTID",
                    varying = list(names(prs)[3:9]),
                    v.names="value",
                    timevar="PRS",
                    times=names(prs)[3:9])
rownames(prs_long) <- NULL

prs_long$PRS <- factor(prs_long$PRS,
                       levels = paste("PRS",c("AD","ADnoAPOE","FTD",
                                              "EEAA","IEAA","LONGEVITY",
                                              "LONGEVITYnoAPOE"),sep="_"))

stat.test <- prs_long %>%
  group_by(PRS) %>%
  t_test(value ~ DX.bl) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x="DX.bl",scales = "free_y",step.increase = 0.22)

stat.test$p.adj.signif <- cut(stat.test$p.adj,
                              breaks=c(-Inf,0.001,0.01,0.05,1),
                              labels = c("***","**","*","ns"))


ggplot(data=prs_long, aes(x=DX.bl,y=value)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1,aes(color=DX.bl),show.legend = F, size=0.1,
              alpha=0.5) +
  stat_pvalue_manual(stat.test,label="p.adj.signif",hide.ns = T) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1))) +
  scale_color_manual(values = c("#1B9E77","#D95F02","#7570B3")) +
  facet_wrap(~PRS,scales = "free") + 
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank()) +
  
  ggtitle("Polygenic Risk scores distribution by diagnostic group") +
  labs(subtitle = "<span style = 'color: #1B9E77;'>**Controls**</span>, <span style = 'color: #D95F02;'>**Mild cognitive impairment**</span> and <span style = 'color: #7570B3;'>**Alzheimer`s disease**</span>" ,
       caption="**Source**: ADNI project | **Diagnostic groups** CN: controls, AD: Alzheimer`s disease, MRI: mild cognitive impairment | **Units** PRS values have been standarized <br/>
       Pairwise t-test have benn performed, and p-values have been adjusted for multiple comparison using FDR method. P-values: \\***<.001 ; \\**<.01; \\*<.05") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = -.05))

ggsave(filename = "PRSDistBaseByDxbl.png",
       path = here("results",paramT,"1_DescriptiveAnalysis",
                   "1_baselineCharacteristics"),
       width = 25,height = 15,units = "cm")


#-----------------------------------------------------------------------------
# v) Boxplots PRS segons Sexe
#-----------------------------------------------------------------------------
prs <- basal[,c("PTID","PTGENDER",prsVars)]
prs_long <- reshape(prs, direction="long",
                    idvar="PTID",
                    varying = list(names(prs)[3:9]),
                    v.names="value",
                    timevar="PRS",
                    times=names(prs)[3:9])
rownames(prs_long) <- NULL
prs_long$PRS <- factor(prs_long$PRS,
                       levels = paste("PRS",c("AD","ADnoAPOE","FTD",
                                              "EEAA","IEAA","LONGEVITY",
                                              "LONGEVITYnoAPOE"),sep="_"))


stat.test <- prs_long %>%
  group_by(PRS) %>%
  t_test(value ~ PTGENDER) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x="PTGENDER",scales = "free_y",step.increase = 0.22)

stat.test$p.adj.signif <- cut(stat.test$p.adj,
                              breaks=c(-Inf,0.001,0.01,0.05,1),
                              labels = c("***","**","*","ns"))


ggplot(data=prs_long, aes(x=PTGENDER,y=value)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1,aes(color=PTGENDER),show.legend = F, size=0.1,
              alpha=0.5) +
  stat_pvalue_manual(stat.test,label="p.adj.signif",hide.ns = T) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1))) +
  scale_color_manual(values = c("#D95F02","#7570B3")) +
  facet_wrap(~PRS,scales = "free") + 
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank()) +
  
  ggtitle("Polygenic Risk scores distribution by sex") +
  labs(subtitle = "<span style = 'color: #D95F02;'>**Female**</span> and <span style = 'color: #7570B3;'>**Male**</span>" ,
       caption="**Source**: ADNI project | **Diagnostic groups** CN: controls, AD: Alzheimer`s disease, MRI: mild cognitive impairment | **Units** PRS values have been standarized <br/>
       Pairwise wilcoxon test have been performed, and p-values have been adjusted for multiple comparison using FDR method. P-values: \\***<.001 ; \\**<.01; \\*<.05") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = -.05))

ggsave(filename = "PRSDistBaseBySex.png",
       path = here("results",paramT,"1_DescriptiveAnalysis",
                   "1_baselineCharacteristics"),
       width = 25,height = 15,units = "cm")





#-----------------------------------------------------------------------------
# iii) Boxplots volums hipocampals segons grup diagnòstic
#-----------------------------------------------------------------------------

volumes_basal <- basal[,c("PTID","DX.bl","ICV_cm3",volVars)]
volumes_long <- reshape(volumes_basal, direction="long",
                       idvar="PTID",
                       varying = list(names(volumes_basal)[3:16]),
                       v.names="value",
                       timevar="subfield",
                       times=names(volumes_basal)[3:16])
rownames(volumes_long) <- NULL

stat.test <- volumes_long %>%
  group_by(subfield) %>%
  t_test(value ~ DX.bl) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x="DX.bl",scales = "free_y",step.increase = 0.22)

stat.test$p.adj.signif <- cut(stat.test$p.adj,
                              breaks=c(-Inf,0.001,0.01,0.05,1),
                              labels = c("***","**","*","ns"))

volLabels <- unique(volumes_long$subfield)
pos <- grep("Whole|ICV",volLabels)
volumes_long$subfield <- factor(volumes_long$subfield,
                                levels = c(volLabels[-pos],
                                           rev(volLabels[pos])))

ggplot(data=volumes_long, aes(x=DX.bl,y=value)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1,aes(color=DX.bl),show.legend = F, size=0.05,
              alpha=0.5) +
  stat_pvalue_manual(stat.test,label="p.adj.signif",hide.ns = T) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1))) +
  scale_color_manual(values = c("#1B9E77","#D95F02","#7570B3")) +
  facet_wrap(~as.factor(subfield),scales = "free") + 
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank()) +
  
  ggtitle("Hippocampal subfield volume distributions at baseline by diagnostic group") +
  labs(subtitle = "<span style = 'color: #1B9E77;'>**Controls**</span>, <span style = 'color: #D95F02;'>**Mild cognitive impairment**</span> and <span style = 'color: #7570B3;'>**Alzheimer`s disease**</span>" ,
       caption="**Source**: ADNI project | **Diagnostic groups** CN: controls, AD: Alzheimer`s disease, MRI: mild cognitive impairment | **Units** Volume measurements in mm3 <br/>
       Multiple comparisons t-test have been performed, adjusting their p-values with FDR method. P-values: \\***<.001 ; \\**<.01; \\*<.05") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = -.05))

  

ggsave(filename = "volumeDistBaseByDxbl.png",
       path = here("results",paramT,"1_DescriptiveAnalysis",
                   "1_baselineCharacteristics"),
       width = 25,height = 15,units = "cm")

#-----------------------------------------------------------------------------
# iv) Boxplots volums hipocampals segons sexe biològic
#-----------------------------------------------------------------------------
volumes_basal <- basal[,c("PTID","PTGENDER","ICV_cm3",volVars)]
volumes_long <- reshape(volumes_basal, direction="long",
                        idvar="PTID",
                        varying = list(names(volumes_basal)[3:16]),
                        v.names="value",
                        timevar="subfield",
                        times=names(volumes_basal)[3:16])
rownames(volumes_long) <- NULL

stat.test <- volumes_long %>%
  group_by(subfield) %>%
  t_test(value ~ PTGENDER) %>% 
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position(x="PTGENDER",scales = "free_y",step.increase = 0.22)

stat.test$p.adj.signif <- cut(stat.test$p.adj,
                              breaks=c(-Inf,0.001,0.01,0.05,1),
                              labels = c("***","**","*","ns"))

volLabels <- unique(volumes_long$subfield)
pos <- grep("Whole|ICV",volLabels)
volumes_long$subfield <- factor(volumes_long$subfield,
                                levels = c(volLabels[-pos],
                                           rev(volLabels[pos])))



ggplot(data=volumes_long, aes(x=PTGENDER,y=value)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1,aes(color=PTGENDER),show.legend = F, size=0.05,
              alpha=0.5) +
  stat_pvalue_manual(stat.test,label="p.adj.signif",hide.ns = T) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1))) +
  scale_color_manual(values = c("#D95F02","#7570B3")) +
  facet_wrap(~as.factor(subfield),scales = "free") + 
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank()) +
  
  ggtitle("Hippocampal subfield volume distributions at baseline by sex") +
  labs(subtitle = "<span style = 'color: #D95F02;'>**Female**</span> and <span style = 'color: #7570B3;'>**Male**</span>" ,
       caption="**Source**: ADNI project |  | **Units** Volume measurements in mm3 | T-tests were performed to check for significance in the observed differences.P-values: \\***<.001 ; \\**<.01; \\*<.05") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = -.05))


ggsave(filename = "volumeDistBaseBySex.png",
       path = here("results",paramT,"1_DescriptiveAnalysis",
                   "1_baselineCharacteristics"),
       width = 25,height = 15,units = "cm")



#-----------------------------------------------------------------------------
# vii) Histogrames volums
#-----------------------------------------------------------------------------
volumes_basal <- basal[,c("PTID",volVars,"ICV_cm3")]
volumes_long <- reshape(volumes_basal, direction="long",
                        idvar="PTID",
                        varying = list(names(volumes_basal)[2:15]),
                        v.names="value",
                        timevar="subfield",
                        times=names(volumes_basal)[2:15])
rownames(volumes_long) <- NULL

ggplot(data=volumes_long,aes(x=value)) +
  geom_histogram(aes(y = ..density..), fill="#7570B3",alpha=0.6)+
  geom_density(color="#D95F02",linewidth=.5)+
  facet_wrap(~subfield,scales = "free") +
  
  ggstatsplot::theme_ggstatsplot() +
  theme(strip.background = element_blank(),
        plot.margin = margin(1,1,.5,1,"cm"),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggtitle("Hippocampal subfield volume distributions at baseline") +
  labs(caption="**Source**: ADNI project |  | **Units** Volume measurements in mm3 ") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = -.05))

ggsave(filename = "volumeDistBaseHist.png",
       path = here("results",paramT,"1_DescriptiveAnalysis",
                   "1_baselineCharacteristics"),
       width = 25,height = 15,units = "cm")


#-----------------------------------------------------------------------------
# vii) Correlació volums hipocampals i edad basal
#-----------------------------------------------------------------------------
aux <- basal[,c("AGE",volVars,"ICV_cm3")]
c1 <- cor(aux)
pval <- cor.mtest(aux,adjust = "none")$p

png(filename = here("results",paramT,"1_DescriptiveAnalysis",
                    "1_baselineCharacteristics","pearsonCorr.png"),
    width = 696,height = 527,units = "px")

corrplot(c1,type="lower",cl.ratio=0.2,tl.srt = 0.45,
         p.mat=pval,insig = "label_sig",sig.level = 0.05,pch="*",pch.cex = .9,
         title = "Pearson correlation coefficients",mar=c(1,2,2,2))


dev.off()

# Correlació negativa signiciativa entre edat i tots els volums hipocampals.
# A més, també correlació positiva de ICV i tots els volums.

# Per altra banda, ICV en homes era major que en dones. Cal ajustar per ICV
# si es vol eliminar-ne l'efecte abans d'avaluar l'efecte del sexe o Risc PRS


# PANEL SCATTER PLOT PER PRS NUMERIQUES ------------------------------------
my_fn <- function(data,mapping,...){
  p <- ggplot(data=data,mapping=mapping) +
    geom_point(color="orange",size=0.1,alpha=0.6) +
    geom_smooth(method="gam",color="lightblue",linewidth=1,se=F,...) 
  p
}

ggpairs(basal[prsVars],
        columnLabels = paste("PRS",c("AD","AD without APOE","FTD",
                                     "IEAA","EEAA","Longevity",
                                     "LONGEVITYnoAPOE")),
        lower=list(continuous=my_fn)) + 
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(face="bold",size=12)) +
  ggtitle("Panel plot for numeric variables") +
  labs(subtitle = "Polygenic Risk Scores",
       caption="**Source**: ADNI project |") +
  theme(plot.caption = element_markdown(hjust = 0, lineheight = 1.5),
        plot.subtitle = element_markdown(size = 12, hjust = 0),
        plot.title = element_text(size = 16, hjust = 0))


ggsave(filename = "panelPRSscatterplot.png",
       path = here("results",paramT,"1_DescriptiveAnalysis",
                   "1_baselineCharacteristics"),
       width = 25,height = 25,units = "cm")





