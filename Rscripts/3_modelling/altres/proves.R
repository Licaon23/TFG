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
load(here(cleandataDir,"dataForAnalysis.RData")) # load data
library(lmerTest)
library(mixedup)
library(sjPlot)
library(ggeffects)

dd$ICV_cm3 <- dd$ICV/1000
dd$time_ageCenter <- scale(dd$time_age,scale=F)


volVars <- sort(names(basal)[14:26])
PRSVars <- names(basal)[7:12]
RiskVars <- names(basal)[27:32]
adjustVars <- c("ICV_cm3","DX.bl","PTGENDER","PTEDUCAT")

#set.seed(533)
ptid <- tapply(dd[,"PTID"], dd$DX.bl,sample,size=10,replace=F)
ptid <- Reduce("c",ptid)
ggplot(subset(dd,PTID%in%ptid),
       aes(x=time_ageCenter, y=CA1,group=PTID))+
  geom_point(size=1)+
  geom_line()+
  facet_wrap(~DX.bl)


# Model de mitjanes per individu
m1 <- lmer(CA1 ~ (1|PTID), data=dd)



prova <- subset(dd,PTID %in% ptid,
                select = c(PTID,time_ageCenter,DX.bl,CA1))

prova$pred <- predict(m1,newdata=prova)

p <- ggplot(prova,
       aes(x=time_ageCenter, y=CA1,group=PTID))+
  geom_point(size=1)+
  geom_line()+
  geom_line(aes(y=pred), col ="red")+
  facet_wrap(~DX.bl)
p



# Afegim random slope

m2 <- lmer(CA1 ~ (1+time_ageCenter|PTID),data=dd)
summary(m2)

prova2 <- subset(dd,PTID %in% ptid)

prova2$pred <- predict(m2,newdata=prova2)

p %+% prova2
ranova(m2)

# Terme quadràtic?

m21 <- lmer(CA1 ~ (1+time_ageCenter + I(time_ageCenter^2)|PTID),data=dd)
summary(m21)
# Singular fit. La variància del terme quadràtic s'ha estimat gairebé nula.


# Ens quedem amb intercept i slope random
# Afegim efecte fix del temps

m3 <- update(m2, ~time_ageCenter + .)
summary(m3)
anova(m3)

prova3 <- subset(dd,PTID %in% ptid)

prova3$pred <- predict(m3,newdata=prova3)
#prova$fix <- predict(m3,newdata=prova,re.form=~0)
p%+%prova3

#Afegim terme quadràtic
m31 <- lmer(CA1 ~ time_ageCenter + I(time_ageCenter^2) + 
              (1+time_ageCenter|PTID),data=dd)
anova(m31)

# Afegim PRS i interacció
m4 <- lmer(CA1 ~ time_ageCenter + I(time_ageCenter^2) + PRS_AD + 
             PRS_AD:time_ageCenter + 
             (1+time_ageCenter|PTID),data=dd)
anova(m4)

# Afegin altres covariables
m5 <- lmer(CA1 ~ time_ageCenter + I(time_ageCenter^2) + PRS_AD + 
             PRS_AD:time_ageCenter + ICV_cm3 +
             (1+time_ageCenter|PTID),data=dd)
anova(m5)

m51 <- lmer(CA1 ~ time_ageCenter + I(time_ageCenter^2) + PRS_AD + 
             PRS_AD:time_ageCenter + ICV_cm3 + DX.bl +
             (1+time_ageCenter|PTID),data=dd)
anova(m51)

m52 <- lmer(CA1 ~ time_ageCenter + I(time_ageCenter^2) + PRS_AD + 
              PRS_AD:time_ageCenter + ICV_cm3 + DX.bl + PTGENDER +
              (1+time_ageCenter|PTID),data=dd)
anova(m52)

m53 <- lmer(CA1 ~ time_ageCenter + I(time_ageCenter^2) + PRS_AD*time_ageCenter + 
              ICV_cm3 + DX.bl + PTGENDER + PTEDUCAT +
              (1+time_ageCenter|PTID),data=dd)
anova(m53)
summary(m53)



# FINAL EFECTES FIXOS
volume <- volVars[1]
PRS <- PRSVars[1]

makeFormula <- function(volume,PRS,adjustVars=NULL, random=c(T,T)){
  
  # Random terms
  randomTerms <- c("(1|PTID)","(1+time_ageCenter)")
  randomTerms <- paste(randomTerms[random],collapse = " + ")
  
  # Formula building
  time <- paste0(c("time_ageCenter","I(time_ageCenter^2)"),collapse=" + ")
  prs <- paste0(c(time,paste(PRS,"time_ageCenter",sep="*")),collapse=" + ")
  covars <- paste0(c(prs,adjustVars),collapse=" + ")
  random <- paste0(c(covars,randomTerms),collapse=" + ")
  expr <- paste0(c(volume,random), collapse = " ~ ")

  return(expr)  
}

# Trajectòria parabòlica + no efecte PRS sobre decline
exp <- makeFormula(volume = "parasubiculum",
                   PRS = "PRS_ADnoAPOE",
                   adjustVars=adjustVars,
                   random = c(T,F))


mDef <- lme4::lmer(formula = exp, data=dd)
extract_fixed_effects(mDef)
g <- ggpredict(mDef,terms=c("time_ageCenter[all]","PRS_ADnoAPOE[-1,0,1]"))
plot(g)

# Trajectòria parabòlica + efecte PRS sobre decline
exp <- makeFormula(volume = "Hippocampal_tail",
                   PRS = "PRS_AD",
                   adjustVars=adjustVars,
                   random = c(T,F))


mDef <- lme4::lmer(formula = exp, data=dd)
extract_fixed_effects(mDef)


# Trajectòria no parabòlica + no efectePRS sobre decline
exp <- makeFormula(volume = "HATA",
                   PRS = "PRS_EEAA",
                   adjustVars=adjustVars,
                   random = c(T,F))


mDef <- lme4::lmer(formula = exp, data=dd)
extract_fixed_effects(mDef)
g <- ggpredict(mDef,terms=c("time_ageCenter[all]","PRS_EEAA[-1,0,1]"))
plot(g)



# TOTS ELS MODELS ---------------------------

modGrid <- expand.grid(PRS = PRSVars,
                       Volume = sort(volVars),
                       #DX.bl = c("All",levels(basal$DX.bl)),
                       stringsAsFactors = F)

formulas <- mapply(FUN=makeFormula,
                   volume=modGrid$Volume,
                   PRS=modGrid$PRS,
                   MoreArgs = list(adjustVars = adjustVars,
                                   random = c(T,F)), 
                   SIMPLIFY = F)

models <- lapply(formulas,lme4::lmer,data=dd)
lapply(models,function(x)extract_fixed_effects(x)[c(2:4,10),])

res <- lapply(models,function(x)extract_fixed_effects(x)[10,])
res <- do.call(rbind,res)

res <- cbind(modGrid,res[-1])
signif <- subset(res,p_value<0.05)
signif[order(signif$PRS),]

#25 a 30 # GC.ML.DG
#49 a 54 #molecular layer
#55 a 60














































ranova(mDef)
# No treiem el random intercept

library(ggeffects)


mdif <- ggpredict(mDef,terms=c("time_ageCenter[all]","PRS_AD[-1,0,1]"))

mdef2 <- update(mDef,~.+PRS_AD:DX.bl + time_ageCenter:DX.bl +
                  time_ageCenter:PRS_AD:DX.bl)
mdif2 <- ggpredict(mdef2,terms=c("time_ageCenter[all]","PRS_AD[-1,0,1]",
                                 "DX.bl"))

plot(mdif)
plot(mdif2,add.data = T,dot.size = .5)

mdif
hypothesis_test(m)
resid_panel(mDef)
lapply(ranef(mDef)$PTID,FUN=hist)
# hi ha un individu amb un efecte aleatori molt allunyat de la mitjana!
# estudiar-lo


### L'efecte de PRS depen de l'estat diagnostic?
mDef2 <- update(mDef, ~.+time_ageCenter:DX.bl + PRS_LONGEVITY:DX.bl + 
                  time_ageCenter:PRS_LONGEVITY:DX.bl)
anova(mDef2)

mdif2 <- ggpredict(mDef2, terms=c("time_ageCenter[all]", "PRS_LONGEVITY[-1,0,1]",
                                  "DX.bl"))
plot(mdif2)


plot_model(mDef2,
           type="pred",
           terms=c("time_ageCenter[all]","PRS_AD[-1,0,1]","DX.bl"))


### L'efecte del PRS depen del sexe?
mDef3 <- update(mDef, ~.+time_age2:PTGENDER + PRS_AD:PTGENDER)

anova(mDef3)
plot_model(mDef3,
           type="pred",
           terms=c("time_age2[all]","PRS_AD[-2,0,2]","PTGENDER"))



# 
# # Model mixte pel evolució volum hippocampal
# expr1 <- paste0(c("time_age","I(time_age^2)",adjustVars,"(1|PTID)"),collapse=" + ")
# expr1 <- paste0(c(volume,expr1), collapse = " ~ ")
# 
# expr2 <- paste0(c("time_age",adjustVars,"(1|PTID)"),collapse=" + ")
# expr2 <- paste0(c(volume,expr2), collapse = " ~ ")
# 
# expr3 <- paste0(c("bs(time_age, df=3)",adjustVars,"(1|PTID)"),collapse=" + ")
# expr3 <- paste0(c(volume,expr3), collapse = " ~ ")
# 
# lmer_out_timequadratic <- lme4::lmer(expr1, data=dd)
# lmer_out_timeLinear <- lme4::lmer(expr2, data=dd)
# lmer_out_timeSpline <- lme4::lmer(expr3, data=dd)
# 
# # Variancia
# extract_VarCorr(lmer_out_timequadratic)
# extract_VarCorr(lmer_out_timeLinear)
# 
# # Efectes fixes
# extract_fixed_effects(lmer_out_timequadratic)
# extract_fixed_effects(lmer_out_timeLinear)
# 
# #Plots
# p1 <- sjPlot::plot_model(lmer_out_timequadratic,type="pred",terms="time_age[all]")
# p2 <- sjPlot::plot_model(lmer_out_timeLinear,type="pred",terms="time_age") 
# p1+scale_y_continuous(limits = c(600,1400))+ 
# p2+scale_y_continuous(limits = c(600,1400))
# 
# # Model comparisons
# anova(lmer_out_timeLinear,lmer_out_timequadratic) 
# anova(lmer_out_timeLinear,lmer_out_timeSpline) 
# 
# AIC(lmer_out_timeLinear)
# AIC(lmer_out_timequadratic)
# AIC(lmer_out_timeSpline)

# El terme quadràtic és estadísticament significatiu.
# I segons test ANOVA el model amb el terme quadràtic explica significativament
# més la resposta

# Tanmateix, els AIC no son gaire diferents. De fet, amb el terme quadràtic
# l'AIC és lleugerament major.

# Falla convergir amb random slope per PTID
# Alerta de escala de variables al afegir time_age al

# Hi ha algun individu extrem (analitzar)


#########################
# Efecte PRS
volume <- volVars[8]
PRS <- PRSVars[6]

expr <- paste0(c(adjustVars,"time_age"),collapse=" + ")
expr <- paste0(c(expr,PRS), collapse = " * ")
expr <- paste0(c(expr,"(1|PTID)"),collapse=" + ")
expr <- paste0(c(volume,expr), collapse = " ~ ")


lmer_out <- lme4::lmer(expr,data=dd)
extract_fixed_effects(lmer_out)
extract_VarCorr(lmer_out)
#summary(lmerTest::lmer(expr,data=dd))

# sjPlot::plot_model(lmer_out, type="pred", 
#                    terms= c("time_age[60,90]","PTGENDER","DX.bl")) +
#   theme_sjplot()

plot_model(lmer_out,type="int")+theme_bw()


