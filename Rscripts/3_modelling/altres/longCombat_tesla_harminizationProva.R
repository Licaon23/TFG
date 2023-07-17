source("Rscripts/libraries.R")
library(FactoMineR)


mri_meta3 <- read.csv(here(rawdataDir,"MRI_metadata","MRI3META.csv"))
mri_meta1 <- read.csv(here(rawdataDir,"MRI_metadata","MRIMETA.csv"))

load(here(cleandataDir,"dataForAnalysis_mixed.RData"))


names(basal)

aux <- by(dd,INDICES = dd$PTID, FUN=function(x){unique(x$tesla)},simplify = F)
basal$tesla <- as.numeric(aux)
head(basal)
table(basal$tesla)
res.pca <- PCA(basal[,15:26],scale=T,ncp=5,graph=F)

eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])

barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eigenvalues), eigenvalues[, 2], 
      type="b", pch=19, col = "red")
fviz_pca_var(res.pca)

library(factoextra)
fviz_pca_var(res.pca)
plotPC <- fviz_pca_ind(res.pca,label="none",habillage = as.factor(basal$tesla),addEllipses = T,
                       palette = c("#00AFBB", "#E7B800"))

df <- as.data.frame(res.pca$ind$coord[,1])
df <- data.frame(PTID=basal$PTID,df,basal$tesla)
names(df) <- c("PTID","PC1","tesla")
t.test(PC1~tesla,data=df)



####################
## HARMONITZAR AMB LONGCOMBAT
library(longCombat)
names(dd)

volVars <- sort(names(dd)[22:34])
dd$time_visit <- dd$time_visit/12
dd$ICV <- dd$ICV/1000

batchBoxplot(idvar='PTID', 
             batchvar='tesla', 
             feature=volVars, 
             formula='time_visit*DX.bl + AGE + PTGENDER + I(ICV/1000)',
             ranef='(1+time_visit|PTID)',
             data=dd)


dd_combat <- dd[,c("PTID","time_visit","tesla","AGE","DX.bl",
                   "PTGENDER","ICV",volVars)]
head(dd_combat)

dd_combat <- 
  longCombat(idvar='PTID', 
            timevar = 'time_visit',
            batchvar='tesla', 
            feature=volVars, 
            formula='time_visit*DX.bl + AGE + PTGENDER + ICV',
            ranef='(1|PTID)',
            data=dd_combat)

head(dd_combat)

dd2 <- dd_combat$data_combat
head(dd2)

basal.combat <- lapply(split(dd2,dd2$PTID), FUN="[",1,)
basal.combat <- do.call(rbind,basal.combat)
head(basal.combat)

res.pca.combat <- PCA(basal.combat[,4:15], scale=T, graph = F)
plotPC.combat <- fviz_pca_ind(res.pca.combat, label="none",
                              habillage = as.factor(basal$tesla),
                              addEllipses = T,
                              palette = c("#00AFBB", "#E7B800"))

df.combat <- as.data.frame(res.pca.combat$ind$coord[,1])
df.combat <- data.frame(PTID=basal$PTID,df.combat,basal.combat$tesla)
names(df.combat) <- c("PTID","PC1","tesla")
head(df.combat)
t.test(PC1~tesla,data=df.combat)

# Grafics
p1 <- ggplot(df, aes(x=as.factor(tesla), y=PC1,
                     color = as.factor(tesla))) +
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width = 0.2), alpha=0.2) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(color = "Tesla parameter",x="")

p2 <- p1 %+% df.combat


my_comparisons <- list( c("1.5", "3"))

p1 <- p1 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 15,method = "t.test")  

p2 <- p2 + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 15, method = "t.test")  


library(gridExtra)
compositPlot <- grid.arrange(plotPC,plotPC.combat,p1,p2)
ggsave(here("results","Others","longCombatHarmonization_tesla.png"),
       plot = compositPlot)
