rm(list=ls())

# Paths --------------------------------------
library(here) # work with relative paths
source(here("Rscripts","paths.R")) # load directory paths

# Graphics ----------------------------------------------
library(ggplot2) # grpahics
library(ggstatsplot) # add-ons to ggplot
#library(RColorBrewer) # colors ggplot2
library(ggtext) # texts to ggplot
library(ggResidpanel) # residuals diagnostic plots
library(GGally)
library(rstatix) # pairwise test
library(ggpubr) # pvalue anotation
library(ggrepel)

# Flux diagrams----------------------------------------
library(DiagrammeR) # build a box diagram
library(DiagrammeRsvg)
library(rsvg)

# Descriptive analysis ----------------------------------
library(doBy) # perfomr group statistics
library(corrplot)
library(Hmisc) # smean.cl.normal

# Latex and HTML tables ----------------------------------
library(kableExtra)

# Export Excel
library(openxlsx)

# Modelling ---------------------------------------------
library(lme4) # linear mixed effect models
library(mixedup) # handy package to extract info from lmer objets
library(lmerTest) # add-on to lme4. It calculates p-values using 
                  # satterthwaite method for degrees of freedom

library(emmeans) # Estimated marginal means and trends.
library(multcomp)
