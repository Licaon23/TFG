###############################################################################
### EXPLORATORY ANALYSIS                                                      #
###   ADNI patient paths diagram                                              #
###   Version 1                                                               #
###   December 2022                                                           #  
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###                                                                           #
###############################################################################
source("Rscripts/libraries.R")
load(here(cleandataDir,"ADNIPathsInfo.RData"))




#---------------------------------------------------------------------
# iii) ADNI project patient diagram: input, follow-ups and abandonment
#---------------------------------------------------------------------

nodes <- c("ADNI1","ADNIGO","ADNI2","ADNI3")
edges <- c("ADNI1-ADNIGO","ADNIGO-ADNI2", "ADNI1-ADNI2" ,"ADNI2-ADNI3")


# FUNCTION flow(): given an edge and the number of patients for a unique
# trajectory, get the number of patients that follow the project from
# one ADNI phase to another.
flow <- function(edge,trajectories){
  pos <- grep(pattern = edge, x=names(trajectories))
  return(sum(trajectories[pos]))
}

# FUNCTION out(): given an ADNI phase, get the number of patients that
# abandon the project at that point.
out <- function(node,trajectories){
  pos <- grep(paste0(".*(",node,"$)"),names(trajectories))
  return(sum(trajectories[pos]))
}

# FUNCTION input(): given an ADNI phase, get the number of new patients that
# enroll to ADNI project at that time.
input <- function(node,trajectories){
  pos <- grep(paste0("^",node,".*$"),names(trajectories))
  return(sum(trajectories[pos]))
}

# Get these values
inputs <- sapply(nodes,FUN=input,trajectories=unique.paths)
outputs <- sapply(nodes,FUN=out,trajectories=unique.paths)
flows <- sapply(edges, FUN=flow, trajectories=unique.paths)

# Get number of patients on each cohort
present <- sapply(split(patientCohorts,patientCohorts$COLPROT),
                  FUN=nrow)

# inputs
# outputs
# flows
# present

diagramValues <- list(new_enrolls =inputs, 
                      lost_pacients = outputs, 
                      followUp = flows, 
                      pacients_cohort = present)


# PLOT NETWORK DIAGRAM---------------------------------

#png(here("figures","ADNIdiagram.png"),width = 800,height = 480)
diagrama <-
  grViz("
digraph a_nice_graph {

label     = 'ADNI project evolution'
labelloc  =  t 
fontsize  = 30 
fontcolor = black

#NODES -----------------------------------
# Cohorts
node [fontname = Helvetica, shape=box]
a [label = '@@1-1']
b [label = '@@1-2']
c [label = '@@1-3']
d [label = '@@1-4']

#Inputs
node [shape=plaintext]
e [label = '@@2-1']
f [label = '@@2-2']
g [label = '@@2-3']
h [label = '@@2-4']

#Outputs
node [shape=plaintext]
i [label = '@@3-1']
j [label = '@@3-2']
k [label = '@@3-3']

# EDGES ---------------------------------------------
# edge definitions with the node IDs
a -> b [label = @@4-1]
b -> c [label = @@4-2]
a -> c [label = @@4-3]
c -> d [label = @@4-4]

# Inputs
edge [color=palegreen, arrowhead=vee,arrowtail=none]
e->a ; f->b; g->c; h->d

#Outputs
edge [color=Red, arrowhead=vee,arrowtail=none]
a->i ; b->j; c->k

}
#values
[1]: paste0(nodes,'\\n',present)
[2]: inputs
[3]: outputs
[4]: flows
")

diagrama %>% export_svg %>% charToRaw %>% 
  rsvg_png(here("results","Others","ADNIdiagram.png"))
rm(list=ls())
