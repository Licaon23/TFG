########################
#Auxiliary functions ###
########################
#-----------------------------------------------------------------------
fix_nobl <- function(dd){
  aux <- split(dd,droplevels(dd$PTID))
  res <- lapply(aux,function(x){
    x$enrolDate <- x$measureDate[1] 
    x$AGE <- round(x$AGE + x$timelapse[1]/12,1)
    return(x[,c("PTID","AGE","enrolDate","measureDate")])
  })
  res <- do.call(rbind,res)
  rownames(res) <- NULL
  return(res)
}

#-----------------------------------------------------------------------
get_repeatedViscode <- function(label,dd){
  splitVol <- split(dd,dd$PTID)
  target <- names(which(sapply(splitVol,
                               FUN=function(x) table(x$VISCODE2)[[label]]>1)))
  attributes(target)$VISCODE2 <- label
  return(target)
}

#-----------------------------------------------------------------------
dropCases <- function(ptid,dd){
  label <- attributes(ptid)$VISCODE2
  
  duplicateMeasures <- sapply(ptid,FUN=function(x){
    which(dd$PTID == x & dd$VISCODE2==label)
  })
  
  dropMeasures <- duplicateMeasures[1,]
  return(dropMeasures)
}

#-----------------------------------------------------------------------
get_viscode <- function(timelapse){
  breaks <- c(-Inf,1.4,4.4,seq(10.4,130.4,by=6),138,150,+Inf)
  labels <- levels(demo$VISCODE2)
  viscode <- cut(timelapse, 
                 breaks=breaks,
                 labels=labels,
                 include.lowest=T,
                 right=F)
  return(viscode)
}


#-----------------------------------------------------------------------
get_overallVolume <- function(subfieldName,dd){
  regex <- paste0("^(left|right)_",subfieldName,"$")
  vol <- grepl(regex,colnames(dd))
  vol <- dd[vol]
  volumeOverall <- do.call("+",vol)
  
  return(volumeOverall)
}


#-----------------------------------------------------------------------
f <- function(label,dd){
  target <- names(which(sapply(split(dd,dd$PTID),
                               FUN=function(x) table(x$VISCODE2)[[label]]>1)))
  length(target)
}
#-----------------------------------------------------------------------

set_attributes <- function(x,dataset,codebook){
  
  var <- dataset[,codebook[x,"Variable"]]
  attributes(var)$Description <- codebook[x,"description"]
  attributes(var)$Label <- codebook[x,"Label"]
  
  return(var)
}
