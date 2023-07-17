###############################################################################
### EXPLORATORY DATA ANALYSIS                                                 #
###   Descriptive statistics functions                                        #
###   Version 1                                                               #
###   January 2023                                                            #
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed:                                                        #
###     -                                                                       #
###############################################################################


#====================================================================
# AUXILIARY FUNCTION
# catSummary
#   INPUTS
#     - x: categorical variable name
#     - xname: label 
#     - dig: number of decimals
#     - data: dataset
#   OUTPUT:
#     Formatted count and frequency n(p%)
#====================================================================
catSummary <- function(x, xname, dig=1, data){
  
  if(is.null(xname)){xname <- x}
  # Extraure variable factor del dataset
  x <- data[[x]]
  
  
  # Compatges i percentatges
  count <- table(x)
  freq <- prop.table(count) * 100
  
  # Format resultat
  count_formatted <- format(as.numeric(count), justify="right")
  
  freq_formatted <- sprintf(paste0("%.",dig,"f"),freq)
  freq_formatted <- format(freq_formatted,justify = "right")
  
  # Result
  res <- paste0(count_formatted," (",freq_formatted,"%)")
  
  # Dataframe resultat
  
  if(length(levels(x))!=2){
    res <- data.frame(variable=paste0("",levels(x)),
                      summary = res)
    res <- rbind(c(xname,rep("",dim(res)[2]-1)),res)
    
  }else{
    res <- res[-1]
    xname <- paste0(xname," (",levels(x)[-1],")")
    res <- data.frame(variable = xname,
                      summary = res)
  }
  
  rownames(res) <- NULL
  return(res)
}


#====================================================================
# AUXILIARY FUNCTION
# numSummary
#   INPUTS
#     - x: numerical variable
#     - xname: label 
#     - dig: number of decimals
#     - data: dataset
#   OUTPUT:
#     Formatted output: mean(sd) for normal variables
#                       median [Q1,Q3] for non-normal variables
#====================================================================
numSummary <- function(x, xname, dig=1, data){
  if(is.null(xname)){xname <- x}
  
  # Extraure variable factor del dataset
  x <- data[[x]]
  
  # Estadístics
  m <- mean(x, na.rm=T)
  s <- sd(x, na.rm=T)
  q <- quantile(x, probs=c(0.25,0.5,0.75), na.rm=T)
  
  
  # Format resultat
  if(min(x,na.rm = T)>100){dig<-0}
  if(min(x,na.rm = T)<10){dig <- 2}
  
  normalityTest <- shapiro.test(x)
  
  if(normalityTest$p.value < 0.05){
    q <- sprintf(paste0("%.",dig,"f"),q)
    res <- paste0(q[2]," [",q[1],", ",q[3],"]")
    
  }else{
    res <- sprintf(paste0("%.",dig,"f"),c(m,s))
    res <- paste0(res[1]," (",res[2],")")
  }
  
  res <- data.frame(variable = xname,
                    summary = res)
  
  return(res)
}


#====================================================================
# AUXILIARY FUNCTION
# summaryWrapper
#   INPUTS
#     - x: variable
#     - xname: label 
#     - dig: number of decimals
#     - data: dataset
#   OUTPUT:
#     Depending on class variable, execute catSummary() or numSummary()
#====================================================================
summaryWrapper <- function(x, xname=NULL, dig=1, data){
  
  #Numeric
  if(is.numeric(data[[x]])){res <- numSummary(x,xname,dig,data)}
  
  #Factor
  if(is.factor(data[[x]])){res <- catSummary(x,xname,dig,data)}
  
  return(res)
}



#====================================================================
# MAIN FUNCTION
# descSummary
#   INPUTS
#     - data: dataset
#     - vars: variable names (vector)
#     - printNames: labels 
#     - dig: number of decimals
#
#   OUTPUT:
#     Dataframe with both variable name and descriptive statistics
#====================================================================
descSummary <- function(data, vars, printNames=NULL, dig=1){
  if(is.null(printNames)){printNames <- vars}
  
  if(length(printNames)!=length(vars)){
    stop("Incorrect input")
  }
  
  
  varNames <- data.frame(x=vars, xnames=printNames)
  res <- lapply(varNames$x, 
                function(x){summaryWrapper(x=x,
                            xname = varNames[which(varNames$x==x),"xnames"],
                            dig=dig,
                            data=data)})
  
  res <- do.call(rbind,res)
  return(res)
}


#====================================================================
# MAIN FUNCTION
# groupContrast
#   INPUTS
#     - x: variable name
#     - by: group factor to be split by 
#     - data: dataset

#   OUTPUT:
#     pvalue of the corresponding test
#====================================================================
groupContrast <- function(x,by,data){
  
  #Factor
  if(is.factor(data[[x]])){
    tab <- table(data[[x]],data[[by]])
    lowCounts <- any(tab<5)
    test <- rbind(c("fisher.test",TRUE),
                  c("chisq.test",FALSE))[c(lowCounts,!lowCounts),]
    pval <- eval(parse(text=paste0(test[1],"(tab,simulate.p.value=",
                                   test[2],")")))$p.value
    
    pval <- ifelse(pval<0.001,"<0.001",sprintf("%.3f",pval))
  }
  
  #Numeric
  if(is.numeric(data[[x]])){
    multiGroup <- length(unique(data[[by]])) > 2 
    normTest <-sapply(split(data[[x]],data[[by]]),
                      function(x) shapiro.test(x)$p.value)
    nonNormal <- min(normTest) < 0.05
    testOptions <- c("wilcox.test","t.test","kruskal.test","aov")
    useTest <- testOptions[c(nonNormal&!multiGroup,
                             !nonNormal&!multiGroup,
                             nonNormal&multiGroup,
                             !nonNormal&multiGroup)]
    expr <- paste(x,"~",by)
    
    
    test <- eval(parse(text=paste0(useTest,"(as.formula(",expr,"),
                                   data=data)")))
    pval <- ifelse(useTest=="aov", pval<-anova(test)$`Pr(>F)`, test$p.value)
    
    #Format
    pval <- ifelse(pval<0.001,"<0.001",sprintf("%.3f",pval))
    #    pval <- sprintf("%.4f",pval)
    
  }
  
  return(pval)
}



