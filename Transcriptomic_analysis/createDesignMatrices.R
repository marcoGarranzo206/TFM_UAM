# each design has its own function
# lots of code repetition for the special cases
# generateSpecialDesign.R has one for all the special types
# metaDataFile = "~/Marco_TFM/Marco_TFM/Data/metaData/GSE111073sample_info.txt"
# isPaired = T
# intercept = F
library(stringr)

generateFormula <- function(dfName = "", varnames, intercept = F){
  
  if(intercept){
    formula <- "~1+"
  }else{
    formula <- "~0+"
  }
  
  i = 0
  
  for(var in varnames){
    
    if( i == 0){
      formula <- paste(formula,  dfName,"$", var,sep = "")
    }else{
      formula <- paste(formula, " + ",dfName,"$", var, sep = "")
    }
    i <- i + 1
    
  }
  return(as.formula(formula))
}
generateStandardDesignMatrix <- function(meta, isPaired = F, intercept = F, ...){
 
  if(is.character(meta) & length(meta) == 1){
  varData <- read.table(meta, header = TRUE, ...)
  }else if (is.data.frame(meta)){
  varData <- meta
  }else{
  
	print("meta must be a file path or data frame")
	return  
  }
  
  if(isPaired | "Id" %in% colnames(varData)){
    varData$Id <- as.factor(varData$Id)
  }
  n <- grep("treatment", colnames(varData))
  if (n != 1){
    
    if (n == ncol(varData)){
      varData <- varData[ , rev(colnames(varData))]
    }else if (n == (ncol(varData) - 1)){
      varData <- data.frame(varData[,n,drop = F],   varData[,1:n-1,drop = F], varData[,n+1,drop = F])
    }else{
      varData <- data.frame(varData[,n,drop = F],   varData[,1:n-1,drop = F], varData[,n+1:ncol(varData),drop = F])
    }
     
    
  }
    
    
  assign("varDataG", varData, envir = .GlobalEnv) #needed because model.matrix looks in the global environment
  designMatrix <- model.matrix( generateFormula("varDataG", colnames(varData), intercept = intercept))
  colnames(designMatrix)
  rownames(designMatrix) <- rownames(varData)
  colnames(designMatrix) <- vapply(strsplit(colnames(designMatrix), "\\$" ), function(x) make.names(x[2]), FUN.VALUE = character(1))
  
  #generate design and contrast matrix of what we will call "simple" studies
  #That is:
  #studies, which, for each compound we only have one dose/time point
  #other categorical values can exist, and their effect will be removed with limma
  
  treatmentNames <- colnames(designMatrix)[grep("^treatment.*", colnames(designMatrix), ignore.case = T)]
  t <- treatmentNames[treatmentNames != "treatmentcontrol"]
  contMatrix <- matrix(0, nrow = ncol(designMatrix), ncol = length(t),dimnames = list(colnames(designMatrix), t) )
  contMatrix["treatmentcontrol", ] <- -1
  
  for(i in t){
    
    contMatrix[i, i] <- 1
    
    }

  ret <-list()
  ret$designMatrix <- designMatrix
  ret$contMatrix <- contMatrix
  return(ret)
}



# diff time points, each w/a control

generateDesignTimePoint <- function(metaDataFile, isPaired = F, intercept = F, varName= "time point",...){
  

  ## diff time points
  ## no 0 time point
  ## each time point has its own control
  varData <- read.table(metaDataFile, header = TRUE,...)
  varData <- data.frame(vapply(X=varData, FUN = function(x) make.names(x), FUN.VALUE = character(nrow(varData)))
                        ,row.names = row.names(varData))
  varName <- make.names(varName)
  varData$treatment <- str_replace_all(varData$treatment, pattern = "\\.", "_") 
  
  if(isPaired | "Id" %in% colnames(varData)){
    varData$Id <- as.factor(varData$Id)
  }
  treatmentVar <- paste(varData$treatment, varData[ , varName], sep = ".")
  newVarData <- data.frame( treatment = treatmentVar,  varData[,!colnames(varData)  %in% c("treatment", varName), drop = F],
                            stringsAsFactors = F)
  
  
  assign("newVarDataG", newVarData, envir = .GlobalEnv)
  designMatrix <- model.matrix( generateFormula("newVarDataG", colnames(newVarData), intercept = F) )
  colnames(designMatrix) <- vapply(strsplit(colnames(designMatrix), "\\$" ), function(x) make.names(x[2]), FUN.VALUE = character(1))
  rownames(designMatrix) <- rownames(newVarData)
  
  treatment <- unique(varData$treatment[varData$treatment != "control"])
  timepoints <- unique( varData[, varName])
  contMatrix <- matrix( data = 0, nrow = ncol(designMatrix), ncol = length(treatment)*length(timepoints))
  colnames(contMatrix) <- paste(rep(treatment, each = length(timepoints)), timepoints,sep = ".")
  rownames(contMatrix) <- colnames(designMatrix)
  
  for(co in colnames(contMatrix)){
    
    col_split <- strsplit(co, "\\.")[[1]]
    tp <- paste(col_split[2:length(col_split)], collapse = ".")
    
    control <- paste("treatmentcontrol", tp, sep = ".")
    treat <- paste("treatment", co, sep = "")
    
    contMatrix[treat, co] <- 1
    contMatrix[control, co] <- -1
  
  }
  
  ret <-list()
  ret$designMatrix <- designMatrix
  ret$contMatrix <- contMatrix
  return(ret)
  
}

generateDesign0asControl <- function(metaDataFile, isPaired = F, varName = "timepoint",...){
  
  #diff time points/conc, 0 is control
  #control can be shared across treatments or not
  varData <- read.table(metaDataFile, header = TRUE, ...)
  varData <- data.frame(vapply(X=varData, FUN = function(x) make.names(x), FUN.VALUE = character(nrow(varData)))
                        ,row.names = row.names(varData))
  varData$treatment <- str_replace_all(varData$treatment, pattern = "\\.", "_") 
  varName <- make.names(varName)
  
  if(isPaired | "Id" %in% colnames(varData)){
    varData$Id <- as.factor(varData$Id)
  }
  
  treatmentVar <- paste(varData$treatment, varData[ , varName], sep = ".")
  newVarData <- data.frame( treatment = treatmentVar,  varData[,!colnames(varData)  %in% c("treatment", varName), drop = F],
                           stringsAsFactors = F)
  
  assign("newVarDataG", newVarData, envir = .GlobalEnv)
  designMatrix <- model.matrix( generateFormula("newVarDataG", colnames(newVarData), intercept = F) )
  colnames(designMatrix) <- vapply(strsplit(colnames(designMatrix), "\\$" ), function(x) make.names(x[2]), FUN.VALUE = character(1))
  rownames(designMatrix) <- rownames(newVarData)
  
  treatmentCompounds <- unique(varData$treatment)
  print(treatmentCompounds)
  print("control" %in% treatmentCompounds)

  if("control" %in% treatmentCompounds){
  
	  treatmentCompounds <- treatmentCompounds[treatmentCompounds != "control"]
	  treatmentVars <- unique(varData[, varName][varData[, varName] != "X0"])
	  contMatrix <- matrix( data = 0, nrow = ncol(designMatrix), ncol = length(treatmentCompounds)*length(treatmentVars))
  	  colnames(contMatrix) <- paste(rep(treatmentCompounds, each = length(treatmentVars)), treatmentVars,sep = ".")
          rownames(contMatrix) <- colnames(designMatrix) #already has
	  contMatrix["treatmentcontrol.X0", ] <- -1

  	  for(co in colnames(contMatrix)){

   	 	treat <- paste("treatment", co, sep = "")
   	 	contMatrix[treat, co] <- 1
	 }
  }else{

  	#diff conc for each compound: redo for everything!
        contMatrix <- matrix( data = 0, nrow = ncol(designMatrix), ncol = length(unique(treatmentVar)) - length(treatmentCompounds))
  	colnames(contMatrix) <- unique(treatmentVar[-grep(pattern = "X0$", x = treatmentVar)])
 	rownames(contMatrix) <- colnames(designMatrix) #already has
  
  	for(co in colnames(contMatrix)){
    
    		control <- sub("\\..*", ".X0", co)
    		control <- paste("treatment", control, sep = "")
    		treat <- paste("treatment", co, sep = "")
    		contMatrix[treat, co] <- 1
   	 	contMatrix[control, co] <- -1
  	}
  }

  ret <-list()
  ret$designMatrix <- designMatrix
  ret$contMatrix <- contMatrix
  return(ret)
}


  
generateDesignBeforeAfter <- function(metaDataFile, isPaired = F, varName = "time point", ...){
  
  #before after studies
  varData <- read.table(metaDataFile, header = TRUE, ...)
  varData <- data.frame(vapply(X=varData, FUN = function(x) make.names(x), FUN.VALUE = character(nrow(varData)))
                        ,row.names = row.names(varData), stringsAsFactors = F)
  
  #in case some treatment has a dot due to bad naming, change to _, otherwise will interfere
  #with value assingmen in contMatrix
  varData$treatment <- str_replace_all(varData$treatment, pattern = "\\.", "_") 
  varName <- make.names(varName)
  
  
  if(isPaired | "Id" %in% colnames(varData)){
    varData$Id <- as.factor(varData$Id)
  }
  
  treatmentVar <- paste(varData$treatment, varData[ , varName], sep = ".")
  newVarData <- data.frame( treatment = treatmentVar,  varData[,!colnames(varData)  %in% c("treatment", varName), drop = F],
                            stringsAsFactors = F)
  
  assign("newVarDataG", newVarData, envir = .GlobalEnv)
  designMatrix <- model.matrix( generateFormula("newVarDataG", colnames(newVarData), intercept = F) )
  colnames(designMatrix) <- vapply(strsplit(colnames(designMatrix), "\\$" ), function(x) make.names(x[2]), FUN.VALUE = character(1))
  rownames(designMatrix) <- rownames(newVarData)
  
  
  treatment <- unique(newVarData$treatment)
  treatment <- treatment[-c(grep(pattern = "^control", x = treatment), grep(pattern = "X0$", x = treatment))]
  
  timepoints <- unique( varData[, varName])
  timepoints <- timepoints[which(timepoints != "X0")]
  
  contMatrix <- matrix( data = 0, nrow = ncol(designMatrix), ncol = length(treatment))
  colnames(contMatrix) <- treatment
  rownames(contMatrix) <- colnames(designMatrix)
  
  contMatrix["treatmentcontrol.X0", ] <- 1 #here control means placebo
  
  for(co in colnames(contMatrix)){
    
    col_split <- strsplit(co, "\\.")[[1]]
    tr <- col_split[1]
    tp <- paste(col_split[2:length(col_split)], collapse = ".")
    
    control <- paste("treatment", tr, ".X0", sep = "") #here control refers to time point 0
    case <- paste("treatment", co, sep = "")
    contMatrix[case, co] <- 1

    if(control %in% row.names(contMatrix)){
	#sometimes we dont have a time point 0 for a given condition
	#which os not ideal
	#but here we just treat those case as if it were a diff time point 
	#study, with the benefit that control has its time
	# effect 0 subtracted
    	contMatrix[control, co] <- -1
    }

    contMatrix[paste("treatmentcontrol", tp, sep = "."), co] <- -1
  }
  
  ret <-list()
  ret$designMatrix <- designMatrix
  ret$contMatrix <- contMatrix
  
  return(ret)
}

contMatrixFromDesign <- function(designMatrix, varsep = "_"){
  
  treatmentColumns <- colnames(designMatrix)[grep(pattern = paste("^treatment", varsep, sep = ""), colnames(designMatrix))]
  t <- treatmentColumns[grep("control$", treatmentColumns, invert = T)]
  contMatrix <- matrix(rep(0,times = ncol(designMatrix)*length(t)), nrow =  ncol(designMatrix) )
  colnames(contMatrix) <- t
  rownames(contMatrix) <- colnames(designMatrix)
  colTreat <- paste("treatment",varsep,"control", sep = "")
  if( colTreat %in% colnames(designMatrix)){
    contMatrix[colTreat, ] <- -1
  }

  for(condition in t){
    #print(condition)
    contMatrix[condition, condition] <- 1
  }
  
  ret <- list()
  ret$designMatrix <- designMatrix
  ret$contMatrix <- contMatrix
  return(ret)
}





