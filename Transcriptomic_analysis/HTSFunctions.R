library(stringr)
library(edgeR)
concatFiles <- function(directory, mergeBy = 1, counts = "",log = "",...){
  
  cat(paste("#",directory, sep= ""), sep = "\n", log = log, append = T)
  names <- list.files(directory, full.names = T)
  dataFrame <- data.frame()
  
  
  for(i in 1:length(names)){
    file <- names[i]
    GSM <- str_extract_all( file, "GSM[0-9]*")[[1]] 
    isgz <- length(grep("\\.gz$", file)) == 1
    if (isgz){
      file <- gzfile(file, "rt")
      
    }
    if(i == 1){
      dataFrame <- read.csv(file,...)
      if(class(mergeBy) == "numeric"){
        mergeBy = colnames(dataFrame)[mergeBy]
      }
      if (counts == ""){
        counts <- colnames(dataFrame)[length(colnames(dataFrame))]
      }
      names(dataFrame)[names(dataFrame) == mergeBy ] <- "geneids"
      names(dataFrame)[names(dataFrame) == counts ] <- GSM
      if(isgz){close(file)}
    }else{
      df <- read.csv(file,...)
      dataFrame[[GSM]] <- df[[counts]][match(dataFrame[["geneids"]], df[[mergeBy]])]
      if(isgz){close(file)}
      
    }
  }
  return(dataFrame)
}


proccessHTS <- function(counts, designMatrix, contMatrix, log = "", direc, file ){
  
  d0 <- DGEList(counts)
  d0 <- calcNormFactors(d0)
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  designMatrix <- designMatrix[order(match(rownames(designMatrix), colnames(counts))), , drop = F]
  #reorder the design matrix so row ordering matches column ordering of annotatedGseset
  #limma will not look for rownames
  y <- voom(d, designMatrix, plot = F)

  fit <- lmFit(y, designMatrix)
  fit<-contrasts.fit(fit, contMatrix)
  fit<-eBayes(fit)
  
  for(coeff in colnames(contMatrix)){
    
    tt <-topTable(fit, coef = coeff, number= nrow(y), adjust.method="BH", sort.by="logFC")
    #tt$entrezIDs <- rownames(y$E)
    signGenes <- tt[ tt$adj.P.Val < 0.05, c("ID","logFC"), drop = F]
    write.table(tt, file = paste( direc,"/topTable/",file, "_top_table_", coeff ,sep = "")  )
    
    if(nrow(signGenes) == 0){
      cat(paste("No significant genes after p value adjustment found for", coeff), sep = "\n"
          , file = log, append = T)
      
    }else{
      getBINGSpaceGenes(signGenes = signGenes, probeToEntrezIDs = F, name = paste(file,coeff, sep = "_"),log = log, direc = direc, col = "ID" )
    }
  }
}
  
