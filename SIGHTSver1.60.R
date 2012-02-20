# use try command to assess whether MASS is installed.
x <- try(library(MASS), silent=TRUE)

if(length(grep("Error", x[1])) > 0){
   install.packages("MASS")
   library(MASS)
}
# use try command to assess whether lattice is installed.
x <- try(library(lattice), silent=TRUE)

if(length(grep("Error", x[1])) > 0){
   install.packages("lattice")
   library(lattice)
}

# use try command to assess whether qvalue is installed.
x <- try(library(qvalue), silent=TRUE)

if(length(grep("Error", x[1])) > 0){
    source("http://www.bioconductor.org/biocLite.R")
    biocLite("qvalue")
}

#
# Apply normalization method to data
#
# method - type of normalization method
# dat - matrix of data. rows are compounds/genes and columns are plates.
# trimm - trim factor used in trimmed polish
# doW - if T apply individual well correction
# params - optional spatial bias vector. if NULL then the median of
#          each row of data is used.
# plates - column index of plates to use to calculate median of rows for params
# repIndex - vector of labels denoting replicate groups
# numRows - number of rows in plate
# numCols - number of columns in plate
#
applyNormalization <- function(dat, trimmy=trimm, params=NULL, plates=spatialBiasEstimatePlates,
                  repIndex=replicateIndex, rows=numRows, cols=numCols, methh=method) {



   if(methh=="Z") {
      mat <- doZ(dat)
   }
   else if (methh=="robZ")   {
      mat <- doRobZ(dat, doW=F)
   }
   else if (methh=="robZW") {
      mat <- doRobZ(dat, doW=T)
   }
   else if(methh=="SPAWN") {
      mat <- doTrimmedPolish(dat, doW=F)
   }
   else if(methh=="SPAWNW") {
      mat <- doTrimmedPolish(dat, doW=T)
   }
   else if(methh=="R") {
      mat <- doRobustRegressionModel(dat, doW=F)
   }
   else if(methh=="RW") {
      mat <- doRobustRegressionModel(dat, doW=T)
   }
   else if(methh=="LMF") {
      mat <- doBarysh(dat)
   }
   else if(methh=="Well Correction") {
      mat <- doMak(dat)
   }
   else if(methh=="Median Filter") {
      mat <- doMedFilter(dat)
   }
    else if(methh=="Loess") {
      mat <- doLoess(dat)
   }
   else {
       output <- cat(methh, "is not a valid normalization method\n")
       message(output)
   }

   return(mat)

} # end applyNormalization



#
# Apply robust Z method to data
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
# params - optional spatial bias vector. if NULL then the median of
#          each row of data is used.
# plates - column index of plates to use to calculate median of rows for params
#
doRobZ <- function(myData, doW=F, params=NULL, plates=spatialBiasEstimatePlates) {

   myRobZ <- function(x) {
      return( (x - median(x, na.rm=T))/mad(x, na.rm=T))
   }

   outMat <- apply(myData, 2, myRobZ)

   # apply individual well correction if required
   if(doW) {

       # calculate spatial bias if not supplied
       if(is.null(params)) {
          params <- apply(outMat[,plates], 1, median, na.rm=T)
       }
       
       message("Completed robust Z score with well normalization")
       return(outMat - params)
   } # end if doW
           
   dimnames(outMat) <- dimnames(myData)
   
   message("Completed robust Z score without well normalization")

   return(outMat)
}



#
# Apply  Z method to data
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
# params - optional spatial bias vector. if NULL then the median of
#          each row of data is used.
# plates - column index of plates to use to calculate median of rows for params
#
doZ <- function(myData, doW=F, params=NULL, plates=1:dim(dat)[[2]]) {

   myZ <- function(x) {
      return( (x - mean(x, na.rm=T))/sd(x, na.rm=T))
   }

   outMat <- apply(myData, 2, myZ)

   # apply individual well correction if required
   if(doW) {

       # calculate spatial bias if not supplied
       if(is.null(params)) {
          params <- apply(outMat[,plates], 1, median, na.rm=T)
       }

       return(outMat - params)
   } # end if doW

   dimnames(outMat) <- dimnames(myData)

   if(doW) {
      message("Completed Z score with well normalization")
   }
   else {
      message("Completed Z score without well normalization")
   }

   return(outMat)
}


#
# Apply trimmed polish to data.
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
# trimm - trim factor used in trimmed polish
# doW - if T apply individual well correction
# params - optional spatial bias vector. if NULL then the median of
#          each row of data is used.
# plates - column index of plates to use to calculate median of rows for params
# numRows - number of rows in plate
#
doTrimmedPolish <- function(myData, trimmy=trimm, doW=F, params=NULL, plates=spatialBiasEstimatePlates,
                  rows=numRows) {

   #source(paste(PROJECT.ROOT,"utilities//scripts//HTS//trimmedPolish.r",
   #             sep="//"))

   len <- dim(myData)[[2]]

   # calculate trimmed mean polish
   tpDat <- NULL
   for(i in 1:len) {
     tp <- c(trimpolish(matrix(myData[,i], nrow=rows, byrow=F),trace.iter=F,
                          trim=trimmy)$residuals)
      tpDat <- cbind(tpDat, tp/mad(tp, na.rm=T))
   }
   
   dimnames(tpDat) <- dimnames(myData)

   # apply individual well correction if required
   if(doW) {

      # calculate  if not supplied
       if(is.null(params)) {
          params <- apply(tpDat[,plates], 1, median, na.rm=T)
       }
       tmp <- tpDat - params
       message("Completed Spatial Polish, trim ", trimmy, " with well normalization")
       return(tmp)
   } # end if doW

   message("Completed Spatial Polish, trim ", trimmy, " without well normalization")

   return(tpDat)
}


#
# apply robust robust regression model separately to each plate.
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
# doW - if T then apply individual well correction (default - all plates)
# params - optional covariate which models the spatial bias. if NULL then the median of
#          each row of data is used.
# plates - vector of column indexes indicating which plates to be used in estimating the spatial
#           bias.
# replicateIndex - vector of labels indicating replicate group.  Each index in vector
#                   matches the corresponding column of mat.
#
doRobustRegressionModel <- function(myData, doW=FALSE, rows= numRows, cols=numCols,
                                  params=NULL, plates=spatialBiasEstimatePlates, repIndex=replicateIndex) {

    library(MASS)      # contains rlm

    len <- dim(myData)[[2]]
    rowIndex <- rep(1:8, cols)
    colIndex <- rep(1:10, each=rows)

    # calculate covariate if not supplied
    if(is.null(params)) {
        params <- apply(myData[,plates], 1, median, na.rm=T)
    }

    # apply rlm model to x with  coefficient params including rows
    # and columns.  Normalize by dividing by variance estimate.
    rlmRCME <- function(x, rows, cols, params) {
       ret <- rlm(x~rows+cols+params)
       return(ret$residuals/ret$s)
    }

    # apply rlm model to x with coefficient params rows and cols.
    # Normalize by dividing by variance estimate.
    rlmRC <- function(x, rows, cols) {
       ret <- rlm(x~rows+cols)
       return(ret$residuals/ret$s)
    }

    outMat <- NULL
    for(i in 1:len) {
       tmpDatt <- myData[,i]
       lenn <- length(tmpDatt)
       naInd <- is.na(tmpDatt)
       if(doW) {
             ret <- rlmRCME(tmpDatt, as.factor(rowIndex),
                            as.factor(colIndex), params)
       }
       else {
          ret <- rlmRC(tmpDatt, as.factor(rowIndex), as.factor(colIndex))
       }
       tmpRet <- rep(NA, lenn)
       tmpRet[!naInd] <- ret
       outMat <- cbind(outMat, tmpRet)
   }

   dimnames(outMat) <- dimnames(myData)

   if(doW) {
      message("Completed Robust Regression model with well normalization")
   }
   else {
      message("Completed Robust Regression model without well normalization")
   }
    return(outMat)
}



#
# Apply loess normalization to data.  ref Baryshnikova, 2011
# NA removal is based on 1 NA in data due to CMBA titration series
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
doLoess <- function(dat) {

   coll <- 1:numCols
   roww <- 1:numRows

   # do loess by plate
   loessMe <- function(plate, nummers=numRows) {
   
         coll <- 1:numCols
         roww <- 1:numRows

       tmp <- matrix(plate, nrow=nummers, byrow=F)
       numCols <- dim(tmp)[[2]]

       rloess <- cloess <- NULL

       # get row loess fit
       for(i in 1:nummers) {
           nas <- rep(NA, numCols)
           naInd <- coll[!is.na(tmp[i,])]
           x <- tmp[i,]
           tmp1 <- loess(x~coll)$fitted
           nas[naInd] <- tmp1
           rloess <- rbind(rloess,nas)
       }

       # get column loess fit
       for(i in 1:numCols) {
           nas <- rep(NA, numRows)
           naInd <- roww[!is.na(tmp[,i])]  
           x <- tmp[,i]
           tmp1 <- loess(x~roww)$fitted
           nas[naInd] <- tmp1
           cloess <- cbind(cloess,nas)
         }

       # get row and column averages of smoothers
       rmeans <- apply(rloess, 1, mean, na.rm=T)
       cmeans <- apply(cloess, 2, mean, na.rm=T)

       # calculate normalized scores in inefficient fashion
       outMat <- matrix(NA, nrow=nummers, ncol=numCols)
       for(i in 1:nummers) {
          for(j in 1:numCols) {
             #cat(i, " ", j, "\n")
             adj <- (rmeans[i]/rloess[i,j])*(cmeans[j]/cloess[i,j])
             outMat[i,j] <- tmp[i,j]*adj
          } # end for j

       } # end for i

       return(c(outMat))
   }

   outie <- apply(dat, 2, loessMe, numRows)

   myRobZ <- function(x) {
      return( (x - median(x, na.rm=T))/mad(x, na.rm=T))
   }

   message("Completed Loess")
   return(apply(outie,2, myRobZ))

} # end doLoess




#
# Apply Baryshnikova median filter normalization, med filter of window size 1, then
# average filter of window size 1.
# # WARNING - neigbhourhood matrices assume plate sizes of 8x10, 16x24, 32x40
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
#
BaryshFilter <- function(dat) {

   currentDir <- getwd()

   # determine size of plate by number of wells
   if(dim(dat)[[1]] <= 80) {
       plate=1
   }
   else if(dim(dat)[[1]] <= 384) {
       plate=2
   }
   else {
       plate=3
   }

   # Get neighbours for large plate (windows of 3 and 4)
   if(plate==3) {
         buddies1 <- createNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=3)
         buddies2 <- createNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=4)
    }
    # get neighbours for small and medium plates (windows of 2)
    else if(plate==2) {
       buddies1 <- buddies2 <- createNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=2)
    }
   # get neighbours for small and medium plates (windows of 1)
    else {
      buddies1 <- buddies2 <- createNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=1)
    }

    # assumes unlogged data
   doBMedFilter <- function(plate, buds=buddies1, useMed=T) {

      #print(dim(buds))
       outVec <- NULL

      for(i in 1:length(plate)) {
          #print(i)
          hood <-  buds[i,!is.na(buds[i,])]
          medd <- median(plate[hood], na.rm=T)
          outVec <- c(outVec, log(plate[i]/medd))
      }
      return(outVec)
   }

   # input is logged residuals from first filter
   doBAvgFilter <- function(plate, buds=buddies2) {

      #print(dim(buds))
       outVec <- NULL

      for(i in 1:length(plate)) {
          #print(i)
          hood <-  buds[i,!is.na(buds[i,])]
          medd <- mean(plate[hood], na.rm=T)
          outVec <- c(outVec, medd)
      }
      return(outVec)
   }

   doFilter <- function(plate, filt) {

      #print(dim(buds))
       outVec <- NULL

      for(i in 1:length(plate)) {
          #print(i)
          outVec <- c(outVec, plate[i]/exp(filt[i]))
      }
      return(outVec)
   }


   # calculate spatial filter
   filt1 <- apply(dat, 2, doBMedFilter)
   filt2 <- apply(filt1, 2, doBAvgFilter)

   # apply spatial filter
   outie <- NULL
   for(i in 1:dim(dat)[[2]]) {
       outie <- cbind(outie, doFilter(dat[,i], filt2[,i]))
   }

   return(outie)
}





#
# Apply Baryshnikov normalization to data.   median filter then loess
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
#
doBarysh <- function(dat) {

   tmp1 <- BaryshFilter(dat)
   tmp2 <- doLoess(tmp1)
   #tmp2 <- doRobZ(tmp1)
   message("Completed LMF")
   return(tmp2)

} # end doBarysh

 #
# Apply Makaranekov normalization to data.  Fit least squares across well location
# and subtract then Z score across well location.
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
#
doMak <- function(dat) {

   makAttack <- function(x) {
         indVar <- 1:length(x)
         #print(length(x)
         ret <- lm(y~ind1, data=data.frame(y=x, ind1=indVar))  # calculate linear fit
         tmp <- ret$residuals
         #plot(tmp)
         out <- (tmp - mean(tmp))/sd(tmp)
         return(out)
   }

   myZ <- function(x) {
      return( (x - mean(x, na.rm=T))/sd(x, na.rm=T))
   }

   zDat <- apply(dat,2,myZ)
   #par(mfrow=c(2,2))
   outie <- t(apply(zDat, 1, makAttack))
   message("Completed Well Correction")
   return(outie)

}



#
# Apply median filter normalization to data.  ref Bushway, 2011
#
# dat - matrix of data. rows are compounds/genes and columns are plates.
# params - optional spatial bias vector. if NULL then the median of
#          each row of data is used.
# plates - column index of plates to use to calculate median of rows for params
# doSeq - if T apply initial row median filter than standard filter, else just apply
#          standard filter
#
doMedFilter <- function(dat, doSeq=T) {

  # numCols <- dim(dat)[[1]]/numRows

   # determine size of plate by number of wells
   if(dim(dat)[[1]] <= 80) {
       plate=1
   }
   else if(dim(dat)[[1]] <= 384) {
       plate=2
   }
   else {
       plate=3
   }

   # Get neighbours for large plate (windows of 3 and 4)
   if(plate==3) {
         buddies <- createNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=3)
         rowBuddies <- createSequentialNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=3)
    }
    # get neighbours for small and medium plates (windows of 2)
    else if(plate==2) {
       buddies <- createNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=2)
       rowBuddies <- createSequentialNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=2)

    }
   # get neighbours for small and medium plates (windows of 1)
    else {
      buddies <- createNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=1)
      rowBuddies <- createSequentialNeighbourhoodMatrix(nrows=numRows, ncols=numCols, wind=2)
    }


   #cat("neighbs", dim(buddies), "\n")
   applyMedFilter <- function(plate, buds=buddies) {

      #print(dim(buds))
      grandMed <- median(plate, na.rm=T)
      outVec <- NULL

      for(i in 1:length(plate)) {
          #print(i)
          hood <-  buds[i,!is.na(buds[i,])]
          medd <- median(plate[hood], na.rm=T)
          outVec <- c(outVec, plate[i]*(grandMed/medd))
      }
      return(outVec)
   }  # end applyMedFilter

   mat1 <- dat
   if(doSeq) {
      mat1 <- apply(dat, 2, applyMedFilter, rowBuddies)
   }
   mat2 <- apply(mat1, 2, applyMedFilter, buddies)

   myRobZ <- function(x) {
      return( (x - median(x, na.rm=T))/mad(x, na.rm=T))
   }

   message("Completed Median Filter")
   return(apply(mat2,2,myRobZ))

} # end doMedFilter



#
# Apply t-test to data
# mat - matrix of data.  each column in the matrix is a single HTS plate and the rows
#       are the well location
# side - paramter describing type of t-test. (two.sided, less, greater)
# muu - number indicating true value of mean
# replicateIndex - vector of labels indicating replicate group.  Each index in vector
#                   matches the corresponding column of mat.
#
applyTtest <- function(mat, side=testType, muu=0, repIndex=replicateIndex) {

   myT <- function(x, side, muu) {
      tmp <- t.test(x, alternative=side, mu=muu)
      return(c( tmp$statistic, mean(x, na.rm=T)-muu, (1/sqrt(tmp$parameter+1))*sd(x, na.rm=T),  tmp$parameter, tmp$p.value))
   }


    outMat <- NULL
    for(i in unique(repIndex)) {
       tmpDatt <- mat[,repIndex==i]
       lenn <- dim(tmpDatt)[[2]]
       colNames <-  paste(i, c(":T-statistic", ":Mean_Difference", ":Standard_Error", ":DegreesOfFreedom", ":p-value"), sep="")
       tmpRes <- t(apply(tmpDatt, 1, myT, side, muu))
       dimnames(tmpRes)[[2]] <- colNames
       outMat <- cbind(outMat, tmpRes)
    }

   #dimnames(outMat)[[1]] <- dimnames(mat)[[1]]
  
    if(!is.null(dimnames(mat)[[1]])) {
       dimnames(outMat)[[1]] <- dimnames(mat)[[1]]
    }

   message("Completed T test")
   return(outMat)
}



#
# Apply RVM test to data
# mat - matrix of data.  each column in the matrix is a single HTS plate and the rows
#       are the well location
# side - paramter describing type of test. (two.sided, less, greater)
# muu - number indicating true value of mean
# replicateIndex - vector of labels indicating replicate group.  Each index in vector
#                   matches the corresponding column of mat.
#
applyRVM <- function(mat, side= testType, repIndex=replicateIndex) {

   numRow <- dim(mat)[[1]]

   # apply RVM test and return p value
   myRVM <- function(x, side="less") {
      if(side=="two.sided") {
          res <- RVT1(x)
      }
      else {
         res <- RVT1Sided(x, side=side)
      }
      return(c(res$v, res$mn, sqrt(1/res$vvr), res$n-1+2*res$a, res$vp))
      #return(list(pval=res$vp, ab=cbind(res$a[1], res$b[1])))
   }

    newDat <- NULL
    for(i in unique(repIndex)) {
        newDat <- rbind(newDat, mat[,repIndex==i])
    }
    tmpRes <-  matrix(myRVM(newDat, side=side), ncol=5, byrow=F)

    rowIndex <- seq(1, dim(tmpRes)[[1]], numRow)

    outMat <- NULL
    startNum <- 1
    for(i in 1:(length(unique(repIndex)))) {
       tmp <- tmpRes[startNum:(startNum+numRow-1),]
       colNames <-  paste(i, c(":RVM T-statistic", ":Mean_Difference", ":Standard_Error", ":DegreesOfFreedom", ":p-value"), sep="")
       dimnames(tmp)[[2]] <- colNames
       outMat <- cbind(outMat, tmp)
       startNum <- startNum +numRow
    }

    if(!is.null(dimnames(mat)[[1]])) {
       dimnames(outMat)[[1]] <- dimnames(mat)[[1]]
    }
    message("Completed RVM statistical test")
    return(outMat)
}



#
# Apply Storey FDR to data
# mat - matrix of p value data.  each column in the matrix is a single HTS plate and the rows
#       are the well location
# FDRMethod - use either smoother or bootstrap to estimate null distribution.
#
applyFDR <- function(mat, method=FDRMethod) {

   myQ <- function(x) {
       return(qvalue(x)$qvalues)
   }

   pvalInd <- seq(5, dim(mat)[[2]],5)
   matData <- c(mat[,pvalInd])
   tmpRet <- qvalue(matData, pi0.method=method)$qvalues
   ret <- matrix(tmpRet, ncol=length(pvalInd), byrow=F)
   colNames <-  paste(dimnames(mat[,pvalInd])[[2]], "FDR", sep="")
   colNames <- sub("p", "q", colNames)
   dimnames(ret) <- list(dimnames(mat)[[1]], colNames)
   message("Completed FDR")
   return(ret)
}

#
# Setter for normalization method
#
setNormalizationMethod <- function(meth=NULL) {
    method <<- meth
    message("Set normalization method: ", method, "\n")
}
#
# Setter for replicateIndex
#
setReplicateIndex <- function(plateIndex=NULL) {
    replicateIndex <<- plateIndex
    outMessage <- NULL
    for(i in 1:length(replicateIndex)) {
        outMessage <- paste(outMessage, replicateIndex[i], sep=" ")
    }
    message("Set replicateIndex: ", outMessage, "\n")
}

#
# Setter for numRows
#
setRowNumber <- function(numberOfRows=NULL) {
    numRows <<- numberOfRows
    message("Set number of rows: ", numRows, "\n")
}

#
# Setter for numCols
#
setColumnNumber <- function(numberOfColumns=NULL) {
    numCols <<- numberOfColumns
    message("Set number of columns: ", numCols, "\n")
}


#
# Setter for spatialBiasEstimatePlates
#
setSpatialBiasEstimatePlates <- function(biasEstimatePlates=NULL) {
    spatialBiasEstimatePlates <<- biasEstimatePlates
    outMessage <- NULL
    for(i in 1:length(spatialBiasEstimatePlates)) {
        outMessage <- paste(outMessage, spatialBiasEstimatePlates[i], sep=" ")
    }

    message("Set spatial bias estimate plates: ", outMessage, "\n")
}

#
# Setter for trimm
#
setTrim <- function(trim=NULL) {
    trimm <<- trim
    message("Set trim: ", trimm, "\n")
}

#
# Setter for testType
#
setTestType <- function(testy="two.sided") {
    testType <<- testy
    message("Set TestType: ", testType, "\n")
}

#
# Setter for FDRmethod
#
setFDRMethod <- function(method="smoother") {
    FDRMethod <<- method
    message("Set FDRMethod: ", method, "\n")
}
########################################################################################
###################################   Graphs   #########################################
########################################################################################


# construct a boxplot in order
#
# mat - matrix of data with each column a single HTS plate
#
boxPlot <- function(mat) {
   boxplot(data.frame(mat), names=dimnames(mat)[[2]], las=2, pars=list(las=2))
}


#
# Scatter plot of all pairwise combinations of replicates
#
#  scatDat - matrix of data with each column a single HTS plate
#
scatterPlot <- function(scatDat) {

    len <- dim(scatDat)[[2]]
    numGraphs <- sum(1:(len-1))
    parArg <- ceiling(sqrt(numGraphs))

    par(mfrow=c(parArg, parArg), pty="s")
    for(i in 1:(len-1)) {
       for(j in (i+1):len) {
          ran <- range(scatDat[,c(i,j)], na.rm=T)
          plot(scatDat[,i], scatDat[,j], xlab=paste("Plate", dimnames(scatDat)[[2]][i]),
               ylab=paste("Plate",dimnames(scatDat)[[2]][j]), xlim=ran,
               ylim=ran, main=paste("Plates:", dimnames(scatDat)[[2]][i], ",",
                dimnames(scatDat)[[2]][j]))
          lines(supsmu(scatDat[,i], scatDat[,j]), col=2)
       }
    }
}

#
# Plot histogram of p values.
#
#  pvals - matrix of data with each column a single HTS plate
#
pValueDistributionPlot <- function(pvals) {

    len <- dim(pvals)[[2]]
    numGraphs <- len
    parArg <- ceiling(sqrt(numGraphs))

    par(mfrow=c(parArg, parArg))

     for(i in 1:len) {
         hist(pvals[,i], breaks=20, main=paste("Plate", i), xlab="Bins", ylab="Count")
     }
}


#
# Plot autocorrelation for each plate
#
# mat - matrix of data with each column a single HTS plate
# numR - number of rows in HTS plate
# numC - number of columns in HTS plate
#
autocorrelationPlot <- function(acDat=mat, numR=numRows, numC=numCols)  {

   rowlag <- 3*numR
   collag <- 3*numC
   colInd <- rep(1:numC, numR)
   colSort <- order(colInd)

   # if acDat is a vector, turn it into a matrix with one column
   if(is.null(dim(acDat)[[2]])) {
       acDat <- matrix(acDat, ncol=1)
       colDat <-  matrix(acDat[colSort,1], ncol=1)
   }
   else {
       len <- dim(acDat)[[2]]
       colDat <- as.matrix(acDat[colSort,], ncol=len, byrow=F)
   }



   len <- dim(acDat)[[2]]
   parArg <- ceiling(sqrt(len))

   rowMat <- colMat <- NULL
   for(i in 1:len) {
      rowMat <- cbind(rowMat, acf(acDat[,i], lag.max=rowlag, plot=F, na.action=na.pass)$acf[-1])
      colMat <- cbind(colMat, acf(colDat[,i], lag.max=collag, plot=F, na.action=na.pass)$acf[-1])
   }

   par(mfrow=c(parArg, parArg))
   yran <- range(c(c(rowMat), c(colMat)))
   yran[1] <- ifelse(yran[1] > 0, 0, yran[1])
   for(i in 1:len) {
       plot(rowMat[,i], xlim=c(1, collag), ylim=yran, xlab="Lag", ylab="AutoCorrelation",
            main=paste("Plate", i), type="l")
       lines(1:collag, colMat[,i], col=2)
       abline(h=0, lty=2)
   }

}

# plot three dimensional plot of a single HTS plate
#
# visDat - matrix of data of a single HTS plate (one column)
# numR - number of rows in HTS plate
# numC - number of columns in HTS plate
#
visualPlot <- function(visDat, numR=numRows, numC=numCols) {

   score <- matrix(visDat, nrow=numR, ncol=numC)
   #score <- matrix(1:80, nrow=8, ncol=10)
   wireframe(score, drape=T)
}

# plot heatmap of a single HTS plate
#
# visDat - matrix of data of a single HTS plate (one column)
# numR - number of rows in HTS plate
# numC - number of columns in HTS plate
#
heatMapPlot <-  function(visDat, numR=numRows, numC=numCols) {

   rows <- rep(1:numR, numC)
   #rows <- (numR+1) - rows
   cols <- rep(1:numC, each=numR)
   revInd <- numR:1
    tmp <- matrix(visDat, nrow=numR, byrow=F)
    tmp1 <- tmp[revInd,]
   heatmap(matrix(tmp1, nrow=numR, byrow=F), xlab="Columns", ylab="Rows",Rowv=NA, Colv=NA,
           labRow=8:1)

    #levelplot(visDat~cols+rows, data=data.frame(visDat, rows, cols),
    #scales=list(y=list(at=c(8,6,4,2), labels=c(1,3,5,7))))
}


#
# plot inverse gamma fit
#
# dat - matrix of data, rows are elements and columns are replicates
#
IGFitPlot <- function(mat) {


    dat <- NULL
    for(i in unique(replicateIndex)) {
        dat <- rbind(dat, mat[,replicateIndex==i])
    }

   degFreedom1 <- dim(dat)[[2]] - 1
   vars <- rowVars(dat, na.rm=T)
   ab <- getab(vars, rep(degFreedom1, length(vars)))
   adj <- ab[1]*ab[2]
   adjVars <- vars*adj[1]
   scum <- my.cum(adjVars, q=1)
   probF <- pf(scum[,1], degFreedom1, 2*ab[1])
   lineWidth <- 2

    #kRes <- ks.gof(adjVars, distribution="f", df1=degFreedom1, df2=2*ab[1])$p.value

   plot(scum[,1], scum[,2], type="l", lwd=lineWidth, xlab="Var", ylab="cdf",
        main="RVM Inverse Gamma Fit")
   lines(scum[,1], scum[,2], col=2)
   lines(scum[,1], probF, col=5, lwd=lineWidth)
   legend("bottomright",  legend=c("Empirical", "Theo F"),lty=c(1,1), col=c(2,5))

}


#######################################################################################
###############################  Normalization utility functions ######################
#######################################################################################


#
# create a row neighbourhood list for matrix of data.  Assumes data will be
# sequentially ordered first by column,
# then row.  Each row defines row and column neighbours for that well.  Index values define which data
# points are neighbours in sequential data.  The neighboorhood includes the middle well
#
# nrows - number of rows in matrix
# ncols  - number of columns in matrix
# wind - window size, defines max distance of neighbours. wind of 2 will create a 5x5
#         neighbourhood along the wells row.
#
createSequentialNeighbourhoodMatrix <- function(ncols=10, nrows=8, wind=2) {

    #createSequentialNeighbourhoodMatrix(namey="E://projects//HTS//Inglese//spatial//data//IngleseConNeighboursRowWindow3.txt",
    #                                   ncols=40, nrows=32, wind=3)

    # total number of wells in matrix
    numWells <- nrows*ncols

    # create wasteful output matrix for max neighbourhood
    outMat <- matrix(NA, nrow=numWells, ncol=numWells)

    # reference matrix to identify sequential index of each well in matrix
    indMat <- matrix(1:numWells, nrow=nrows, ncol=ncols)

    # keeps track of where we are in sequential vector of matrix
    seqWellInd <- 0

    # keeps track of maximum number of buddies a well can have
    maxBuddies <- 0

   # iterate through each well location and find neighbours
   for(coll in 1:ncols) {
      for(roww in 1:nrows) {
         #cat("Finding buddies for", roww, " ", coll, "\n")

         seqWellInd <- seqWellInd + 1

         # vector to store neighbours
         buddies <- NULL

         # iterate through neighbourood by row
         for(j in (coll-wind):(coll+wind)) {

               # only keep valid neighbours (within range of matrix)
               if ((j >=1) & (j <= ncols)) {
                    #cat(roww, " ", j, "\n")
                    buddies <- c(buddies, indMat[roww,j])
               } # end if i

         } # end for j

          # assign neighbours to storage matrix
          outMat[seqWellInd, 1:length(buddies)] <- sort(buddies)

          # check if number of neighbours exceeds max neighbours yet encountered
          if(length(buddies) > maxBuddies) {maxBuddies <- length(buddies)}

       } # end for roww
   } # end for coll

    # remove excess columns
    outMat <- outMat[,1:maxBuddies]
    return(outMat)

} # end createSequentialNeighbourhoodMatrix



#
# create neighbourhood list for matrix of data.  Assumes data will be sequentially ordered
# first by column,
# then row.  Each row defines neighbours for that well.  Index values define which data
# points are neighbours in sequential data.  The neighboorhood includes the middle well
#
# nrows - number of rows in matrix
# ncols  - number of columns in matrix
# wind - window size, defines max distance of neighbours. wind of 2 will create a 5x5
#         neighbourhood.
#
createNeighbourhoodMatrix <- function(ncols=10, nrows=8, wind=2) {
                                       
    #createNeighbourhoodMatrix(namey="E://projects//HTS//Inglese//spatial//data//IngleseConNeighboursWindow4.txt",
    #                                   ncols=40, nrows=32, wind=4)

    # total number of wells in matrix
    numWells <- nrows*ncols

    # create wasteful output matrix for max neighbourhood
    outMat <- matrix(NA, nrow=numWells, ncol=numWells)

    # reference matrix to identify sequential index of each well in matrix
    indMat <- matrix(1:numWells, nrow=nrows, ncol=ncols)

    # keeps track of where we are in sequential vector of matrix
    seqWellInd <- 0

    # keeps track of maximum number of buddies a well can have
    maxBuddies <- 0

   # iterate through each well location and find neighbours
   for(coll in 1:ncols) {
      for(roww in 1:nrows) {
         #cat("Finding buddies for", roww, " ", coll, "\n")

         seqWellInd <- seqWellInd + 1

         # vector to store neighbours
         buddies <- NULL

         # iterate through neighbourood , first by row, then column
         for(i in (roww-wind):(roww+wind)) {
            for(j in (coll-wind):(coll+wind)) {

               # only keep valid neighbours (within range of matrix)
               if ((i >= 1) & (i <= nrows) & (j >=1) & (j <= ncols)) {
                    #cat(i, " ", j, "\n")
                    buddies <- c(buddies, indMat[i,j])
               } # end if i

             } # end for j
          } # end for i

          # assign neighbours to storage matrix
          outMat[seqWellInd, 1:length(buddies)] <- sort(buddies)

          # check if number of neighbours exceeds max neighbours yet encountered
          if(length(buddies) > maxBuddies) {maxBuddies <- length(buddies)}

       } # end for roww
   } # end for coll

    # remove excess columns
    outMat <- outMat[,1:maxBuddies]
    
    return(outMat)

}  # end createNeighbourhoodMatrix







#
# Get rid of the very large data points so graphs scale better
#
# quan - trimming quantile
#
getLen <- function(data, quan=.90) {
	return(trunc(quan*length(data)))

}





#
# Generate empirical probability density function for data.
#
# data - vector of data points
# q - trimming value.  remove 1-q points as outliers from greater tail.
#
my.cum <- function(data, q=.9) {

	#x <- seq(0.001,3.5,0.001)
	len <- getLen(data, quan=q)
	#len <- 1000
	maxi <- sort(data)[len]
	x <- seq(min(data), maxi, length=len)
	#x <- quantile(data[1:len], ppoints(100))
	p <- rep(0,len-1)
        lenny <- length(data)

	for(i in 1:len)
          p[i] <- (sum(data<x[i]))/lenny

       
	return(cbind(x, p))
}




#
# Trimmed polish.  Modified from median polish code
#
#
# Apply a trimmed median polish
#
# x: a numeric matrix.
# eps: real number greater than 0. A tolerance for convergence
# maxiter: the maximum number of iterations
# trace.iter: logical. Should progress in convergence be reported?
# na.rm: logical. Should missing values be removed?
# trim: trim value when applying trimmed mean
#
trimpolish <- function (x, eps = 0.01, maxiter = 10, trace.iter = FALSE, na.rm = TRUE, trim=0.5) 
{
    z <- as.matrix(x)
    nr <- nrow(z)
    nc <- ncol(z)
    t <- 0
    r <- numeric(nr)
    c <- numeric(nc)
    oldsum <- 0
    for (iter in 1:maxiter) {
        rdelta <- apply(z, 1, mean, na.rm = na.rm, trim=trim)
        z <- z - matrix(rdelta, nr = nr, nc = nc)
        r <- r + rdelta
        delta <- mean(c, na.rm = na.rm, trim=trim)
        c <- c - delta
        t <- t + delta
        cdelta <- apply(z, 2, mean, na.rm = na.rm, trim=trim)
        z <- z - matrix(cdelta, nr = nr, nc = nc, byrow = TRUE)
        c <- c + cdelta
        delta <- mean(r, na.rm = na.rm, trim=trim)
        r <- r - delta
        t <- t + delta
        newsum <- sum(abs(z), na.rm = na.rm)
        converged <- newsum == 0 || abs(newsum - oldsum) < eps * 
            newsum
        if (converged) 
            break
        oldsum <- newsum
        if (trace.iter) 
            cat(iter, ":", newsum, "\n")
    }
    if (converged) {
        if (trace.iter) 
            cat("Final:", newsum, "\n")
    }
    else warning(gettextf("medpolish() did not converge in %d iterations", 
        maxiter), domain = NA)
    names(r) <- rownames(z)
    names(c) <- colnames(z)
    ans <- list(overall = t, row = r, col = c, residuals = z,
        name = deparse(substitute(x)))
    class(ans) <- "medpolish"
    ans
}

#
# R implementation of Simon RVM algorithm.  The only modification (I believe)
# is in the getab function where flik is passed in differently.
#
# Wright and Simon. A random variance model for detection of differential gene expression
#                   in small microarrray experiments.
#


# do a one sided one sample test (greater/less than 0)
"RVT1Sided"<-
function(data, side="less")  #one sample T-test, detects average gene expression different from 0
                # data has genes as rows arrays as collumns
				  #
				  # function output is data.frame with a row for each gene
				  # mn=mean value for gene
				  # vr=unadjusted variance for gene
				  # n is number of non-missing samples for gene
				  # t is unadjusted t-statistic
				  # tp is unadjusted t p-value
				  # v is RVM model statistic
				  # vp is RVM model p-value
				  # vvr is adjusted variance
{	vr<-rowVars(data,na.rm=T)
	n<-rowSums(!is.na(data))
	mn<-rowMeans(data,na.rm=T)
	
	# NOTE: Originally this function used n-2 for df.  I believe this is a copy and
	# paste error from RVT2.  The correct degree of freedom is n-1 for a one-sample test
	a<-getab(vr,n-1)
	b<-a[2]
	a<-a[1]
	cat("alpha=",a,"beta=",b,"\n")
	t<-mn/sqrt(vr/n)
	tp<-1-pt(t,df=n-1)
	vvr<-n*(1+2*a/(n-1))/(vr+2/((n-1)*b))
	v<-mn*sqrt(vvr)
	degf = n-1+2*a

        if(side=="less") {
           tp<-pt(t,df=n-1)
           vp <- pt(v,df=degf)

        }
        else if(side=="greater") {
           tp<-1-pt(t,df=n-1)
           vp <- 1-pt(v,df=degf)
        }
        
	return(data.frame(n, mn, vr, t, tp, vvr, v, vp, a, b))
}




"RVT1"<-
function(data)  #one sample T-test, detects average gene expression different from 0
                # data has genes as rows arrays as collumns
				  #
				  # function output is data.frame with a row for each gene
				  # mn=mean value for gene
				  # vr=unadjusted variance for gene
				  # n is number of non-missing samples for gene
				  # t is unadjusted t-statistic
				  # tp is unadjusted t p-value
				  # v is RVM model statistic
				  # vp is RVM model p-value
				  # vvr is adjusted variance
{	vr<-rowVars(data,na.rm=T)
	n<-rowSums(!is.na(data))
	mn<-rowMeans(data,na.rm=T)
	
	# NOTE: Originally this function used n-2 for df.  I believe this is a copy and
	# paste error from RVT2.  The correct degree of freedom is n-1 for a one-sample test
	a<-getab(vr,n-1)
	b<-a[2]
	a<-a[1]
	cat("alpha=",a,"beta=",b,"\n")
	t<-mn/sqrt(vr/n)
	tp<-2*(1-pt(abs(t),df=n-1))
	vvr<-n*(1+2*a/(n-1))/(vr+2/((n-1)*b))
	v<-mn*sqrt(vvr)
	degf = n-1+2*a
	vp <- 2*(1-pt(abs(v),df=degf))
	
	return(data.frame(n, mn, vr, t, tp, vvr, v, vp, a, b))
}


#
# Returns vp - posterior p value
#
"RVT2" <-
function(data, labels, ab=NULL)
{	data<-data[,!is.na(labels)]
	labels<-labels[!is.na(labels)]
	a<-unique(labels)
	a<-a[order(a)]
	labels<-as.integer(labels==a[1])
	dat1<-data[,labels==1]
	dat2<-data[,labels==0]

	vr1<-rowVars(dat1,na.rm=T)
	vr2<-rowVars(dat2,na.rm=T)
	n1<-rowSums(!is.na(dat1))
	n2<-rowSums(!is.na(dat2))
	mn1<-rowMeans(dat1,na.rm=T)
	mn2<-rowMeans(dat2,na.rm=T)
	mn<-mn2-mn1
	n<-n1+n2
	vr<-(n1-1)*vr1+(n2-1)*vr2
	vr<-vr/(n-2)
	vr[vr==0]<-0.001

	if (is.null(ab)) {
  	    a<-getab(vr,n-2)
	    b<-a[2]
	    a<-a[1]
	}
	else {
	    a <- ab[1]
	    b <- ab[2]
	}
	#cat("alpha=",a,"beta=",b,"\n")
	newn<-1/(1/n1+1/n2)
	t<-mn/sqrt(vr/newn)
	n[n<3]<-3
	tp<-2*(pt(-abs(t),df=n-2))
	vvr<-newn*(1+2*a/(n-2))/(vr+2/((n-2)*b)) # inverse of standard error

        # TEST
        #vr1 <- 1/((1+2*a/(n-2))/(vr+2/((n-2)*b)))
        #t1 <- mn/sqrt(vr1/newn)
        #se1 <- sqrt(vr/newn)
        #se2 <- sqrt(vr1/newn)
	v<-mn*sqrt(vvr)
	degf = n-2+2*a
	vp<-2*(1-pt(abs(v),df=degf))

	if (!is.null(ab)) {
	    # only interested in values that were influenced by the given a and b
	    return(data.frame(sqrt(1/vvr), vp, degf))
	}

	vvr=newn/vvr
	dataTable<-data.frame(n1,n2,mn1,mn2,vr1,vr2,mn,vr,n,t,tp,v,vp,vvr,degf)
	return(list(dataTable=dataTable, ab=c(a,b)))
}



# Random variance F-test, labels is set of class lables,
# data has genes as rows arrays as collumns 
#
# Fval unadjusted F statistic
# VFval  RVM F statistic
# Fp unadjusted F p-value
# VFval  RVM F p-value
# function output is data.frame with a row for each gene
# ni, n is number of non-missing samples in subclasses and total data
# mni = mean value for gene in each subclass
# vri unadjusted variance for gene in each subclass
RVF <- function(data,lables) {
  
        lab <- unique(lables)
	gnum<-dim(data)[1]
	clnum<-length(lab)
	mni<-ni<-vri<-matrix(NA,gnum,clnum)
	for(i in 1:clnum) {
          mni[,i]<-rowMeans(data[,lables==lab[i]],na.rm=T)
          vri[,i]<-rowVars(data[,lables==lab[i]],na.rm=T)
          ni[,i]<-rowSums(!is.na(data[,lables==lab[i]]),na.rm=T)
	}		
	n<-rowSums(ni)
	vr<-rowVars(data,na.rm=T)
	residssq<-rowSums((ni-1)*vri)
	totssq<-vr*(n-1)
	num<-(totssq-residssq)/(clnum-1)
	resdf<-n-clnum
	denom<-(residssq)/(resdf)
	Fval<-num/denom
	a<-getab(denom,resdf)
	#cat(a,resdf[1],"\n")
	b<-a[2]
	a<-a[1]
	denom2<-((resdf)*denom+2/b)/(resdf+2*a)
        resdf2<-resdf+2*a
	VFval<-num/denom2
	Fp<-1-pf(Fval,clnum-1,resdf)
	VFp<-1-pf(VFval,clnum-1,resdf2)

	dataTable <-data.frame(1:gnum,Fval,VFval,Fp,VFp,denom2,mni,vri,ni,n)
        return(list(dataTable=dataTable, ab=c(a,b)))
}



# Calculate parameters of prior distribution
"getab"<-   #returns estiamtes of a and b
function(sig,n)
{	set<-(!is.na(sig)&n>0&sig>0)
	sig<-sig[set]
	n<-n[set]
	set<-n>4
	#cat("here")
	if (sum(set)>0)
	{  	m1<-(n[set]-2)/((n[set])*sig[set])
		m2<-(n[set]-2)*(n[set]-4)/((n[set])*sig[set])^2
		m1<-mean(m1,na.rm=T)
		m2<-mean(m2,na.rm=T)
		b<-m2/m1-m1
		a<-m1^2/(m2-m1^2)
	}
	else{ a<-b<-1}
	#cat(a,b,"initial ab\n")
	strt<-c(a,b)
	#g<-function(p,y) flik(p,y)
	#g$y<-c(sig,n)
###PATCH
	g<-function(p,yunq) flik(p,yunq)
###
	####a<-nlm(g,strt,initstep=.1,xc.tol=0.01)
	a<-nlm(g,strt, yunq=c(sig,n))
	#a$x<-abs(a$x)
	#cat(a$x, "final ab\n")
	a$estimate<-abs(a$estimate)
}




"flik"<-
function(p,y)
# log liklihood for a*b*x from an F distribution with m and 2*a degrees of freedom
# y is a vector containing data and the m values, p contains a and b.
{	x<-y[1:(length(y)/2)]
	m<-y[(length(y)/2+1):length(y)]
	p<-abs(p)
	a<-p[1]
	b<-p[2]
	x<-x*(a*b)
	n<-2*a
	out<-base::log(df(x,m,n))+base::log(a*b)
	sum(-out)
}


rowVars<-function(x, ...) { apply(x, 1, var, ...) }


