#############################
## Identify Range Overlaps ##
#############################
## Author: Thomas Girke
## Last update: 8-Feb-11
## Details on usage and use cases are available here:
## http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-Analysis-Routines-with-IRanges-Geno

## Utility: identify overlaps in range data sets, such as annotation or alignment positions defined
## by two IRanges/GRanges objects.  

## Overlap types
## olup: startup & endin
## Q --------------
## S       -------------

## oldown: startin & enddown
## Q       -------------
## S --------------

## inside: startin & endin 
## Q     -----
## S --------------

## contained: startup & enddown
## Q --------------
## S     -----

###########################################################
## (A) olRanges Function for IRanges and GRanges Objects ##
###########################################################
olRanges <- function(query, subject, output="gr", ...) {
  require(GenomicRanges); require(IRanges)
  
  ## Input check
  if(!((class(query)=="GRanges" & class(subject)=="GRanges") | (class(query)=="IRanges" & class(subject)=="IRanges"))) {
    stop("Query and subject need to be of same class, either GRanges or IRanges!")
  }
  
  ## Find overlapping ranges
  if(class(query)=="GRanges") {
    seqlengths(query) <- rep(NA, length(seqlengths(query)))
    seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
  }
  olindex <- as.matrix(findOverlaps(query, subject, ...))
  query <- query[olindex[,1]]
  subject <- subject[olindex[,2]]
  olma <- cbind(Qstart=start(query), Qend=end(query), Sstart=start(subject), Send=end(subject))
  
  ## Pre-queries for overlaps
  startup <- olma[,"Sstart"] < olma[,"Qstart"]
  enddown <- olma[,"Send"] > olma[,"Qend"]
  startin <- olma[,"Sstart"] >= olma[,"Qstart"] & olma[,"Sstart"] <= olma[,"Qend"]
  endin <- olma[,"Send"] >= olma[,"Qstart"] & olma[,"Send"] <=  olma[,"Qend"]
  
  ## Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin 
  contained <- startup & enddown
  
  ## Overlap types in one vector
  OLtype <- rep("", length(olma[,"Qstart"]))
  OLtype[olup] <- "olup"
  OLtype[oldown] <- "oldown"
  OLtype[inside] <- "inside" 
  OLtype[contained] <- "contained"
  
  ## Overlap positions
  OLstart <- rep(0, length(olma[,"Qstart"]))
  OLend <- rep(0, length(olma[,"Qstart"]))
  OLstart[olup] <- olma[,"Qstart"][olup]
  OLend[olup] <- olma[,"Send"][olup]
  OLstart[oldown] <- olma[,"Sstart"][oldown]
  OLend[oldown] <- olma[,"Qend"][oldown]
  OLstart[inside] <- olma[,"Sstart"][inside]
  OLend[inside] <- olma[,"Send"][inside]
  OLstart[contained] <- olma[,"Qstart"][contained]
  OLend[contained] <- olma[,"Qend"][contained]
  
  ## Absolute and relative length of overlaps
  OLlength <- (OLend - OLstart) + 1
  OLpercQ <- OLlength/width(query)*100
  OLpercS <- OLlength/width(subject)*100
  
  ## Output type
  oldf <- data.frame(Qindex=olindex[,1], Sindex=olindex[,2], olma, OLstart, OLend, OLlength, OLpercQ, OLpercS, OLtype)
  if(class(query) == "GRanges") {
    oldf <- cbind(space=as.character(seqnames(query)), oldf)
  }
  if(output=="df") {
    return(oldf)
  }
  if(output=="gr") {
    if(class(query)=="GRanges") {
      elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf)
    }
    if(class(query)=="IRanges") {
      query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start=oldf[,"Qstart"], end=oldf[,"Qend"]), strand = Rle(strand(rep("+", length(query)))), oldf)  
    }
    return(query)
  }
}

## Run olRanges function
## Sample Data Sets
# grq <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), 
#                ranges = IRanges(seq(1, 100, by=10), end = seq(30, 120, by=10)), 
#                strand = Rle(strand(c("-", "+", "-")), c(1, 7, 2)))
# grs <- shift(grq[c(2,5,6)], 5)
# olRanges(query=grq, subject=grs, output="df") 
# olRanges(query=grq, subject=grs, output="gr") 


################################################
## (B) Older rangeOL Function for Data Frames ##
################################################
rangeOL <- function(featureDF=featureDF, label="ID1", start=25, end=35) {
  ## Pre-queries for overlaps
  qstartup <- start < featureDF[,1]
  qenddown <- end > featureDF[,2]
  qstartin <- start >= featureDF[,1] & start <= featureDF[,2]
  qendin <- end >= featureDF[,1] & end <= featureDF[,2]
  
  ## Overlap types
  qolup <- qstartup & qendin
  qoldown <- qstartin & qenddown
  qinside <- qstartin & qendin 
  qcontains <- qstartup & qenddown
  
  ## Overlap types in one vector
  OLtype <- rep("", length(featureDF[,1]))
  OLtype[qolup] <- "olup"
  OLtype[qoldown] <- "oldown"
  OLtype[qinside] <- "inside" 
  OLtype[qcontains] <- "contained"
  
  ## Overlap Positions
  OLstart <- rep("", length(featureDF[,1]))
  OLend <- rep("", length(featureDF[,1]))
  OLstart[qolup] <- featureDF[,1][qolup]
  OLend[qolup] <- end
  OLstart[qoldown] <- start
  OLend[qoldown] <- featureDF[,2][qoldown]
  OLstart[qinside] <- start
  OLend[qinside] <- end
  OLstart[qcontains] <- featureDF[,1][qcontains]
  OLend[qcontains] <- featureDF[,2][qcontains]
  
  ## Length of overlaps
  OLlength <- (as.numeric(OLend) - as.numeric(OLstart)) + 1 
  
  ## Feature label
  OLquery <- rep("", length(featureDF[,1]))
  OLquery[nchar(OLtype)>0] <- paste(label, " (", start, ":", end, ")", sep="")
  
  ## Output
  featureDF <- cbind(featureDF, OLquery, OLtype, OLstart, OLend, OLlength)
  return(featureDF[nchar(OLtype)>0,])
}

## Run rangeOL function
# featureDF <- data.frame(Start=c(1, 20, 90, 150, 180), End=c(40, 30, 105, 160, 200), Feature="mRNA")
# rangeOL(featureDF=featureDF, label="ID1", start=25, end=35)
# lapply(seq(along=featureDF[,1]), function(x) rangeOL(featureDF=featureDF, label=featureDF[x,3], start=as.numeric(featureDF[x,1]), end=as.numeric(featureDF[x,2])))

