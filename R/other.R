# get Rle coverage for a single range:
coverageSingleRange <- function(bigWigViews, i) {
  idx <- structure(seq_len(ncol(bigWigViews)), names=names(bigWigViews))
  bwr <- bigWigRanges(bigWigViews[i,])
  lapply(idx, function(j) {
    cvr <- coverage(bigWigViews[i,])[[j]][[as.character(seqnames(bwr))]]
    if (end(ranges(bwr)) > length(cvr)) {
      stop("ranges in BigWigViews possibly extend beyond the end defined in the BigWig")
    }
    Views(cvr, ranges(bwr))[[1]]
  })
}

# some ridiculous functions
rowSumsListRles <- function(l) {
  Reduce("+",l)
}
rowMeansListRles <- function(l) {
  Reduce("+",l) / length(l)
}

# hackety hackety t-test
rleTTest <- function(rles, idx1, idx2, s0=1) {
  n1 <- length(idx1)
  n2 <- length(idx2)
  mean1 <- rowMeansListRles(rles[idx1])
  mean2 <- rowMeansListRles(rles[idx2])
  sse1 <- rowSumsListRles(lapply(rles[idx1], function(x) (x - mean1)^2))  
  sse2 <- rowSumsListRles(lapply(rles[idx2], function(x) (x - mean2)^2))
  s <- sqrt( ( sse1 + sse2 ) / (n1 + n2 - 2) )
  t <- (mean1 - mean2) / ((s + s0) * sqrt( 1/n1 + 1/n2 ))
  t
}

# get integer coverage over the GRanges specified by BigWigViews
intCoverageMatrix <- function(bwv) {
  # get the SimpleList of RleLists
  cvrOverBigWigs <- coverage(bwv)
  bwr <- bigWigRanges(bwv)
  charRangesNames <- as.character(seqnames(bwr))
  rangesList <- split(ranges(bwr),seqnames(bwr))
  cvrList <- lapply(cvrOverBigWigs, function(cvr) {
    listOfLists <- viewApply(RleViewsList(rleList=cvr[names(cvr) %in% charRangesNames],
                                          rangesList=rangesList), as.integer, simplify=FALSE)
    # turn the list of lists into a long vector
    do.call(c, lapply(listOfLists, function(x) do.call(c, as.list(x))))
  })
  # bind the vectors from each sample
  do.call(cbind, cvrList)
}
