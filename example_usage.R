# some Encode RNA-Seq BigWigs each 160 Mb
ftpPath <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/"
ftpFiles <- c("wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig",
              "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep2.bigWig",
              "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig",
              "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig")

# download if they don't exist
for (f in ftpFiles[!file.exists(ftpFiles)]) download.file(paste0(ftpPath,f),destfile=f)
fls <- ftpFiles

library(rtracklayer)
source("AllClasses.R")
source("AllGenerics.R")
source("methods-BigWigViews.R")

# construct a BigWigViews instance
gr <- tileGenome(c("chr1"=249e6),cut.last.tile.in.chrom=TRUE,tilewidth=5e6)
bwv <- BigWigViews(bigWigPaths=fls, bigWigRanges=gr)

# show
bwv

# the carry-over methods from BamViews
names(bwv)
dim(bwv)
dimnames(bwv)
bwv[1:2,2]

# try the coverage method
# this gives a SimpleList for each sample
# with elements: RleList of coverage over each seq
z <- coverage(bwv)
print(object.size(z),unit="Mb")

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

z <- coverageSingleRange(bwv,1)
print(object.size(z),unit="Mb")

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

# stream along the genomic ranges and calculate t tests
# this should also go and grab the scaling factor from the bigWigSamples DataFrame
system.time({ts <- lapply(seq_len(nrow(bwv)), function(i) {
  cvr <- coverageSingleRange(bwv,i)
  t <- rleTTest(cvr, 1:2, 3:4)
  t
})})


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

# this gives the matrix of coverage
z <- intCoverageMatrix(bwv)
print(object.size(z),unit="Mb")

# plot coverage across replicates
idx <- rowSums(z) > 0
plot(z[idx,1],z[idx,2],cex=.5)