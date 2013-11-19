# some Encode RNA-Seq BigWigs each 160 Mb
ftpPath <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/"
ftpFiles <- c("wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig",
              "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep2.bigWig")

# download if they don't exist
for (f in ftpFiles[!file.exists(ftpFiles)]) download.file(paste0(ftpPath,f),destfile=f)
fls <- ftpFiles

library(rtracklayer)
source("AllClasses.R")
source("AllGenerics.R")
source("methods-BigWigViews.R")

# just make a SimpleList of coverage
gr <- GRanges(c("chr1","chr1","chr2"),IRanges(c(43e6,147e6,74e6)+1,width=2e5))
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

# get integer coverage over the GRanges specified by BigWigViews
intCoverageMatrix <- function(bmv) {
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
z <- intCoverageMatrix(bmv)
print(object.size(z),unit="Mb")

# plot coverage across replicates
idx <- rowSums(z) > 0
plot(z[idx,1],z[idx,2],cex=.5)