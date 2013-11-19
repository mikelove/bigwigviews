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
gr <- GRanges(c("chr1","chr1","chr2"),IRanges(c(43e6,147e6,74e6)+1,width=1e6))
bwv <- BigWigViews(bigWigPaths=fls, bigWigRanges=gr)
z <- coverage(bwv)
print(object.size(z),unit="Mb")

# if the bigWigRange is over a single chromosome, 
# this should get integer coverage over these ranges
intCoverageMatrix <- function(bmv) {
  cvrOverBigWigs <- coverage(bwv)
  cvrList <- lapply(cvrOverBigWigs, function(cvr) {
    charRangesNames <- as.character(seqnames(bigWigRanges(bwv)))
    cvrOverSeqs <- lapply(charRangesNames, function(seqname) {
      bwr <- bigWigRanges(bwv)
      viewApply(Views(cvr[[seqname]], ranges(bwr[seqnames(bwr) == seqname])), as.integer)
    })
    do.call(c, cvrOverSeqs)
  })
  do.call(cbind, cvrList)
}

# this gives the matrix of coverage
z <- intCoverageMatrix(bmv)
print(object.size(z),unit="Mb")

# plot coverage across replicates
idx <- rowSums(z) > 0
plot(z[idx,1],z[idx,2],cex=.1,log="xy",col=rgb(0,0,0,.2))