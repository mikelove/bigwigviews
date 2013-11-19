# download some Encode RNA-Seq BigWigs each 160 Mb
ftpPath <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/"
ftpFiles <- c("wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig",
              "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep2.bigWig")

for (f in ftpFiles[!file.exists(ftpFiles)]) {
  download.file(paste0(ftpPath,f),destfile=f)
}
fls <- ftpFiles

# just make a SimpleList of coverage
library(rtracklayer)
gr <- GRanges("chr1",IRanges(15e6+1,width=5e5))
bwv <- BigWigViews(bigWigPaths=fls, bigWigRanges=gr)
z <- coverage(bwv)

# if the bigWigRange is over a single chromosome, 
# this should get integer coverage over these ranges
intCovMatrix <- function(bmv) {
  cvrOverBigWigs <- coverage(bwv)
  cvrList <- lapply(cvrOverBigWigs, function(cvr) {
    as.integer(Views(cvr[[as.character(seqnames(bigWigRanges(bwv))[1])]], 
                     ranges(bigWigRanges(bwv)))[[1]])
  })
  do.call(cbind, cvrList)
}

# this gives the matrix of coverage
z <- intCovMatrix(bmv)