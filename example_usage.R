# download some Encode RNA-Seq BigWigs 4-16 Mb
ftpPath <- "http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/rna_signal/jan2011/hub/"
ftpFiles <- c("wgEncodeCshlShortRnaSeqK562CellShortMinusRaw.bigWig",
              "wgEncodeCshlShortRnaSeqK562CellShortPlusRaw.bigWig",
              "wgEncodeCshlShortRnaSeqK562CellShorttotalTapMinusRawRep1.bigWig",
              "wgEncodeCshlShortRnaSeqK562CellShorttotalTapMinusRawRep2.bigWig",
              "wgEncodeCshlShortRnaSeqK562CellShorttotalTapPlusRawRep1.bigWig", 
              "wgEncodeCshlShortRnaSeqK562CellShorttotalTapPlusRawRep2.bigWig")
for (f in ftpFiles[!file.exists(ftpFiles)]) {
  download.file(paste0(ftpPath,f),destfile=f)
}
fls <- ftpFiles

library(rtracklayer)
gr <- GRanges("chr1",IRanges(15e6+1,width=5e5))
bwv <- BigWigViews(bigWigPaths=fls, bigWigRanges=gr)
z <- coverage(bwv)

# if the bigWigRange is over a single chromosome, 
# this will get integer coverage over these ranges
intCovMatrix <- function(bmv) {
  cvrOverBigWigs <- coverage(bwv)
  cvrList <- lapply(cvrOverBigWigs, function(cvr) {
    as.integer(Views(cvr[[as.character(seqnames(bigWigRanges(bwv))[1])]], 
                     ranges(bigWigRanges(bwv)))[[1]])
  })
  do.call(cbind, cvrList)
}

z <- intCovMatrix(bmv)