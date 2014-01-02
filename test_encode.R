
# some Encode RNA-Seq BigWigs each 160 Mb
ftpPath <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/"
ftpFiles <- c("wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig",
              "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep2.bigWig",
              "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig",
              "wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep2.bigWig")

localPath <- system.file("testData", package="bigwigviews")
library(rtracklayer)
fls <- paste0(localPath, "/", ftpFiles)

# construct a BigWigViews instance
gr <- GRanges("chr1",IRanges(11:20 * 1e6 + 1,width=1e4))
bwv <- BigWigFileViews(filePaths=fls, fileRanges=gr)

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

z <- coverageSingleRange(bwv,1)
print(object.size(z),unit="Mb")

# stream along the genomic ranges and calculate t tests
# this should also go and grab the scaling factor from the bigWigSamples DataFrame
t <- system.time({ts <- lapply(seq_len(nrow(bwv)), function(i) {
  cvr <- coverageSingleRange(bwv,i)
  t <- rleTTest(cvr, 1:2, 3:4)
  t
})})

print(t)

# subset as we are going making a dense matrix
gr <- tileGenome(c("chr1"=249e6),cut.last.tile.in.chrom=TRUE,tilewidth=1e5)
bwv <- BigWigViews(bigWigPaths=fls, bigWigRanges=gr)
bwv <- bwv[101:110,]

# this gives the matrix of coverage
z <- intCoverageMatrix(bwv)
print(object.size(z),unit="Mb")

# plot coverage across replicates
# idx <- rowSums(z) > 0
# plot(z[idx,1],z[idx,2],cex=.5)

