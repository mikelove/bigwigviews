setMethod(BigWigViews, c(bigWigRanges="GRanges"), 
          function(bigWigPaths=character(0),
                   bigWigSamples=DataFrame(row.names=make.unique(basename(bigWigPaths))),
                   bigWigRanges,
                   bigWigExperiment=list(), ...)
          {
            new("BigWigViews", ..., bigWigPaths=bigWigPaths,
                bigWigSamples=bigWigSamples, bigWigRanges=bigWigRanges,
                bigWigExperiment=bigWigExperiment)
          })

setMethod(.validity, "BigWigViews", function(object) {
  msg <- NULL
  if (length(bigWigPaths(object)) != nrow(bigWigSamples(object)))
    msg <- c(msg,
             "length(bigWigPaths(object)) != nrow(bigWigSamples(object))")
  if (is.null(msg)) TRUE else msg
})

bigWigPaths <-
  function(x) setNames(slot(x, "bigWigPaths"), names(x))

`bigWigDirname<-` <-
  function(x, ..., value)
  {
    initialize(x,
               bigWigPaths=file.path(value, basename(bigWigPaths(x))))
  }

bigWigSamples <-
  function(x) slot(x, "bigWigSamples")

`bigWigSamples<-` <-
  function(x, value) initialize(x, bigWigSamples=value)

bigWigRanges <-
  function(x) slot(x, "bigWigRanges")

`bigWigRanges<-` <-
  function(x, value) initialize(x, bigWigRanges=value)

bigWigExperiment <-
  function(x) slot(x, "bigWigExperiment")

setMethod(dim, "BigWigViews", function(x) {
  c(length(bigWigRanges(x)), length(bigWigPaths(x)))
})

setMethod(names, "BigWigViews", function(x) {
  rownames(bigWigSamples(x))
})

.BigWigViews_delegate <-
  function(what, bigWigViews, fun, ...)
  {
    idx <- structure(seq_len(ncol(bigWigViews)), names=names(bigWigViews))
    result <- lapply(idx, fun, bigWigViews, ...)
    if (length(result) != ncol(bigWigViews)) {
      stop(sprintf("'%s' failed on '%s'", what,
                   paste(setdiff(names(bigWigViews), names(result)),
                         collapse="' '")))
    }
    names(result) <- names(bigWigViews)
    do.call(new, list("SimpleList", listData=result,
                      elementMetadata=bigWigSamples(bigWigViews)))
  }

# TODO check on this split() copied over from the BamViews code
.BigWigViews_which <- function(file)
{
  grange <- bigWigRanges(file)
  which <- split(ranges(grange), seqnames(grange))
  which
}

setMethod(coverage, "BigWigViews",
          function(x, ...)
          {
            which <- .BigWigViews_which(x)
            fun <- function(i, BigWigViews, ..., verbose) {
              bigWigGRanges <- import(BigWigFile(bigWigPaths(BigWigViews)[i]), which=which)
              cvr <- coverage(bigWigGRanges, weight=mcols(bigWigGRanges)$score, ...)
            }
            .BigWigViews_delegate("coverage", x, fun, ...)
          })