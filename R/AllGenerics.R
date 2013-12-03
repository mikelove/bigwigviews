#' Constructor for BigWigViews objects
#'
#' @param bigWigPaths character vector
#' @param bigWigSamples DataFrame
#' @param bigWigRanges GRanges
#' @param bigWigExperiment a list
#'
#' @export
setGeneric("BigWigViews",
           function(bigWigPaths=character(0),
                    bigWigSamples=DataFrame(row.names=make.unique(basename(bigWigPaths))),
                    bigWigRanges,
                    bigWigExperiment=list(), ...)
             standardGeneric("BigWigViews"),
           signature="bigWigRanges")
