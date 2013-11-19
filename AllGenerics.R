setGeneric("BigWigViews",
           function(bigWigPaths=character(0),
                    bigWigSamples=DataFrame(row.names=make.unique(basename(bigWigPaths))),
                    bigWigRanges,
                    bigWigExperiment=list(), ...)
             standardGeneric("BigWigViews"),
           signature="bigWigRanges")