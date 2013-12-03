setGeneric(".validity", function(object) standardGeneric(".validity"))

#' BigWigViews class
#'
#' A lightweight pointer to a series of BigWigs and a set
#' of GRanges over which to view the data.
setClass("BigWigViews",
         representation=representation(
           bigWigPaths="character",
           bigWigSamples="DataFrame",
           bigWigRanges="GRanges",
           bigWigExperiment="list"),
         validity=.validity)
