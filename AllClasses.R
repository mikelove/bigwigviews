setGeneric(".validity", function(object) standardGeneric(".validity"))

setClass("BigWigViews",
         representation=representation(
           bigWigPaths="character",
           bigWigSamples="DataFrame",
           bigWigRanges="GRanges",
           bigWigExperiment="list"),
         validity=.validity)