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

setMethod(BigWigViews, c(bigWigRanges="RangedData"), 
          function(bigWigPaths=character(0),
                   bigWigSamples=DataFrame(row.names=make.unique(basename(bigWigPaths))),
                   bigWigRanges,
                   bigWigExperiment=list(), ...)
          {
            bigWigRanges <- as(bigWigRanges, "GRanges")
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

#' @export
bigWigPaths <-
  function(x) setNames(slot(x, "bigWigPaths"), names(x))

#' @export
`bigWigDirname<-` <-
  function(x, ..., value)
  {
    initialize(x,
               bigWigPaths=file.path(value, basename(bigWigPaths(x))))
  }

#' @export
bigWigSamples <-
  function(x) slot(x, "bigWigSamples")

#' @export
`bigWigSamples<-` <-
  function(x, value) initialize(x, bigWigSamples=value)

#' @export
bigWigRanges <-
  function(x) slot(x, "bigWigRanges")

#' @export
`bigWigRanges<-` <-
  function(x, value) initialize(x, bigWigRanges=value)

#' @export
bigWigExperiment <-
  function(x) slot(x, "bigWigExperiment")

setMethod(dim, "BigWigViews", function(x) {
  c(length(bigWigRanges(x)), length(bigWigPaths(x)))
})

setMethod(names, "BigWigViews", function(x) {
  rownames(bigWigSamples(x))
})

setReplaceMethod("names", "BigWigViews", function(x, value) {
  rownames(bigWigSamples(x)) <- value
  x
})

setMethod(dimnames, "BigWigViews", function(x) {
  list(names(bigWigRanges(x)),
       rownames(bigWigSamples(x)))
})

setReplaceMethod("dimnames", "BigWigViews", function(x, value) {
  names(bigWigRanges(x)) <- value[[1]]
  rownames(bigWigSamples(x)) <- value[[2]]
  x
})

setMethod("[", c("BigWigViews", "ANY", "missing"),
          function(x, i, j, ..., drop=TRUE)
          {
            initialize(x, bigWigRanges=bigWigRanges(x)[i,])
          })

setMethod("[", c("BigWigViews", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE)
          {
            if (is.character(j))
              j <- match(j, colnames(x))
            if (any(is.na(j)))
              stop("subscript 'j' out of bounds")
            initialize(x, bigWigPaths=bigWigPaths(x)[j],
                       bigWigSamples=bigWigSamples(x)[j,,drop=FALSE])
          })

setMethod("[", c("BigWigViews", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE)
          {
            if (is.character(i))
              j <- match(i, rownames(x))
            if (is.character(j))
              j <- match(j, colnames(x))
            if (any(is.na(i)))
              stop("subscript 'i' out of bounds")
            if (any(is.na(j)))
              stop("subscript 'j' out of bounds")
            initialize(x, bigWigRanges=bigWigRanges(x)[i,],
                       bigWigPaths=bigWigPaths(x)[j],
                       bigWigSamples=bigWigSamples(x)[j,,drop=FALSE])
          })


# NOTE: here I am just using lapply() instead of srapply / fapply as in BigWigViews
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

setMethod(coverage, "BigWigViews",
          function(x, ...)
          {
            fun <- function(i, bigWigViews, ..., verbose) {
              import(BigWigFile(bigWigPaths(bigWigViews)[i]), 
                     which=bigWigRanges(bigWigViews), 
                     asRle=TRUE, ...)
            }
            .BigWigViews_delegate("coverage", x, fun, ...)
          })

setMethod(show, "BigWigViews", function(object) {
  cat(class(object), "dim:",
      paste(dim(object), c("ranges", "samples"), collapse=" x "),
      "\n")
  cat("names:", BiocGenerics:::selectSome(names(object)), "\n")
  cat("detail: use bigWigPaths(), bigWigSamples(), bigWigRanges(), ...",
      "\n")
})
