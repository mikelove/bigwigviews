BigWigViews
===========

R version of desired BigWigViews functionality

The BigWig paper is here:

http://bioinformatics.oxfordjournals.org/content/26/17/2204.long

I just grabbed the necessary classes, generics and methods from Rsamtools::BamViews, removing the bamIndices slot, because sample and index are combined in BigWig files.

Two novel bits of code:

- the coverage method for BigWigViews, added to methods-BigWigViews.R
- which *apply function to use in .BigWigViews_delegate(). Currently I use lapply:

```R
    result <- lapply(idx, fun, bigWigViews, ...)
```

An example function for getting an integer matrix is intCoverageMatrix() in example_usage.R

# Design questions

- The current implementation of BamViews parallelizes over samples using .BamViews_delegate. But we might think about parallelization over ranges instead? I have used lapply() in .BigWigViews_delegate() rather than ShortRead::srapply. Or do we want to parallelize over both?

- BigWigFile instead of character vector for specifying bigWigPaths. BigWigFile has methods seqinfo() and summary() which make use of the information contained in the BigWig. E.g. then you can warn about conflicts between bigWigRanges and the BigWig files.

- What about the BigWigSelection argument to import.bw

