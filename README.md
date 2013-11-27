BigWigViews
===========

R version of desired BigWigViews functionality

Note: being developed with Bioconductor version 2.14

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

- BamViews offers parallelization over samples using srapply() in .BamViews_delegate(). We are focusing with BigWigViews on chunking over ranges instead. But we have to worry about I/O. Are we assuming users put bigwigs in scratch?

- coverage on a BigWigViews with a single range currently gives the RleList of all seqnames. should it subset to just the coverage over the single range of interest?

- a standard place for scaling factors: another slot, or a column name in the bigWigSamples DataFrame

- allow BigWigFile to specify bigWigPaths. BigWigFile has methods seqinfo() and summary() which make use of the information contained in the BigWig. E.g. then you can warn about conflicts between bigWigRanges and the BigWig files.

- use BigWigSelection argument to import.bw


