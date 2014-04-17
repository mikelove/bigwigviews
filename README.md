**Note**: this functionality is part of Bioconductor package `GenomicFiles` as of April 2014

Code is hosted at https://github.com/Bioconductor/GenomicFileViews

---

The BigWig paper is here:

http://bioinformatics.oxfordjournals.org/content/26/17/2204.long

# Design questions

- BamViews offers parallelization over samples using srapply() in .BamViews_delegate(). We are focusing with BigWigViews on chunking over ranges instead. But we have to worry about I/O. Are we assuming users put bigwigs in scratch?

- coverage on a BigWigViews with a single range currently gives the RleList of all seqnames. should it subset to just the coverage over the single range of interest?

- a standard place for scaling factors: another slot, or a column name in the bigWigSamples DataFrame

- allow BigWigFile to specify bigWigPaths. BigWigFile has methods seqinfo() and summary() which make use of the information contained in the BigWig. E.g. then you can warn about conflicts between bigWigRanges and the BigWig files.

- use BigWigSelection argument to import.bw

- check the seqlevels!!!
