bigwigviews
===========

R version of desired BigWigViews function

The BigWig paper is here:

http://bioinformatics.oxfordjournals.org/content/26/17/2204.long

For the moment, I just grabbed the necessary classes, generics and methods from BamViews, and added:

```
setMethod(coverage, "BigWigViews",
```

which makes the import() and coverage() calls. This might already not work over multiple chromosomes.
