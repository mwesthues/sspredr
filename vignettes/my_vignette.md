---
title: "Single Step Prediction"
author: "Matthias Westhues"
date: "2017-09-04"
output:
  rmarkdown::html_vignette:
    fig_caption: yes
    theme: readable
highlight:tango
vignette: >
  %\VignetteIndexEntry{Single Step Prediction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


The main functionality of `sspred` entails:

1.   Convert between digit- and string-representations of SNPs.
2.   Check for SNP heterozygosity, minor allele frequency and call frequency.
3.   Create ETA objects for single-predictor predictions using
     [BGLR](https://github.com/gdlc/BGLR-R).
4.   Create ETA objects for single-step prediction using BGLR.


## SNP conversion

```r
library("sspredr")
data("imp_snps")

# Return the names of all markers with a heterozygosity of less than 10%.
het_loci_nms <- compute_het(
  x = imp_snps,
  output = "marker_names",
  het_threshold = 0.1
  )
```



```r
# Return the heterozygosity at each locus.
loci_het_rate <- compute_het(
  x = imp_snps,
  output = "marker_heterozygosity",
  het_threshold = 0
  )
```
