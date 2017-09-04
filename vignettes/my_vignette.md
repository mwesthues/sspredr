---
title: "Single Step Prediction"
author: "Matthias Westhues"
date: "2017-09-04"
output:
  rmarkdown::html_vignette:
    fig_caption: yes
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
library("dplyr")
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
data("imp_snps")

# Return the names of all markers with a heterozygosity of less than 10%.
het_loci_nms <- imp_snps %>%
  sspredr::compute_het(
    output = "marker_names",
    het_threshold = 0.1
    )
```



```r
# Return the heterozygosity at each locus.
loci_het_rate <- imp_snps %>%
  sspredr::compute_het(
    output = "marker_heterozygosity",
    het_threshold = 0
    )
```


## Missing values
In many cases your 'raw' SNP data will contain missing values.
Since throwing them away would be a considerable waste of your resources you
should strive to impute them.
I recommend to use an imputation software such as
[BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) for this
purpose, but here I've also provided a simple function which imputes missing
values with the major allele at any locus:


```r
data("marker_numeric")

# Examine the data.
print(marker_numeric)
```

```
##       col1 col2 col3 col4 col5 col6 col7 col8 col9 col10
## row1     2    2   NA    2    0    2    0    0    0    NA
## row2     2    2    0    2    2    1    2    0    2     2
## row3     0    0    2    2    2    2    2    0    2     0
## row4    NA    1    0    1    2    0    2    2    0     2
## row5     2    2    0    0    2    0    0    2    1     2
## row6     0    2    0   NA    0    0    1    0    0     2
## row7     0    0    2    2    2    0    0    0    2     0
## row8     2    0   NA   NA    0    0    0    0    1     2
## row9     0    0    2    1    2    2    2    0    0     2
## row10    0    2    2    2    0    0    0    2    2     2
```

```r
marker_numeric %>%
  is.na() %>%
  sum()
```

```
## [1] 6
```

```r
# First, compute the major genotype at each locus.
major_genotype <- marker_numeric %>%
  sspredr::compute_maf(
    x = .,
    output = "geno_list",
    missing = NA_real_,
    maf_threshold = 0
  ) %>%
  .[["major_genotype"]]

# Second, replace all missing genotypes with the missing allele.
imp_numeric <- marker_numeric %>%
  sspredr::replace_na_with_major_genotype(
    dat = .,
    missing_value = NA_real_,
    major_genotype = major_genotype
    )
print(imp_numeric)
```

```
##       col1 col2 col3 col4 col5 col6 col7 col8 col9 col10
## row1     2    2    0    2    0    2    0    0    0     2
## row2     2    2    0    2    2    1    2    0    2     2
## row3     0    0    2    2    2    2    2    0    2     0
## row4     0    1    0    1    2    0    2    2    0     2
## row5     2    2    0    0    2    0    0    2    1     2
## row6     0    2    0    2    0    0    1    0    0     2
## row7     0    0    2    2    2    0    0    0    2     0
## row8     2    0    0    2    0    0    0    0    1     2
## row9     0    0    2    1    2    2    2    0    0     2
## row10    0    2    2    2    0    0    0    2    2     2
```

```r
# Check whether there are still missing values in the marker matrix.
imp_numeric %>%
  is.na() %>%
  sum()
```

```
## [1] 0
```


## SNP quality filtering
In any genomic prediction it is crucial to filter SNPs based on the criteria
call frequency, minor allele frequency and heteroyzgosity rate.
Moreover, particular in the case of small samples of genotypes, the set of SNPs
may be characterized by identical marker loci.
Let's load a set of marker data and recreate such a scenario by creating a
duplicate of a marker locus and appending it to the marker matrix.


```r
data("marker_numeric")
snp <- marker_numeric

# Add a duplicated marker locus.
snp <- cbind(snp, snp[, 1])

# Check whether any marker locus is a duplicate of another locus.
snp %>%
  duplicated(., MARGIN = 2) %>%
  sum()
```

```
## [1] 1
```


