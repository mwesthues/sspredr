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


## Genotype encoding
Depending on the technology for calling SNP reads and depending on the data
format that you're using for storing SNP data you will have entries in your
SNP matrix that are either coded as strings, numeric values or integer values.
The function `recode_snps` can translate between the string- and the
number-representation while taking into account which value encodes a major or
a minor allele at any locus.

First determine the major and minor allele at each locus of SNP-genotypes.
We set the argument 'maf_threshold' to zero because at this point we do not
want to apply a quality check for minor allele frequency.


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
data("marker_character")

# Check the original data format.
typeof(marker_character)
```

```
## [1] "character"
```

```r
print(marker_character[, 1:5])
```

```
##       col1 col2 col3 col4 col5
## row1  "BB" "BB" "??" "BB" "AA"
## row2  "BB" "BB" "AA" "BB" "BB"
## row3  "AA" "AA" "BB" "BB" "BB"
## row4  "??" "AB" "AA" "AB" "BB"
## row5  "BB" "BB" "AA" "AA" "BB"
## row6  "AA" "BB" "AA" "??" "AA"
## row7  "AA" "AA" "BB" "BB" "BB"
## row8  "BB" "AA" "??" "??" "AA"
## row9  "AA" "AA" "BB" "AB" "BB"
## row10 "AA" "BB" "BB" "BB" "AA"
```

```r
geno_lst <- compute_maf(
  marker_character,
  output = "geno_list",
  missing_value = "??",
  maf_threshold = 0
  )
```

Next, we want to recode the major allele as `2`, the minor allele as `0` and
heterozygotes as `1`.


```r
marker_numeric <- sspredr::recode_snps(
  marker_character,
  major = geno_lst[["major_genotype"]],
  minor = geno_lst[["minor_genotype"]],
  major_coding = 2,
  minor_coding = 0,
  het_coding = 1,
  missing_value = "??",
  na_coding = NA_real_
  )

# Check the new data format.
typeof(marker_numeric)
```

```
## [1] "double"
```

```r
print(marker_numeric[, 1:5])
```

```
##       col1 col2 col3 col4 col5
## row1     0    2   NA    2    0
## row2     0    2    2    2    2
## row3     2    0    0    2    2
## row4    NA    1    2    1    2
## row5     0    2    2    0    2
## row6     2    2    2   NA    0
## row7     2    0    0    2    2
## row8     0    0   NA   NA    0
## row9     2    0    0    1    2
## row10    2    2    0    2    0
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


### SNP heterozygosity
In the following example we have a matrix of SNPs where some loci comprise an
excessively large amount of heterozygosity.
We would like to keep only loci where the heterozygosity rate is less than
say 10%.


```r
data("imp_snps")

# Return the names of all markers with a heterozygosity of less than 10%.
het_loci_nms <- sspredr::compute_het(
  imp_snps,
  output = "marker_names",
  het_threshold = 0.1
  )

low_het_snps <- imp_snps[, colnames(imp_snps) %in% het_loci_nms]
print(paste("The original number of loci was:", nrow(imp_snps)))
```

```
## [1] "The original number of loci was: 90"
```

```r
print(paste("The number of low-heterozygosity loci is:", nrow(low_het_snps)))
```

```
## [1] "The number of low-heterozygosity loci is: 90"
```



```r
# Return the heterozygosity at each locus.
loci_het_rate <- imp_snps %>%
  sspredr::compute_het(
    output = "marker_heterozygosity",
    het_threshold = 0
    )
```






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


