Allele-Specific Expression in a Population (ASEP)
======================

`ASEP`, overcomes three analytical challenges, utilizes cross-subject and cross-SNP information to detect gene-level ASE under one condition and differential ASE between two conditions (e.g., pre- versus post-treatment) in a population using RNA sequencing (RNA-seq) data.

<p align="center"> 
<img src="https://github.com/Jiaxin-Fan/ASEP/raw/master/Figure.png" width="500">
</p>

How to cite `ASEP`
-------------------
Please cite the following publication:

> *ASEP: gene-based detection of allele-specific expression in a population by RNA-seq*<br />
> <small>J. Fan, J. Hu, C. Xue, H. Zhang, M. Reilly, R. Xiao, M. Li<br /></small>
> (https://www.biorxiv.org/content/10.1101/798124v1) 

Installation
------------

``` r
# install devtools if necessary
if (!"devtools" %in% rownames(installed.packages())) {
  install.packages('devtools')
}
# install the ASEP package
if (!"ASEP" %in% rownames(installed.packages())) {
  devtools::install_github('Jiaxin-Fan/ASEP')
}
# load
library(ASEP)
```

More Information
-----------------
Please see [Tutorial](https://jiaxin-fan.github.io/ASEP/articles/introduction.html).

