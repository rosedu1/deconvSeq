# deconvSeq

R package for performing cell type deconvolution of bulk RNA sequencing, single cell RNA sequencing, and bisulfite sequencing data.

## Installation

The user can download the [tar ball](https://github.com/rosedu1/deconvSeq/tree/master/tarball/current/), and run `R CMD INSTALL` on it, or use the **devtools** package to install from GitHub:

```r
## devtools is required
library(devtools)
install_github("rosedu1/deconvSeq")
```

Note: 
1) Windows users need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and [devtools](http://CRAN.R-project.org/package=devtools) to install from GitHub.
2) Windows users need to install **simpleSingleCell** separately from the source as there is no binary for Windows

```r
source("http://bioconductor.org/biocLite.R")
biocLite("simpleSingleCell", type = "source")
```

3) If some packages do not install, set repositories to include both CRAN and Bioconductor with `setRepositories(ind=c(1:4))`.


## Help

Vignette: [HTML Vignette](https://rosedu1.github.io/deconvSeq/deconvSeq_vignette.html)

## Version notes
5/20/2020 version 0.1.3 1) in withDoc.R, function prep_scrnaseq, replaced defunct functions calculateQCMetrics with perCellQCMetrics and replaced defunct calcAverage with calculateAverage. 2) in function getmethmat, added filetypes to include output from Bismark as well as BSMAP, 3) in function getmethmat, added condition to not filter if there are no cases where all Cs or all Ts are 0.

## Contact

You are welcome to:
* submit suggestions and bug-reports at: <https://github.com/rosedu1/deconvSeq/issues>
* email: <rdu@bwh.harvard.edu>
