# deconvSeq

R package for performing cell type deconvolution of bulk RNA sequencing, single cell RNA sequencing, and bisulfite sequencing data.

## Installation

The user can download the [tar ball](https://github.com/rosedu1/deconvSeq/tree/master/tarball/current/), untar the tarball, and use the **devtools** package to install dependencies then run `R CMD INSTALL` on the command line:

```r
#in R
library(devtools)
install_deps("/PATH/TO/deconvSeq",dependencies=TRUE)
#on command line
R CMD INSTALL deconvSeq_VERSION_NUMBER.tar.gz
```

Alternatively, use the **devtools** package to install from GitHub:

```r
## devtools is required
library(devtools)
install_github("rosedu1/deconvSeq", dependencies=TRUE)
```

Note: 
1) Windows users need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and [devtools](https://CRAN.R-project.org/package=devtools) to install from GitHub.
2) If some packages do not install, set repositories to include both CRAN and Bioconductor with `setRepositories(ind=c(1:4))`.


## Help

Vignette: [HTML Vignette](https://rosedu1.github.io/deconvSeq_vignette.html)

## Version notes
5/20/20 version 0.1.3 1) in withDoc.R, function prep_scrnaseq, replaced defunct functions calculateQCMetrics with perCellQCMetrics and replaced defunct calcAverage with calculateAverage. 2) in function getmethmat, added filetypes to include output from Bismark as well as BSMAP, 3) in function getmethmat, added condition to not filter if there are no cases where all Cs or all Ts are 0.

5/20/20 version 0.1.4 1) Updated mydiffMethFromDesign_matrix.R, line 169, class(res)=="list" changed to class(res)[1]=="list". This does not affect the functionality of the package.

5/21/20 version 0.2.0 1) Reduced size of data_scrnaseq.rda to reduce package size. 2) Updated prep_scrnaseq to include error message for wrong genenametype.

5/23/20 version 0.2.1 1) Removed dependency on simpleSingleCell from package

5/25/20 version 0.2.2 1) Removed examples from getcorr

5/25/20 version 0.2.3 1) Commented out time consuming portions of vignette and added precompiled data

5/25/20 version 0.2.4 1) Added citation to vignette and package

6/17/20 version 0.2.5 1) Reduced vignette run time


## Citation
Du R, Carey V, and Weiss ST. deconvSeq: deconvolution of cell mixture distribution in sequencing data, Bioinformatics, 35:5095-5102, 2019.

## Contact

You are welcome to:
* submit suggestions and bug-reports at: <https://github.com/rosedu1/deconvSeq/issues>
* email: <rdu@bwh.harvard.edu>
