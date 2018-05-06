
<!-- README.md is generated from README.Rmd. Please edit that file -->
updog
=====

[![Build Status](https://travis-ci.org/dcgerard/updog.svg?branch=master)](https://travis-ci.org/dcgerard/updog) [![Appveyor Build status](https://ci.appveyor.com/api/projects/status/c80fdy61ead6s3vr?svg=true)](https://ci.appveyor.com/project/dcgerard/updog-06s8t) [![codecov](https://codecov.io/gh/dcgerard/updog/branch/master/graph/badge.svg)](https://codecov.io/gh/dcgerard/updog) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Updog provides a suite of methods for genotyping autopolyploids from next-generation sequencing (NGS) data. It does this while accounting for many common features of NGS data: allelic bias, overdispersion, sequencing error, and (possibly) outlying observations. It is named updog for "Using Parental Data for Offspring Genotyping" because we originally developed the method for full-sib populations, but it works now for more general populations.

The main function is `flexdog`, which provides many options for the distribution of the genotypes in your sample.

The experimental function `mupdog` allows for correlation between the individuals' genotypes while jointly estimating the genotypes of the individuals at all provided SNPs. The implementation uses a variational approximation. This is designed for samples where the individuals share a complex relatedness structure (e.g. siblings, cousins, uncles, half-siblings, etc).

The original `updog` package is now names `updogAlpha` and may be found [here](https://github.com/dcgerard/updogAlpha).

See also [ebg](https://github.com/pblischak/polyploid-genotyping), [fitPoly](https://cran.r-project.org/package=fitPoly), and [TET](http://www.g3journal.org/content/suppl/2017/01/19/g3.117.039008.DC1). Our best "competitor" is probably [fitPoly](https://cran.r-project.org/package=fitPoly) :smile:.

Installation
------------

You can install updog from Github with:

``` r
# install.packages("devtools")
devtools::install_github("dcgerard/updog")
```

If you want to use the `use_cvxr = TRUE` option in `flexdog` (not generally recommended), you will need to install the [CVXR](https://cran.r-project.org/package=CVXR) package. Before I could install CVXR in Ubuntu, I had to run in the terminal

``` bash
sudo apt-get install libmpfr-dev
```

and then run in R

``` r
install.packages("Rmpfr")
```

How to Cite
-----------

If you use `flexdog`, please cite

> Gerard, D., Ferrão L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Harnessing Empirical Bayes and Mendelian Segregation for Genotyping Autopolyploids from Messy Sequencing Data. *bioRxiv*. doi: [10.1101/281550](https://doi.org/10.1101/281550).

Or, using BibTex:

``` tex
@article {gerard2018harnessing,
    author = {Gerard, David and Ferr{\~a}o, Luis Felipe Ventorim and Garcia, Antonio Augusto Franco and Stephens, Matthew},
    title = {Harnessing Empirical Bayes and Mendelian Segregation for Genotyping Autopolyploids from Messy Sequencing Data},
    year = {2018},
    doi = {10.1101/281550},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2018/03/16/281550},
    eprint = {https://www.biorxiv.org/content/early/2018/03/16/281550.full.pdf},
    journal = {bioRxiv}
}
```

Otherwise, please cite

> David Gerard (2018). mupdog: Extension of the updog Method for Genotyping Autopolyploids. R package version 0.0.2.

Or, using BibTex

``` tex
  @Manual{gerard2018mupdog,
    title = {mupdog: Extension of the updog Method for Genotyping Autopolyploids},
    author = {David Gerard},
    year = {2018},
    note = {R package version 0.0.2},
  }
```

Code of Conduct
---------------

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
