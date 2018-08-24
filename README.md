
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/SMAC-Group/av.svg?branch=master)](https://travis-ci.org/SMAC-Group/av) [![Project Status: Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Licence](https://img.shields.io/badge/licence-CC%20BY--NC--SA%204.0-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN](http://www.r-pkg.org/badges/version/av)](https://cran.r-project.org/package=av) [![packageversion](https://img.shields.io/badge/Package%20version-0.1.0-orange.svg?style=flat-square)](commits/develop) [![Last-changedate](https://img.shields.io/badge/last%20change-2018--08--24-yellowgreen.svg)](/commits/master) <a href="https://twitter.com/intent/follow?screen_name=SMAC_Group"> <img src="https://img.shields.io/twitter/follow/SMAC_Group.svg?style=social&logo=twitter"
  alt="follow on Twitter"></a>

Welcome to `av` package <a href="https://smac-group.com/"><img src="man/figures/logo.png" align="right" style="width: 16%; height: 16%"/></a>
=============================================================================================================================================

This package provides the tools necessary to compute the empirical Allan Variance (AV) and use it to estimate the parameters of (latent) time series models. The estimation of the Allan Variance is based on the estimator proposed by Allan (1966) and, based on this quantity, the Allan Variance Linear Regression approach is often used by engineers to retrieve the parameters of time series models which are assumed to underlie the observed signals (see for example Guerrier, Molinari, and Stebler 2016). These estimators are implemented in this package along with the relevant plotting and summary functions.

To see what `simts` is capable of, please refer to the "Vignettes" tabs above.

Install Instructions
--------------------

To install the `simts` package, there is currently one option: [GitHub](https://github.com/SMAC-Group/av/).

### Installing the package through GitHub

For users who are interested in having the latest developments, this option is ideal. Though, more dependancies are required to run a stable version of the package. Most importantly, users **must** have a compiler installed on their machine that is compatible with R (e.g. Clang).

*The setup to obtain the development version of `simts` is platform dependent.*

### Requirements and Dependencies

**OS X**

Some users report the need to use X11 to suppress shared library errors. To install X11, visit [xquartz.org](http://www.xquartz.org/).

**Linux**

Both curl and libxml are required.

For **Debian** systems, enter the following in terminal:

``` bash
sudo apt-get install curl libcurl3 libcurl3-dev libxml2 libxml2-dev
```

For **RHEL** systems, enter the following in terminal:

``` bash
sudo yum install curl curl-devel libxml2 libxml2-dev
```

``` r
# Install dependencies
install.packages(c("devtools"))

# Install/Update the package from GitHub
devtools::install_github("SMAC-Group/av")

# Install the package with Vignettes/User Guides 
devtools::install_github("SMAC-Group/av", build_vignettes = TRUE)
```

Allan, David W. 1966. “Statistics of Atomic Frequency Standards.” *Proceedings of the IEEE* 54 (2). IEEE: 221–30.

Guerrier, Stéphane, Roberto Molinari, and Yannick Stebler. 2016. “Theoretical Limitations of Allan Variance-Based Regression for Time Series Model Estimation.” *IEEE Signal Processing Letters* 23 (5). IEEE: 597–601.
