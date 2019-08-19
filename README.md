
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build
Status](https://travis-ci.org/SMAC-Group/avar.svg?branch=master)](https://travis-ci.org/SMAC-Group/avar)
[![AppVeyor build
Status](https://ci.appveyor.com/api/projects/status/github/SMAC-Group/avar?branch=master&svg=true)](https://ci.appveyor.com/project/SMAC-Group/avar)
[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![packageversion](https://img.shields.io/badge/Package%20version-0.1.0-orange.svg?style=flat-square)](commits/develop)
[![Last-changedate](https://img.shields.io/badge/last%20change-2019--08--19-yellowgreen.svg)](/commits/master)

# Welcome to `avar` package <a href="https://smac-group.com/"><img src="man/figures/logo.png" align="right" style="width: 16%; height: 16%"/></a>

This package provides the tools necessary to compute the empirical Allan
Variance (AV) and use it to estimate the parameters of (latent) time
series models. The estimation of the Allan Variance is performed through
the estimator proposed by Allan (1966) and, based on this quantity, the
Allan Variance Linear Regression (AVLR) approach is often used by
engineers to retrieve the parameters of time series models which are
assumed to underlie the observed signals (see for example Guerrier,
Molinari, and Stebler 2016). These estimators are implemented in this
package along with the relevant plotting and summary functions.

To see what `simts` is capable of, please refer to the “Vignettes” tabs
above.

## Install Instructions

To install the `simts` package, there is currently one option:
[GitHub](https://github.com/SMAC-Group/av/).

### Installing the package through GitHub

For users who are interested in having the latest developments, this
option is ideal. Though, more dependancies are required to run a stable
version of the package. Most importantly, users **must** have a compiler
installed on their machine that is compatible with R (e.g. Clang).

*The setup to obtain the development version of `simts` is platform
dependent.*

### Requirements and Dependencies

**OS X**

Some users report the need to use X11 to suppress shared library errors.
To install X11, visit [xquartz.org](http://www.xquartz.org/).

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
devtools::install_github("SMAC-Group/avar")

# Install the package with Vignettes/User Guides 
devtools::install_github("SMAC-Group/avar", build_vignettes = TRUE)
```

<div id="refs" class="references">

<div id="ref-allan1966statistics">

Allan, David W. 1966. “Statistics of Atomic Frequency Standards.”
*Proceedings of the IEEE* 54 (2). IEEE: 221–30.

</div>

<div id="ref-guerrier2016theoretical">

Guerrier, Stéphane, Roberto Molinari, and Yannick Stebler. 2016.
“Theoretical Limitations of Allan Variance-Based Regression for Time
Series Model Estimation.” *IEEE Signal Processing Letters* 23 (5). IEEE:
597–601.

</div>

</div>
