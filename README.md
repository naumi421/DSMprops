
# DSMprops: Digital Soil Mapping of numeric properties

<!-- badges: start -->

<!-- badges: end -->

A hybrid repository including an R package for commonly used functions
and project executable scripts (exec folder) documenting current
workflows for project implementation of DSM models.

## Installation

You can install the released version of DSMprops from GitHub for R
functions with:

``` r
# install.packages("devtools")
devtools::install_github("naumi421/DSMprops")
```

Note that [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
(windows) or other pre-requisites (e.g.
[Ubuntu](https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-18-04))
are necessary for installing devtools. Alternatively, you can clone this
repository locally and then build R package to have access to both
package function calls and exec workflows.

In bash:

``` bash
cd /path/where/repo/goes
git clone naumi421/DSMprops.git
```

Rstudio also has options for cloning a repo with good guidance available
at <https://happygitwithr.com/>.

Once repository has been cloned. Open the R project in Rstudio anduse
*Build \> Install and Restart* to install DSMprops package which can
then be used in other sessions outside of the project..

# DSMprops project overview:

The aim is to provide code guidance for all steps of Unites States
National Cooperative Soil Survey digital soil mapping workflows. We
break these out by certain tasks needed to acheive the
mapping.

## The properties we want to map

| Property                 | Units       | Notes                                            |
| ------------------------ | ----------- | ------------------------------------------------ |
| Clay                     | g/kg or %   |                                                  |
| Sand                     | g/kg or %   |                                                  |
| Silt                     | g/kg or %   |                                                  |
| Depth                    | cm          |                                                  |
| pH                       | unitless    |                                                  |
| rock fragments           | (m3 m-3)    |                                                  |
| Available water capacity | (mm mm-1 ?) | probably need to run Rosetta or Saxton and Rawls |

# Potential properties to be predicted

| Property                | Units      | Notes |
| ----------------------- | ---------- | ----- |
| ECEC                    | (cmolc/kg) |       |
| N                       | (%)        |       |
| C                       | (%)        |       |
| electrical conductivity | (dS/m)     |       |

# Semi-detailed Structure of Project.

## A. Pedon dataset cleaning and prep \[per property?\]

### A.1 Forest Service NRM data \[Colby\]

### A.2 NASIS point data

### A.3 KSSL data

### A.4 RaCA

### A.5 Other point data

## B. Depth Harmonization

### (0-5, 5-15, 15-30, 30-60, 60-100, 100-200 cm) \[per property\]

## C. Build Training Matrix

### C.1 Extract DEM derivatives

### C.2 Extract spectral derivatives

### C.3 Extract other geospatial covariates

## D. Build Models \[per property\]

### D.1 Build list to hold each regional model per property \[big point of discussion here\]

### D.2 Set up training routine (cross validation repeated 10 times or something).

### D.3 Training/evaluation data split

### D.4 Build continental model

### D.5 Build regional models

## E. Model Evaluation

### E.1 Evalute continental model

### E.2 Evaulate each regional model

### E.3 Covariate importance

## F. Spatial Predictions \[how to use regional models?\]

## G. Spatial Uncertainty

## H. Locations where the addition of legacy or new sample data would prove most beneficial to model improvement and/or uncertainty reduction.
