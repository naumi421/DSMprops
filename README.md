
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DSMprops: Digital Soil Mapping of numeric properties

<!-- badges: start -->

<!-- badges: end -->

A hybrid repository including an R package for commonly used functions
and project executable scripts (exec folder) documenting current
workflows for project implementation of DSM models in the Soil Landscapes
of the United States (SOLUS) 100m product.

## Installation

You can install the released version of DSMprops from GitHub for R
functions with:

``` r
install.packages("remotes")
remotes::install_github("naumi421/DSMprops")
```

To get both R package functionality and access to all exec workflows,
you can clone this repository locally and then build R package. Although
you can just clone the repository to get all the code and be able to
build the package, devtools is needed to contribute to the repository.
Note that [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
(windows) or other pre-requisites (e.g.
[Ubuntu](https://www.digitalocean.com/community/tutorials/how-to-install-r-packages-using-devtools-on-ubuntu-18-04))
are necessary for installing devtools.

Cloning in bash:

``` bash
$ cd /path/where/repo/goes
$ git clone https://github.com/naumi421/DSMprops.git (or alternate repo URL)
```

Rstudio also has options for cloning a repo with good guidance available
at <https://happygitwithr.com/>.

Once repository has been cloned. Open the R project in Rstudio and use
*Build \> Install and Restart* to install DSMprops package which can
then be used in other sessions outside of the project.

# DSMprops project overview:

The aim is to provide code guidance for all steps of Unites States
National Cooperative Soil Survey digital soil mapping workflows. We
break these out by certain tasks needed to acheive the
mapping.

## Properties included in SOLUS100 dataset

| Property                 | Units       | Notes                                            |
| ------------------------ | ----------- | ------------------------------------------------ |
| Clay                     | % mass      |                                                  |
| Silt                     | % mass      |                                                  |
| Sand                     | % mass      |                                                  |
| Very fine sand           | % mass      |                                                  |
| Fine sand                | % mass      |                                                  |
| Medium sand              | % mass      |                                                  |
| Coarse sand              | % mass      |                                                  |
| Very coarse sand         | % mass      |                                                  |
| Depth to bedrock         | cm          | Depth to lithic of paralithic contact            |
| Depth to restriction     | cm          | Depth to layer that limits air, water, or roots  |
| pH                       | -log[H+]    | 1:1 method                                       |
| Rock fragments           | % vol.      | Particles >2mm, whole soil basis                 |
| Bulk Density             | g cm<sup>-3</sup> | Oven Dry                                   |
| Calcium carbonate        | % mass      |                                                  |
| CEC at pH7               | meq / 100g  | Cation Exchange Capacity (CEC) <2mm fraction     |
| ECEC                     | meq / 100g  | Effective CEC                                    |
| Soil organic carbon      | % mass      |                                                  |
| Sodium Adsorption Ratio  | unitless    | Ratio of [Na]/([Ca]+[Mg]) in Saturated Paste     |
| Electrical conductivity  | dS / m      | Saturated paste                                  |


# Semi-detailed Structure of Project.

<!-- create url links to each part when repo location is finalized -->

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
