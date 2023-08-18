
# covXtreme - repository for the penalised piecewise constant covariate marginal and conditional extreme value models that allows for contour estimation

## Introduction

Piecewise_Covariate_Extremes (PCE) model and software for estimation of environmental design contours using the conditional extremes model of Heffernan and Tawn
[2004]. The sample is composed of peaks over threshold values for both a conditioning variate and its associated conditioned variates. Each pair is allocated to a particular covariate bin; all (joint)
observations with the same covariate bin are assumed to have common extreme value characteristics. The non-stationary marginal extreme value characteristics of each variate is estimated using
roughness-penalised maximum likelihood estimation using a generalised Pareto (GP) model above the threshold and gamma below. The extremal dependence structure between the variates on a transformed standard scale (Gumbel or Laplace) is then estimated using a conditional extremes model, also piecewise non-stationary with respect to covariates. Different approaches to contour estimation,
generally reliant on simulation under the fitted models, are outlined.

More theoretical details of the model can be found in https://www.sciencedirect.com/science/article/abs/pii/S0029801819303798. This work was part-funded by the European Union ERA-NET project entitled “Environmental Contours for SAfe DEsign of Ships and other marine structures (ECSADES). The code is all written in MATLAB. 

## Repository overview

```text
main
├─── develop
```

* `master` is the main branch of the code. 
* `develop` is where any bug fixes and testing will take place.

## User guide

More details of the code can be found in the PCE_User_Guide pdf document: https://github.com/sede-x/Piecewise_Covariate_Extremes/blob/master/PCE_UserGuide.pdf.

## Contact

If you require any more information, please contact ross.towe{at}.shell.com. Please include *Piecewise_Covariate_Extremes* in the subject of the email. 

## Current status
* Assessment of the environmental contours.
* Improving the fit of the Heffernan and Tawn model.

## Feedback
This repository is maintained by David Randell (@davidrandell84), Emma Ross (ERoss0), Philip Jonathan (@ygraigarw) and Ross Towe (@RPTowe). 

