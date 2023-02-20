
# Piecewise_Covariate_Extremes - repository for the penalise piecewise constant covariate marginal and conditional extreme value models that allows for contour estimation

## Introduction

Piecewise_Covariate_Extremes (PCE) model and software for estimation of environmental design contours using the conditional extremes model of Heffernan and Tawn
[2004]. The sample is composed of peaks over threshold values for both a conditioning variate and its associated conditioned variates. Each pair is allocated to a particular covariate bin; all (joint)
observations with the same covariate bin are assumed to have common extreme value characteristics. The non-stationary marginal extreme value characteristics of each variate is estimated using
roughness-penalised maximum likelihood estimation using a generalised Pareto (GP) model above the threshold and gamma below. The extremal dependence structure between the variates on a transformed standard scale (Gumbel or Laplace) is then estimated using a conditional extremes model, also piecewise non-stationary with respect to covariates. Different approaches to contour estimation,
generally reliant on simulation under the fitted models, are outlined.

More theoretical details of the model can be found in https://www.sciencedirect.com/science/article/abs/pii/S0029801819303798. This work was part-funded by the European Union ERA-NET project entitled “Environmental Contours for SAfe DEsign of Ships and other marine structures (ECSADES). The code is all written in MATLAB. 

## Repository overview

```text
master
├─── develop
```

* `master` is the main branch of the code. 
* `develop` is where any bug fixes and testing will take place.. 

## User guide

More details of the code can be found in the PPC_User_Guide pdf document: https://github.com/sede-x/Piecewise_Covariate_Extremes/blob/master/PPC_UserGuide.pdf.

## Contact

If you require any more information, please contact ross.towe{at}.shell.com. Please include *Piecewise_Covariate_Extremes* in the subject of the email. 

## Current status
* Assessment of the environmental contours.
* Improving the fit of the Heffernan and Tawn model.

## Feedback
This repository is maintained by David Randell (@davidrandell84), Emma Ross (ERoss0), Philip Jonathan (@ygraigarw) and Ross Towe (@RPTowe). They are all a part of the MetOcean sede-x group on GitHub. 


# Contributing to Piecewise_Covariate_Extremes

First off, thanks for taking the time to contribute!

The following is a set of guidelines for contributing to *Piecewise_Covariate_Extremes* and its packages. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

## What should I know before I get started?

*Piecewise_Covariate_Extremes* is model and software for estimation of environmental design contours using the conditional extremes model of Heffernan and Tawn [2004].

More theoretical details of the model can be found in https://www.sciencedirect.com/science/article/abs/pii/S0029801819303798. This work was part-funded by the European Union ERA-NET project entitled “Environmental Contours for SAfe DEsign of Ships and other marine structures (ECSADES). The code is all written in MATLAB.

The code is constantly developing and we will ensure to maintin that it is up to date and make clear any changes that may have an impact of the end user.


#### Package Conventions

There are a few conventions that have developed over time around packages:

* The repository has the following structure
    * *Code* this is where the MATLAB code fits that is called in the *UnitTests* folder. There is a separate subfolder the different versions of the Heffernan and Tawn model this allows us to switch between the different assumed distributions for the residual distribution. This folder will need to be on your MATLAB path. 
    * *Paper* this folder contains the user guide that steps through an example application and the stages needed to be undertaken to complete an analysis.
    * *UnitTests* this folder contains a number of unit tests that have been used to determine the robustness of the code. 

* *UnitTests/Stages_** contains a set of scripts that enable us to complete a full analysis:
    * *Stage1_SimulateData* - optional stage that allows you to simulate example data.
    * *Stage1_PeakPicking* - select of storm peaks from the main variable and selection of associated values.
    * *Stage2_SetBinEdges* - set the bin edges for the piecewise model.
    * *Stage3_FitMargin* - fit the marginal model: piecewise constant extreme value analysis model.
    * *Stage4_PeakPicking* - fit the dependence model: Heffernan and Tawn.
    * *Stage5_Contour* - estimation of environmental contours.

## How Can I Contribute?

### Reporting Bugs/Issues and Suggesting Enhancements

This section guides you through submitting an enhancement suggestion or an issue with the *Piecewise_Covariate_Extremes* code, including completely new features and minor improvements to existing functionality. Please submit an issue request with a detailed explanation on the current problems you are facing.

Before creating enhancement suggestions, please include as much detail as possible and any steps the steps that you imagine would be needed for the enhancement to be included into the code. Given budget and time we cannot guarantee that every enhancement will be incorporate into the code. 

Please start each item with the following handle:
* An enhancement: @ENHANCEMENT 
* A bug: @BUG
* An issue: @ISSUE

This will help to identify high priorities problems that 

### Pull Requests

The process described here has several goals:

- Ensure that the code maintains its standards and is not affected by any new additions. 
- Fix problems that are important to users.

### Git Commit Messages

* Use the present tense. 
* Use the imperative mood. 
* Reference issues and pull requests liberally after the first line.

