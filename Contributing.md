<!--
SPDX-FileCopyrightText: 2023 Shell Global Solutions International B.V. All Rights Reserved.

SPDX-License-Identifier: Apache-2.0
-->
# Contributing to covXtreme

First off, thanks for taking the time to contribute!

The following is a set of guidelines for contributing to *covXtreme* and its packages. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request. For further information about how the repository please the governance page: [Governance.md](https://github.com/sede-open/covXtreme/blob/main/Governance.md) 

## What should I know before I get started?

*covXtreme* is model and software for estimation of environmental design contours using the conditional extremes model of Heffernan and Tawn [2004].

More theoretical details of the model can be found in https://www.sciencedirect.com/science/article/pii/S1364815224000963. This work was part-funded by the European Union ERA-NET project entitled “Environmental Contours for SAfe DEsign of Ships and other marine structures (ECSADES). The code is all written in MATLAB.

The code is constantly developing and we will ensure to maintain that it is up to date and make clear any changes that may have an impact of the end user.


#### Package Conventions

There are a few conventions that have developed over time around packages:

* The repository has the following structure
    * *Code* this is where the MATLAB code fits that is called in the *RunScripts* folder. This code folder will need to be on your MATLAB path. 
    * *RunScripts* this folder contains a number of run scripts that can be used as a basis to run the code. 

* *RunScripts/Stages_** contains a set of scripts that enable us to complete a full analysis:
    * *Stage1_SimulateData* - optional stage that allows you to simulate example data.
    * *Stage1_PeakPicking* - select of storm peaks from the main variable and selection of associated values.
    * *Stage2_SetBinEdges* - set the bin edges for the piecewise model.
    * *Stage3_FitMargin* - fit the marginal model: piecewise constant extreme value analysis model.
    * *Stage4_PeakPicking* - fit the dependence model: Heffernan and Tawn.
    * *Stage5_Contour* - estimation of environmental contours.

## How Can I Contribute?

### Getting Started
In order to starting contributing you can fork the repository and make any changes to the code. Once you are satisfied with the changes made to the code, they can be merged into the main branch of the repository through a pull request. More details can be found here: https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project

### Potential Enhancements 
* Multivariate diagnostics such as the coefficient of tail dependence e.g., eta, chi and chibar.
* Choice of optimal bin edges - a potential approach would be to use reversible-jump Markov chain Monte Carlo.
* Implementation of alternative contour methods. 
* Reimplementation of *covXtreme* in Python.  

### Reporting Bugs/Issues and Suggesting Enhancements

This section guides you through submitting an enhancement suggestion or an issue with the *covXtreme* code, including completely new features and minor improvements to existing functionality. Please submit an issue request with a detailed explanation on the current problems you are facing.

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
