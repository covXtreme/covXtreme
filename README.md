<!--
SPDX-FileCopyrightText: 2023 Shell Global Solutions International B.V. All Rights Reserved.

SPDX-License-Identifier: Apache-2.0
-->
# covXtreme - repository for the penalised piecewise constant covariate marginal and conditional extreme value models that allows for contour estimation

## Background

Hazard risk analysis often involves modelling extreme events, for example of natural phenomena such as rainfall, temperature and ocean waves, or of man-made systems such as industrial processes and stock-markets.  We seek to understand the characteristics of the most exceptional (largest or most rare) events we have seen historically, and potentially might see if we look harder.

In the world of statistical modelling, the central limit theorem provides us with solid reasoning to model the average of a sample of data with a Normal (or Gaussian) distribution. For extremes of data, mathematical theory tells us that extreme value distributions are most appropriate to use. This enables us to construct statistical models for extreme events, and to perform quantitative risk analysis using them. Hence, in a natural hazard context, we can design flood defences and sea walls to pre-specified performance standards. 

Underlying factors, referred to as covariates, typically affect the behaviour of large events. For example in the case of extreme ocean waves, important covariates include wave direction and season. To characterise the behaviour of extremes adequately, it is often essential to incorporate non-stationarity with respect to covariates. The covXtreme software provides a way to achieve this. Similarly, we are often interested in extreme events across multiple variables, not just one. So we need a model for the joint distribution of extremes of one or more variables, again varying with covariates. covXtreme provides a way to do this.

covXtreme was developed for oceanographic applications, but it has been used more widely. You will probably find covXtreme interesting, if:
* You are interested in quantifying extremes using statistical analysis of a data set, and might want to estimate extreme quantiles or “return levels”,
* The characteristics of your data are not steady; that is, they vary with covariates which you know about,
* You may be interested in quantifying extremes of multiple variables at the same time, and
* You want to quantify how confident you can be in your modelling.
  
covXtreme will then allow you to perform a pragmatic but statistically sound analysis. 

More motivation behind the development of covXtreme can be found in the following blog: [Characterising extreme environments](https://medium.com/data-centric-engineering-blog/characterising-extreme-environments-cc97b2403fcb).

## Introduction

The covXtreme model and software for estimation of environmental design contours using the conditional extremes model of Heffernan and Tawn
[2004]. The sample is composed of peaks over threshold values for both a conditioning variate and its associated conditioned variates. Each pair is allocated to a particular covariate bin; all (joint)
observations with the same covariate bin are assumed to have common extreme value characteristics. The non-stationary marginal extreme value characteristics of each variate is estimated using
roughness-penalised maximum likelihood estimation using a generalised Pareto (GP) model above the threshold and gamma below. The extremal dependence structure between the variates on a transformed standard scale (Gumbel or Laplace) is then estimated using a conditional extremes model, also piecewise non-stationary with respect to covariates. Different approaches to contour estimation,
generally reliant on simulation under the fitted models, are outlined.

More theoretical details of the code and respective models can be found in the following submitted paper: [covXtreme: MATLAB software for non-stationary penalised piecewise constant marginal and conditional extreme value models](https://arxiv.org/abs/2309.172950). This work was part-funded by the European Union ERA-NET project entitled “Environmental Contours for SAfe DEsign of Ships and other marine structures (ECSADES). The code is all written in MATLAB. 

## Repository overview

```text
main
├─── develop
```

* `main` is the main branch of the code. 
* `develop` is where any bug fixes and testing will take place.

## User guide

More details of the code can be found in the covXtreme user guide document see [covXtreme_UserGuide.pdf](https://github.com/sede-open/covXtreme/blob/main/covXtreme_UserGuide.pdf) in the main repository folder.

## Issues and suggestions

Please see the [Contributing.md](https://github.com/sede-open/covXtreme/blob/main/Contributing.md) in the main folder of the GitHub repository. 

## Contact

If you require any additional information, please contact ross.towe{at}shell.com. Please include *covXtreme* in the subject of the email. 

## Feedback
This repository is maintained by David Randell (@davidrandell84), Emma Ross (ERoss0), Philip Jonathan (@ygraigarw) and Ross Towe (@RPTowe). 

