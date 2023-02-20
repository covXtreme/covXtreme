# ACTIONS

Paper writing => Phil
- Admin: git for paper development, template for article
- Sketch out paper

Analyses for paper => Phil
- Find three good case studies

Updating user guide => Ross/Emma
- What has been changed / added since user guide was written?

Updating code => Ross
- Header comments need to be standardised
- GetSummary explosed

# PAPER

- Env Mod Soft https://www.sciencedirect.com/journal/environmental-modelling-and-software

- Abstract
  = Simple practically-useful approach to non-stationary marginal and conditional extreme value analysis with uncertainty quantification
- Introduction : PPC in context (metocean, environment generally, other software)
  = Modelling context
  = Marginal and joint extremes
  = Covariates
  = Splines, Voronoi, Simon Wood, Ben Youngman
  = Paul Northrop review extremes software!
  = Earlier software reviews
  = Competitor software
  = What is the PPC USP?
  = Surge paper, ECSADES review paper
  = Extensions (PPC, PPL) with Ed
  = Other users: Kaust papers (need to search)
  
- PPC methodology
  = Marginal model
  = Conditional model
  = Emphasise novel bits
    : DL residuals
    : Backfitting for parameter estimation
    : Improved diagnostics 
    : Importance sampling in place of/along side MC simulation
      :: Reference demanning paper, ... there's one other paper with import samp in it (approx 2018)
    : Estimation of return values and associated values
      :: Which options are available?
      :: Which options have been thoroughly tested?

- PPC software
  = High level overview
  = User guide for details

- Case study 1 : Hs-Tp
  = Introductory case: Hs-Tp, single covariate

- Case study 2 : non-metocean
  = hydrology
  = temperature with hour of day, or season
  = rainfall with season
  = needs covariates
  = Emma Eastoe atmospheric pollutants

- Case study 3 : multivariate response, 2D covariate
  = Fancy: OTM, WS, Tp | Hs, directional-seasonal

- Discussion
  = Mention Hs, WS, Tp | Response
  = Mention spatial conditional extremes
  = Mention MEM / heatwave, Stan Amsterdam paper
  = Mention multivaraite MEM, EVAR, Stan recent work

# SOFTWARE

- code
  = code is basically ready
  = header comments need to be standardised; a mess at present, but easy to improve
  = create new repo

- user guide
  = this is in really good state, so we just need to check it's up to date
  = what's not up to date? Ross, what have you modified or added?
 
# BACKGROUND

Resources in Shell
- Z: hot storage
- Y: cold storage
- PPC testing: Z:\project\MetOcean\PPC (this is where all the case studies are)
- Code: \git\MetOcean_PPC_CE (this is also C:\Users\Philip.Jonathan\Git\MetOcean_PPC_CE)
- User Guide Y:\project\Extremes\ECSADES\PPC\PPC_Reportv2 (this is also )
- ECSADES Review paper: Y:\project\Extremes\ReportsPapersConferences\2019_ECSADESReviewPaper

Anything else:
Y:\project\Extremes\ECSADES

Key selling point : simplest practically useful approach with some accommodation of nonstationarity

Application
- Introductory case: Hs-Tp, single covariate
- Fancy: OTM, WS, Tp | Hs, directional-seasonal
- Hs, WS, Tp | Response





