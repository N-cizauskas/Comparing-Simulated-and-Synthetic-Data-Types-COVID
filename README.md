# Comparing-Simulated-and-Synthetic-Data-Types


## Description

This contains R code for creating a simulation trial of COVID-19 vaccine data and corresponding synthetic datasets including analysis and graphs. 


Simulated clinical trial data was generated for three different types for trials (randomized control trial, observational study, and external data) based on probabilities from COVID-19 data. Data comes from February-August 2021 in the United Kingdom. The outcome being looked at is testing positive for COVID-19, and the treatment/exposure variable of interest is the BNT162b2 Pfizer-BioNTech vaccine. Ethnicity and sex are included as potential confounders in the RCT and observational study. Simulated data is tested through a logisitc regression model to confirm previously published findings.  References for original datasets can be found below.

Synthetic data was generated using the R library Synthpop. Synthetic data is created by subsetting the original data to include only variables of interest, and then creating a cookbook of values.  Three methods of synthetic data generation were used, each corresponding to its own R script: CART (Categorical and Regression Tree modelling), random sampling, and linear/logistic regression. More information on how Synthpop generates synthetic data can be found here:
https://www.synthpop.org.uk/about-synthpop.html#methodology

The simulated and synthetic data are compared between data types and sample sizes by measuring the treatment effect using the chi-squared analysis.  The standard mean difference was also calulcated and graphed.

All_Sim_Synth_Script.R contains the pipeline for synthetic data generated using the CART method.
All_Sim_Synth_Script_Rand.R contains the pipeline for synthetic data generated using the random sampling method.
All_Sim_Synth_Script_LL.R contains the pipeline for synthetic data generated using the linear/logisitc regression method.


## Installation

The script can be downloaded and run in R Studio.

## Libraries

The R Libraries used include the following:


- tidyverse
- sandwich
- stargazer
- truncnorm
- tableone
- survival
- ggplot2
- ggmap
- table1
- lubridate
- plyr
- reshape
- MASS
- reshape2
- synthpop
- dplyr
- stddiff
- Matching
- survey
- ggstance
- hrbrthemes
- viridis
- rgp
- rsimsum



## Usage

This code is for comparing the quality of synthetic controls between different data types.  The simulated probabilties can be altered to fit new datasets for further testing.

## Citations

### More information on Synthpop Synthesis Methods

- Nowok, B., Raab, G.M. and Dibben, C. (2016) ‘synthpop: Bespoke Creation of Synthetic Data in R’, Journal of Statistical Software, 74, pp. 1–26. Available at: https://doi.org/10.18637/jss.v074.i11.


### Randomized Control Trial Original Data

- Demographic data for coronavirus (COVID-19) testing (England): 28 May to 26 August (no date) GOV.UK. Available at: https://www.gov.uk/government/publications/demographic-data-for-coronavirus-testing-england-28-may-to-26-august/demographic-data-for-coronavirus-covid-19-testing-england-28-may-to-26-august (Accessed: 8 March 2024).
- Yang, Z.-R. et al. (2023) ‘Efficacy of SARS-CoV-2 vaccines and the dose–response relationship with three major antibodies: a systematic review and meta-analysis of randomised controlled trials’, The Lancet Microbe, 4(4), pp. e236–e246. Available at: https://doi.org/10.1016/S2666-5247(22)00390-1.


### Observational Study Original Data

- Bernal, J.L. et al. (2021) ‘Early effectiveness of COVID-19 vaccination with BNT162b2 mRNA vaccine and ChAdOx1 adenovirus vector vaccine on symptomatic disease, hospitalisations and mortality in older adults in England’. medRxiv, p. 2021.03.01.21252652. Available at: https://doi.org/10.1101/2021.03.01.21252652.
- Demographic data for coronavirus (COVID-19) testing (England): 28 May to 26 August (no date) GOV.UK. Available at: https://www.gov.uk/government/publications/demographic-data-for-coronavirus-testing-england-28-may-to-26-august/demographic-data-for-coronavirus-covid-19-testing-england-28-may-to-26-august (Accessed: 8 March 2024).


### External Original Data
- Cases in England | Coronavirus in the UK (2023). Available at: https://coronavirus.data.gov.uk/details/cases?areaType=nation&areaName=England (Accessed: 28 February 2024).
- Demographic data for coronavirus (COVID-19) testing (England): 28 May to 26 August (no date) GOV.UK. Available at: https://www.gov.uk/government/publications/demographic-data-for-coronavirus-testing-england-28-may-to-26-august/demographic-data-for-coronavirus-covid-19-testing-england-28-may-to-26-august (Accessed: 8 March 2024).
- Statistics » COVID-19 vaccinations archive (no date). Available at: https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/covid-19-vaccinations-archive/ (Accessed: 28 February 2024).

