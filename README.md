
<!-- README.md is generated from README.Rmd. Please edit that file -->

## About

This repository contains R code and data to reproduce the data
processing, analysis, and visualization for the study *“Many non-native
plant species are threatened in parts of their native range.”*

While all data are included, the Red List file
(`rl_combined_REDLIST_updated.xlsx`) required for
`01a-data-prep-rl-exploration.R` is available upon request.

## Repository Structure

- **`R-code/`**: Contains seven R scripts for data processing and
  analysis.

  - **`00-preamble.R`**: Loads necessary packages and custom functions.
  - **`01a-data-prep-rl-exploration.R`**: Explores the Red List
    synthesis provided by Laura Méndez Cuéllar, this products is already
    taxonomically harmonized and country names standardized.
  - **`01b-data-prep-rl-naturalization.R`**: Merges Red List and GloNaF
    naturalization data.
  - **`02-magnitude-paradox.R`**: Quantifies the number of non-native
    plants threatened or near-threatened and examines the influence of
    data completeness and neutral factors.
  - **`03-data-analysis-gbif-range.R`**: Quantifies the Area of
    Occupancy (AOO) of threatened non-native species, distinguishing
    native, threatened, and non-native AOOs.
  - **`04-data-analysis-phylogeny.R`**: Assesses whether threatened
    non-native species are disproportionately represented in certain
    plant families.
  - **`05-data-analysis-geography.R`**: Analyzes spatial patterns of
    threatened non-native species.

- **`Data/`**: Contains all input, output, and intermediary data files.

- **`Figures/`**: Contains all figures generated in the analysis.

## Contact

Please contact me at <ingmar.staude@uni-leipzig.de> if you have further
questions.
