# Georeferenced individual tree measurements and condition observations of African baobabs and comparator trees.

This respository contains the data and code used to process and analyse community-contributed Baobab in Africa Mapped (BAM) dataset. It provides the resources needed to reproduce the quality control checks (QAQC) procedures, the figures, and tables presented in the following publication: 

> Chafika Phiri, Andrew M. Cunliffe, Sarah Venter, Nicholas Bii, Aisha T. Bello, Francis Chitsa, Dawn Cory-Toussaint, Kenneth F. Egbadzor, Tapiwa Gumbo, Dave Hodgson, Michelle Hurley, Steven Khosa, Alan D. Nare, Justice Muvengwi, Okopi Ade, Adler Santana Medeiros, Lena Grieger. (Year). *Paper title*. Journal. DOI: https://doi.org/...

A permanent version of this repository is archived at [![DOI](https://zenodo.org/badge/DOI/...svg)](...)

Use of this code is licensed under [![Licence: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](...)

Contact Chafika Phiri at Email: cp844@exeter.ac.uk

---

## Project overview

The Baobab in Africa Mapped (BAM) project aims to create an open-source, community-contributed dataset on the distribution and condition of African baobab trees. This dataset provides a valuable foundation for geospatial investigations of baobab distribution and condition, enabling assessments of their responses to environmental drivers, comparisons with co-occurring species, and evaluations of conservation effectiveness. Furthermore, it contributes to improving the representation of African landscapes within openly accessible ecological datasets. The dataset is harmonised, quality checked, and made available for public use.

---

## Repository contents

### Scripts

| Script | Description |
|---------|-------------|
| *Harmonise_tree_data.R* | Reads individual contributor files and harmonises their column names and variables against the BAM variable dictionary. |
| *TRY_comparison_check.R* | Identifies potential duplicate records by cross-referencing all *Adansonia* observations in the TRY database (Version 7) against the harmonised BAM baobab records. |
| *QAQC_BAM_Data.qmd* | Performs quality checks including validation, exclusion, correction, privacy and publication rules; generates the public dataset. |
| *Environmental_covariates.R* | Extracts environmental covariates for every baobab location, including climate, soil, topography, and protected-area status. |
| *Summary_report_BAM.qmd* | Generates summary statistics, spatial distribution maps, and visualisations to provide an overview of the BAM datase. |


### Data

| File | Description |
|------|-------------|
| *BAM_data_for_public.csv* | Cleaned dataset for public release. |
| *BAM_variable_dictionary.csv* | Defines each variable, its unit or permitted values and its quality-control limits. |
| *baobab_environmental_covariates.csv* | Contains environmental covariates for each unique baobab location. |
| *crowdsourced_data_collection_protocols* | Data collection guidelines for the crowdsource campaigns. |

---

## Acknowledgements

We gratefully acknowledge financial support from Few and Far Collection within the Oppenheimer Programme in African Landscape Systems (OPALS) (https://opals-exeter.org/), which made this research possible. We sincerely thank landowners at Kremetart, Bahiti, ZZ2, Francois Pauer, South African National Parks (SANParks) and all their neighbours for granting permission to conduct fieldwork on their properties. 

## Getting started

Clone this repository and review the project overview and repository structure before running the analysis. 

TESS Lab projects typically use [renv](https://rstudio.github.io/renv/) to record package dependencies and software versions. 
Where an renv.lock file is included, restore the project environment before running:

```r
renv::restore()
```

## Running the analysis

Scripts should be run in numerical order:

1. *Harmonise_tree_data.R* 
2. *TRY_comparison_check.R*
3. *QAQC_BAM_Data.qmd*
4. *Environmental_covariates.R*
5. *Summary_report_BAM.qmd*

---
