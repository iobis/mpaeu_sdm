Habitat range maps for European marine species
================

## About the project

This work is part of the [MPA Europe project](https://mpa-europe.eu/). OBIS is leading the WP3, which aims to generate distribution maps for marine species and habitats in Europe. This repository contains the code for generating the SDMs (species distribution models) and stacked SDMs (habitat maps).

The core functions behind our modelling framework are in the repository [iobis/mpaeu_msdm](https://github.com/iobis/mpaeu_msdm)(from 'methods' SDM), which contains the package `obissdm`. A more detailed documentation of our framework can be found [here](https://iobis.github.io/mpaeu_docs).

Maps are available through a Shiny application, accessible from [shiny.obis.org/distmaps](https://shiny.obis.org/distmaps). Codes used to produce the platform are also open, and are available through [this repository](https://github.com/iobis/mpaeu_map_platform).

NOTE: This project is still under development. 

## Replicating the project

This GitHub repository contains only the codes/functions. You can clone it to your computer, and download the remaining data from Amazon cloud and from the sources (e.g. Bio-ORACLE, OBIS and GBIF). The easiest way is to run sequentially the 5 main codes (named `p*_`) which organizes the steps.

- **p1_prepare_wd.R**: will install the requirements (R and Python packages) and check the folder structure  
- **p2_download_data.R**: download all necessary data  
- **p3_prepare_data.R**: used to prepare the data for modelling (including standardization and QC)  
- **p4_model_distribution.R**: fit the models and make predictions. Code is run in parallel
- **p5_stack_habitats.R**: stack the SDMs from different groups to produce habitat maps

## Directory structure

Cloning this repository and running the main codes will render a directory with the following structure:


    ├── README.md              : Description of this repository
    ├── LICENSE                : Repository license
    ├── mpaeu_sdm.Rproj        : RStudio project file
    ├── .gitignore             : Files and directories to be ignored by git
    ├── requirements.R         : Project requirements
    ├── check.R                : Check project structure
    ├── sdm_conf.yml           : Configuration file for the models
    │
    ├── data
    │   ├── raw                : Source data obtained from repositories (e.g. OBIS, GBIF)
    │   ├── distances          : Distances with barriers, used for QC steps
    │   ├── log                : Log objects
    │   ├── species            : Processed species data
    │   ├── shapefiles         : Shapefiles
    │   └── environmental      : Environmental data
    │       ├── current        : Data for current period
    │       ├── terrain        : Data for terrain variables (e.g. bathymetry)
    │       └── future         : Data for future period (a folder for each scenario)
    │
    ├── codes                  : All codes
    │
    ├── functions              : Functions used in the project
    │
    ├── results                : Results for the SDMs - see details below
    │
    └── analysis               : Short analysis done over the project

## Main codes

As already noted, most of the components of the SDM framework are provided through the [`obissdm` package](https://github.com/iobis/mpaeu_msdm). The codes and functions of this repository only _operationalize_ the modelling.

Codes are all commented, with additional information provided on headers. In general, codes follow this naming convention:

- check\_\*: check species occurring in an area, etc.
- get\_\*: obtain data for something.
- prepare\_\*: prepare the data to be used.
- model\_\*: modeling the species’ distribution.
- plot\_\*: plotting.
- pre_tests\_\*: tests with virtual species.

The code named `model_subset.R` enable you to pass a subset of species for modelling (instead of the full list). If you want to run models for just a subset of species, don't run the code `p4_model_distributions.R` and run this one instead.

## Running models for the full list of species

Ensure that the working directory is correctly built. From **the root of the working directory** run the 3 first steps:

``` bash
Rscript codes/p1_prepare_wd.R
Rscript codes/p2_download_data.R
Rscript codes/p3_prepare_data.R
```
Then, run the 4th step to obtain the models:

``` bash
Rscript codes/p4_model_distributions.R
```

Or, to run for just a subset, change the `model_subset.R` file and then run it:

``` bash
Rscript codes/model_subset.R
```
## Results structure

The results are organized as:

`taxonid={aphiaID}/model={acronym of model run}/<folder> OR <file>`

With folders being 'figures', 'metrics', 'models' or 'predictions'

All files will contain `taxonid={aphiaID}_model={acronym of model run}` as part of their name.

Two files are saved on the root: 'taxonid={aphiaID}_model={acronym of model run}_what=fitocc.parquet', which contain the points used for model fitting, and 'taxonid={aphiaID}_model={acronym of model run}_what=log.json', a log file containing rich details about model fitting.

## Additional information

Some steps are controlled using [`storr`](https://richfitz.github.io/storr/). This will create `*_storr` folders, which you can later delete.

If you experience any trouble using this resource, contact helpdesk@obis.org