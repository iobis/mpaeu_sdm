Habitat range maps for European marine species
================

## About the project



NOTE: This project is still under development. 

Due to limits of file size, the raw data is not uploaded. However, you can obtain all data using the codes named `pre_download_*.R`

## Directory structure

    ├── README.md              : Description of this repository
    ├── LICENSE                : Repository license
    ├── mpaeu_sdm.Rproj : RStudio project file
    ├── .gitignore             : Files and directories to be ignored by git
    │
    ├── data
    │   ├── raw                : Source data obtained from repositories (e.g. OBIS, GBIF)
    │   ├── species            : Processed species data
    │   │   └── key=*          : Folder for species * (AphiaID)
    │   │       └──date=*      : Folder for download/version date * (YMD format)
    │   ├── shapefiles         : Shapefiles
    │   └── environmental      : Environmental data
    │       ├── current        : Data for current period
    │       └── future         : Data for future period (a folder for each scenario)
    │
    ├── codes                  : All codes
    │
    ├── functions              : Functions used in the project
    │
    ├── models                 : SDM models for each species 
    │                            (format: key={aphiaID}model={algorithm}date={YMDdate}.rds)
    │
    ├── results                : Results for the SDMs
    │                            (format: key={aphiaID}/model={algorithm}/...)
    │
    ├── docs                   : Repository website generated
    │
    └── src                    : Repository website source

## Main codes

Most of the components of the SDM framework are provided through the [`obissdm` package](https://github.com/iobis/mpaeu_msdm). The codes and functions of this repository only _operationalize_ the modelling.

### Codes (in order of execution)

- check\_\*: check species occurring in an area, etc.
- get\_\*: obtain data for something.
- prepare\_\*: prepare the data to be used.
- model\_\*: modeling the species’ distribution.
- predict\_\*: make predictions.
- plot\_\*: plotting.

### Functions

## Basic components

## Running models for the full list of species

```
Rscript source_run.R
```
