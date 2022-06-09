# drb-gw-hw-model-prep
Code repo to prepare groundwater and headwater-related datasets for modeling river temperature in the Delaware River Basin

This repo contains a targets pipeline for compiling datasets and a snakemake workflow for extracting simulated groundwater discharge from a MODFLOW groundwater model

## Extracting the catchment / reach attributes
The scripts to compile the catchment attributes utilize an R targets pipeline that is intialized with the "_targets.r" file in the main directory.


## Extracting the MODFLOW outputs
The scripts to run the MODFLOW extraction are contained within the directory "MODFLOW_extraction." They are python scripts and utilize a Snakemake workflow.

### To run the Snakemake workflow locally:

1. Install the dependencies in the `environment_MODFLOW_CONUS.yaml` file. With conda you can do this with `conda env create -f environment_MODFLOW_CONUS.yaml`
2. Activate your conda environment `source activate mf_data_extraction`
3. Edit the run configuration (including paths for I/O data) in the appropriate `config.yml` (either `config_MODFLOW_DRB.yml` or `config_MODFLOW_CONUS.yml`)
4. Run Snakemake with `snakemake --configfile config_MODFLOW_CONUS.yml -s Snakefile_MODFLOW_CONUS -j1`

### To run the Snakemake Workflow on TallGrass
1. Request a GPU allocation and start an interactive shell

        salloc -N 1 -t 4:00:00 -A <account>
        srun -A <account> --pty bash

2. Follow steps 1-4 above as you would to run the workflow locally. 

___

### Disclaimer
This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the official USGS copyright policy

Although this software program has been used by the U.S. Geological Survey (USGS), no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided “AS IS.”