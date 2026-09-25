[![DOI](https://zenodo.org/badge/315139777.svg)](https://zenodo.org/badge/latestdoi/315139777) [![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIt)
[![GitHub contributors](https://img.shields.io/github/contributors/michaelwenk/sherlock.svg)](https://github.com/michaelwenk/sherlock/graphs/contributors/) [![GitHub issues](https://img.shields.io/github/issues/michaelwenk/sherlock.svg)](https://github.com/michaelwenk/sherlock/issues/) [![GitHub release](https://img.shields.io/github/release/michaelwenk/sherlock.svg)](https://github.com/michaelwenk/sherlock/releases/)

<img width="150" alt="sherlock-logo" src="public/Sherlock.png" align="right">

# Sherlock

Web services for Computer-Assisted Structure Elucidation (CASE).

A [graphical user interface](https://github.com/michaelwenk/sherlock-frontend) and [publication](https://doi.org/10.3390/molecules28031448) are available.

## Core Features

- Dereplication
- Elucidation
- Retrieval of previously generated results

The dereplication, chemical shift prediction and statistical detection of structural constraints are enabled by using entries from spectral knowledge bases, such as [NMRShiftDB](https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2withsignals.sd), containing both structural and spectral properties.

[casekit](https://github.com/michaelwenk/casekit) serves as computational library in Sherlock.

[PyLSD](https://github.com/nuzillard/PyLSD/blob/4db027118acf3d9c77d3a5f8dc9ce51dd3cbd468/pylsd-linux-a8.tar.gz) is used for structure generation.

## Docker and Databases

This project uses Docker containers (https://www.docker.com) and starts them via docker compose. Make sure that docker compose is installed.

NOTE: It is recommended to set the accessible RAM to 8 GB or higher and the number of available CPU cores to four or higher. For example, this can be done in the Docker Desktop application, see [here](/public/Docker_settings.png).

### Download

Clone this repository and change the directory:

    git clone https://github.com/michaelwenk/sherlock.git && \
    cd sherlock && \
    cp env.dist .env

### Steps for Data Import and Statistics

In order to fill the database with datasets and to create different statistics, databases in specific SD format, e.g. [NMRShiftDB](https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2withsignals.sd), can be imported into Sherlock's backend system.

The directory _data/nmrshiftdb_ (default) may contain the NMRShiftDB file but ending with _.sdf_, such as _nmrshiftdb.sdf_. The directory _data/acd_labs_predictions_ (default) may contain database with predicted spectra produced by ACD/Labs. Each file should not contain more than 125.000 entries, due to potential in-memory issues. Multiple files with subsets are possible, e.g. from [LOTUS](https://doi.org/10.5281/zenodo.8175939) or COCONUT, e.g. _lotus_1.sdf_, _lotus_2.sdf_ etc.

The directory _data/lsd/filters_ (default) may contain molecular fragments in LSD format, such as the filter examples in the _data/lsd/PyLSD/LSD/Filters_ folder, which then will be applied as structural constraints (badlist) to not allow such fragments during the structure generation process.

Uncomment und modify _DATASET_NMRSHIFTDB_PATH_, _DATASET_ACD_LABS_PREDICTIONS_PATH_ or _CUSTOM_FILTERS_PATH_ in the _.env_ file if another folder than the default ones should be mounted and used by Docker.

After starting via docker compose, run the following script to delete previous database records, load-in the spectral data and to build the statistics.

    sh scripts/fill_DBs.sh

After the dataset import, this script continues with building different statitics, such as about hybridisation states or HOSE codes. The progress of that process can be tracked via console logging or via a API endpoint with login data, which can be adjusted in the _.env_ file:

    curl -X 'GET' 'http://localhost:8080/database/buildStatistics/status' -u 'admin:password'

### Start

Start the docker compose with pre-built images in detached mode:

    docker compose up -d

Or to re-build the image, if not already given then download [PyLSD](https://github.com/nuzillard/PyLSD/blob/4db027118acf3d9c77d3a5f8dc9ce51dd3cbd468/pylsd-linux-a8.tar.gz) and make sure that the unpacked folder is located in the _data/lsd/_ folder and is re-named to "_PyLSD_". Copy the two modified files from the _data/lsd/PyLSD_Variant_mod_ folder into the _data/lsd/PyLSD/Variant_ folder and overwrite the two existing files there.

To build and start the services in detached mode use:

    docker compose up -d --build

The following URL leads to the Swagger UI and offers an overview of what REST endpoints are available, both public or internal (private with login only) ones:

    http://localhost:8080/swagger-ui.html

To access the generated OpenAPI JSON file just visit or fetch:

    http://localhost:8080/api-docs

### Stop

To shutdown the services and to remove the docker compose network use:

    docker compose down -v

###
