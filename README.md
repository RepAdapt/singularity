# singularity
## Background
This repo contains [Apptainer](https://apptainer.org) (formerly Singularity) commands and image locations for RepAdapt genotyping tools. 

Singularity commands in this repo are based on an existing installation of Apptainer. Installation instructions can be [found here](https://apptainer.org/docs/admin/main/installation.html). It is also possible, and probably eaiser, to install these tools with [conda/mamba](https://github.com/conda-forge/miniforge). 

## Installation and environment activation
A singularity installation is possible using conda/mamba as well, either directly, as:
```
mamba create -n apptainer apptainer -c conda-forge
```
Or as part of the [Nextflow nf-core](https://nf-co.re/docs/nf-core-tools/installation)

```
mamba create -n nf-core nf-core -c bioconda
```

And activated with:
```
mamba activate apptainer
```
Or
```
mamba activate nf-core
```

## RepAdapt images
Singularity/apptainer images can now be downloaded and used for bioinformatic analysis.
Commands for pulling images of all RepAdapt tools are in [RepAdaptSingularity.md](https://github.com/RepAdapt/singularity/blob/main/RepAdaptSingularity.md)
