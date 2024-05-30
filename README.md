# Stability of the oral microbiome in children with HIV

[![DOI](https://zenodo.org/badge/776010265.svg)](https://zenodo.org/doi/10.5281/zenodo.11396311)

This repository describes the analysis performed in the paper *Mann et al. //////*. Raw read processing steps and all figure generation/statistical analysis scripts are found in the folders in this repository. If you are trying to recreate any of the scripts/analyses described in the paper and run into any trouble, please open an issue!

NOTE: UNDER CONSTRUCTION 🚧

## Setup 

This repository assumes you are running a Unix environment (e.g., Mac OSX or Linux) and have conda installed.

To get this repository:

- Install and set up anaconda or miniconda as described at the [bioconda
  documentation](https://bioconda.github.io/user/install.html), including
  setting up channels.
- Clone this repository to your machine and change into the directory with:


```bash
git clone https://github.com/aemann01/long_oral_microbiome && cd long_oral_microbiome
```

- Run the following to install the environment

```bash
conda create --name 2023-Long_oral_microbiome python=3.7
# activate the environment
conda activate 2023-Long_oral_microbiome
# install the following packages
conda install -c bioconda cutadapt=1.18
conda install -c bioconda seqtk=1.4
conda install -c bioconda kraken2=2.1.2
conda install -c bioconda mafft=7.520
conda install -c bioconda fasttree=2.1.11
conda install -c conda-forge r-base=4.1.0
conda install -c conda-forge r-irkernel=1.3.2
conda install -c conda-forge jupyter=1.0.0

# install jupyter notebook with pip
pip install notebook
```

- Add the R kernel to Jupyter by installing a kernel spec

```bash
R -e 'IRkernel::installspec()'
# install stringi
R -e 'install.packages("stringi", repos="http://cran.us.r-project.org")'
```

- To deactivate the environment

```bash
conda deactivate
```
