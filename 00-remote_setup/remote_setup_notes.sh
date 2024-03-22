# on palmetto

#PBS -N env_create
#PBS -l select=1:ncpus=2:mem=50gb
#PBS -l walltime=2:00:00
#PBS -m abe

module add anaconda3/2023.1
conda env create -f long_oral_microbiome/environment.yml

# next activate the environment on palmetto and install R kernel for jupyter notebook
conda activate 2023-Long_oral_microbiome
R -e 'IRkernel::installspec()'
# install stringi
R -e 'install.packages("stringi", repos="http://cran.us.r-project.org")'









# remote notebook on pickles?
ssh -L 8154:localhost:8154 allie@pickles.clemson.edu

# start jupyter notebook on remote machine
conda activate 2023-Long_oral_microbiome
jupyter-notebook
