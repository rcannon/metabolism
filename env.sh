#!/bin/bash

module purge
module load gcc/9.1.0
module load cmake/3.21.4
module load python/miniconda3.9
source /share/apps/python/miniconda3.9/etc/profile.d/conda.sh

# activate the local conda environment located in 
# metab/python/anaconda
export PROJECT_DIR=$(pwd)
conda activate $PROJECT_DIR/python/anaconda