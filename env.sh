#!/bin/sh
module purge
module load gcc/9.1.0
module load python/miniconda3.8
source /share/apps/python/miniconda3.8/etc/profile.d/conda.sh

# if the following fails, make sure spack is installed:
# git clone https://github.com/spack/spack.git ~/spack
source spack/share/spack/setup-env.sh

spack env activate /qfs/people/cann484/metab



