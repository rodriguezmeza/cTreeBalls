#!/bin/bash
#SBATCH --ntasks 2 # was 32
#SBATCH -J tpcf_n2_10millon
#SBATCH --array=1 #Run 28 copies of the code # was --array=1-2
#SBATCH -o standard/standard_output_file.%A.%a.out
#SBATCH -e standard/standard_error_file.%A.%a.err
#SBATCH -p cosma8
#SBATCH -A dp203
#SBATCH -t 12:00:00 # Instead of 72:00:00

source $HOME/.bashrc
module load intel_comp/2018a
module load gsl
# Define the maximum number of threads to use
export OMP_NUM_THREADS=64

tpcf sizeHistN=40 lbox=10000 verb=3 mToPlot=2 rangeN=100 nbody=100000000 search=treeomp

# submit as
# sbatch tpcf.sl
#
# squeue |grep dc-rodr6
# squeue -lu dc-rodr6

