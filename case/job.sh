#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=naca_opt
#SBATCH --mail-user=yukuns@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=10g
#SBATCH --time=36:00:00
#SBATCH --account=yulinpan1
#SBATCH --partition=standard
##SBATCH --begin=2020-11-26T15:00:00

# The application(s) to execute along with its input arguments and options:
#module load openfoam/v1906
source ~/foam2006
./Allrun
