#!/bin/bash
#SBATCH -A p31220
#SBATCH -p short
#SBATCH -t 03:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=28
#SBATCH --job-name=SI_local

julia HBS_Model.jl