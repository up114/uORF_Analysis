#!/bin/bash
#SBATCH -J myjob         # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p normal              # Queue name (use normal, development, etc.)
#SBATCH -N 1                   # Number of nodes
#SBATCH -n 1                   # Number of MPI tasks
#SBATCH -t 6:00:00             # Expected runtime (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=umapaul@outlook.com

# Load required modules
module spider python # Load Python module (adjust version as needed)

# Activate a virtual environment if necessary
conda init --all
conda activate my_env


# Change to the directory where your Python script is located
cd /scratch/09369/umapaul/scripts

# Run your Python script
python uorf_reads_psite.py

# Deactivate the virtual environment if you activated one
deactivate
