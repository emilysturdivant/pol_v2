
#!/bin/bash

#SBATCH --partition=ALL
#SBATCH --nodes=1
#SBATCH --mincpus=30
#SBATCH --mem=60G
#SBATCH --time=10-15:15:00     # 60 days, 15 hours and 15 minutes
#SBATCH --output=%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emilysturdivant@gmail.com
#SBATCH --job-name="Testing"
#SBATCH -e err.%j

source activate R410

srun Rscript /share/Part1/esturdivant/pol_v2/R/02_process_SDM.R

#sleep(10)

echo "Termine"