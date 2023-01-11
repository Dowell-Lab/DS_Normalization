#!/bin/bash 
#SBATCH --job-name=sync # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --partition=long
#SBATCH --mem=15gb # Memory limit
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/sync.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/sync.%j.err # Standard error log




INDIR=/Shares/down/genome/read_mapping_hg38/
WKDIR=/scratch/Shares/dowell/down/temp_genomehg38/
rsync -r $INDIR $WKDIR
