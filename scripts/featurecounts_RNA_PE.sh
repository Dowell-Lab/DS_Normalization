#!/bin/bash
#SBATCH --job-name=featurecounts # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=user@email.com # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=16 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=00:30:00 # Time limit hrs:min:sec
#SBATCH -p short
#SBATCH --mem=4gb # Memory limit
#SBATCH --output=/path/to/e_o/%x_%j.out
#SBATCH --error=/path/to/e_o/%x_%j.err

module load python/3.6.3
module load subread

while test $# -gt 0; do
           case "$1" in
                -h|--help)
                    echo "Uses FeatureCounts to count paired reads over input regions"
                    echo " "
                    echo "options:"
                    echo "-h, --help, show brief help"
                    echo "-b|--bamdir, BAM directory file"
                    echo "-r|--refsaf, SAF file to count over"
                    echo "-n|--name, BAM prefix"
                    echo "-o|--outfile, output file for counts" 
                    exit 0
                    ;;
                -b|--bamdir)
                    shift
                    BAMDIR=$1
                    shift
                    ;;
                -n|--name)
                    shift
                    NAME=$1
                    shift
                    ;;
                -r|--refsaf)
                    shift
                    REFSAF=$1
                    shift
                    ;;
                -o|--outfile)
                    shift
                    OUTFILE=$1
                    shift
                    ;;
                *)
                   echo "$1 is not a recognized flag!"
                   return 1;
                   ;;
          esac
  done

printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

declare -a BamFiles
for file in `ls ${BAMDIR}/${NAME}*bam`
do
    BamFiles=("${BamFiles[@]}" "$file")
done

function logr {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

#logr "Parsed Params: ""$regionFile"

logr "Starting Analysis"

NUM_CORES=8
safFile=$REFSAF
countsOut=${OUTFILE}.sense.txt

featureCounts \
                -T "$NUM_CORES" \
                -s 1 \
		-O \
                -f \
		-p \
                --fracOverlap 0.51 \
                -F 'SAF' \
                -a "$safFile" \
                -o "$countsOut" \
                ${BamFiles[@]} 


countsOut=${OUTFILE}.antisense.txt
featureCounts \
                -T "$NUM_CORES" \
                -s 2 \
		-O \
                -f \
		-p \
                --fracOverlap 0.51 \
                -F 'SAF' \
                -a "$safFile" \
                -o "$countsOut" \
                ${BamFiles[@]}

countsOut=${OUTFILE}.unstranded.txt
featureCounts \
                -T "$NUM_CORES" \
                -s 0 \
                -O \
                -f \
		-p \
                --fracOverlap 0.51 \
                -F 'SAF' \
                -a "$safFile" \
                -o "$countsOut" \
                ${BamFiles[@]} 
