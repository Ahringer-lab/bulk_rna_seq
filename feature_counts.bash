#!/bin/bash
#SBATCH --job-name=CHIPseq
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=12
#SBATCH --mem=50gb
#SBATCH --partition=2004
#SBATCH --output=pipeline_%j.log # Standard output and error log

##########################################################################################################################################################################
############################## feature counts sbatch script froma sample sheet #############################################################################
# This code will gather all the bam file names listed in the bam sample sheet from the input folder into an array to initiate pipeline runs for each one in the HPC
# Feature counts is an alternative to kallisto in the main pipeline
# The script creates a parent directory with the run ID, this is a date/time stamp unless specified as an option
# The script will create a folder for each bam file with the bam id as it's name within the parent folder
# Remember to change the SBATCH options above to configure for your run
# This script should only be run on the HPC using sbatch
# Options include:
#      threads = Will multi-thread any process to this number - Not currently used
#      input = Change the path of the input bam files (from star), default is ~/data
#      id = Change the name of the parent pipeline output folder, the default is a datestamp
# This only uses the default featurecounts settings, options can be added as required
# Author Steve Walsh May 2024
###########################################################################################################################################################################

#Uncomment below for debugging
#set -x

#Set the defaults
outdir=~/out
bam_dir=~/data/
GTF_FILE=/mnt/home3/ahringer/index_files/annotation_files/c.elegans.PRJNA13758.WS285_cannonical/c_elegans.PRJNA13758.WS285.canonical_geneset.gtf
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
WD="$(pwd)"

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --threads <number of threads> --input <input path> --id <Run ID>"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l threads: -l input: -l id: -- "$@") || exit_with_bad_args

#Get the inputs
eval set -- "$options"
while true; do
    case "$1" in
        --threads)
            shift
            THREADS="$1"
            ;;
        --input)
            shift
            bam_dir="$1"
            ;;
        --id)
            shift
            RUNID="$1"
            ;;
         --)
            shift
            break
            ;;
    esac
    shift
done

analysis_out_dir=${outdir}/${RUNID}
mkdir ${analysis_out_dir}/
mkdir ${analysis_out_dir}/Feature_counts

cp sample_sheet_fc.csv ${analysis_out_dir}/sample_sheet_fc_${RUNID}.csv
INPUT=${analysis_out_dir}/sample_sheet_fc_${RUNID}.csv

while IFS= read -r LINE 
do

    # split line into array using , delimitator
    #echo ${var1}
    ARRAYLINE=(${LINE//,/ })
    BAM_FILE=${ARRAYLINE[0]}

    #Make directories for the feature call fileds as we go
    mkdir ${analysis_out_dir}/Feature_counts/${BAM_FILE}
    cd ${analysis_out_dir}/Feature_counts/${BAM_FILE}

    srun featureCounts -p -O -T ${THREADS} -a ${GTF_FILE} -o featureCounts_${BAM_FILE} ${bam_dir}/${BAM_FILE}_Aligned.sortedByCoord.out.bam &

done < ${INPUT}

#Wait for all pipelines to finish
wait
