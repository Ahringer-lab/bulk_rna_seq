#!/bin/bash
#SBATCH --job-name=RNASeq  
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=6 
#SBATCH --mem=25gb
#SBATCH --output=pipeline_%j.log # Standard output and error log

###############################################################################################################################################
############################## bulk rna-seq sbatch job submission script ######################################################################
# This code will gather all the fastq file names from the input folder into an array initiate pipeline runs for each one in the HPC
# The pipeline assumes fastq file have been merged using the bash script on github
# The fastq id is the first part of the standard file name afer merging e.g. for JAtab71-F-1_merged_R1_001.fastq.gz the id is JAtab71-F-1
# The script creates a parent directory with the run ID, this is a date/time stamp unless specified as an option
# The script will create a folder for each pair of fastq files with the fastq id as it's name within the parent folder
# Remember to change the SBATCH options above to configure for your run, ntasks should be the number of fastq pairs
# This script should only be run on the HPC using sbatch
# Options include:
#      threads = Will multi-thread any process to this number
#      input = Change the path of the input fastq files, default is ~/data
#      id = Change the name of the output folder, the default is a datestamp
#      mergeID = If the file names have been merged differently the input can be changed here 'fastqid_<Add the flag here>_R1/R2_001.fastq.gz'
#      star_index = The location of the STAR index
#      kallisto_index = The location of the Kallisto index
# Author Steve Walsh May 2024
################################################################################################################################################


#Set the defaults
outdir=~/out
kallisto_index=/mnt/home3/ahringer/index_files/built_indexes/kallisto/c.elegans.full.april2024/c_elegans.PRJNA13758.WS285.canonical_geneset.g.idx
fastq_dir=~/data/
star_index=/mnt/home3/ahringer/index_files/built_indexes/star/c.elegans.full.2024/
CHROM_SIZES=/mnt/home3/ahringer/index_files/genomes/c_elegans.PRJNA13758.WS285.genomic.chrom.sizes
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
MERGEID=merged
WD="$(pwd)"

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --threads <number of threads> --input <input path> --id <Run ID>  --mergeID <merge ID> --star_index --kallisto_index"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l threads: -l input: -l id: -l mergeID -l star_index -l kallisto_index -- "$@") || exit_with_bad_args

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
            fastq_dir="$1"
            ;;
        --id)
            shift
            RUNID="$1"
            ;;
        --mergeID)
            shift
            MERGEID="$1"
            ;;
        --star_index)
            shift
            star_index="$1"
            ;;
        --kallisto_index)
            shift
            kallisto_index="$1"
            ;;
         --)
            shift
            break
            ;;
    esac
    shift
done

#Set and create the ouput directory based on Run ID (Date/time if not set)
analysis_out_dir=${outdir}/${RUNID}
mkdir $analysis_out_dir
echo "$analysis_out_dir"

#Set up stats folder
mkdir ${analysis_out_dir}/stats

MERGEID=_merged

#Set up an array to store fastq file id's
declare -A FILES

cd ~/data

# Gather the fastq files from the input folder and add to FILES array, the ID's are the hash key
for f in *fastq.gz; do                  # search the files with the suffix
    base=${f%${MERGEID}*}                        # remove after "_MERGID_" To make sample ID the hash key
    if [[ $f == $base* ]] && [[ $f == *"R1"* ]]; then    # if the variable is the current sample ID and is forward
        FILES[$base]=$f                  # then store the filename
    elif [[ $f == $base* ]] && [[ $f == *"R2"* ]]; then # if the variable is the current sample and is reverse
        FILES[$base]+=" $f"
    fi
done

#Loops through the fastq names, make directories for each output, ${base} holds the sample name
for base in "${!FILES[@]}"; do
    echo "${base}${MERGEID}_R1_001.fastq.gz"
    echo "${base}${MERGEID}_R2_001.fastq.gz"

cd ${WD}
          ./bulk_rna_seq_pipeline.bash --fastqid ${base} --threads ${THREADS} --input ${fastq_dir} --id ${RUNID} --mergeID ${MERGEID} --star_index ${star_index} --kallisto_index ${kallisto_index} &
done

#Wait for all pipelines to finish
wait

#Carry out multiqc across all samples
cd ${analysis_out_dir}
multiqc .

#Make the summary stats file
cd ${analysis_out_dir}/stats
echo \#Run ID\: ${RUNID} > summary_stats.txt
echo -n Sample_name, >> summary_stats.txt
awk 'FNR==3{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Forward_fastq, >> summary_stats.txt
awk 'FNR==4{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Reverse_fastq, >> summary_stats.txt
awk 'FNR==5{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Aligned_reads, >> summary_stats.txt
awk 'FNR==6{printf $0 >> "summary_stats.txt"}' *
sed 's/.$//' summary_stats.txt >> summary_stats.csv
rm summary_stats.txt
