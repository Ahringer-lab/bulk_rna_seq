#!/bin/bash
#
#TODO? SET UP SLURM OPTIONS FOR COMMAND LINE
#SBATCH --job-name=Ahringer_bulk_rna-seq_bash_pipeline
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem-per-cpu=4000

#############################################################################################################
############################## bulk rna-seq bash pipeline ###################################################
# This code will carryout a basic bulk rna-seq analysis pipeline for all fastq files in the input repository
# The pipeline is for use on a slurm hpc
# The pipeline assumes fastq file have been merged used using the bash script on github
# The plan going forward is to implement Ahringer pipelines in Nextflow so this will not be developed beyond a basic workflow.
# Options include:
#      threads = Will multi-thread any process to this number
#      input = Change the path of the input fastq files
#      id = Change the name of the output folder, the default is a datestamp
#      holdinput = keep input files where they are
#      mergeID = If the file names have been merged differently the input can be changed here 'Myfile_<Add the flag here>_R1/R2_001.fastq.gz'
# Author Steve Walsh May 2024
##############################################################################################################

#Set the defaults
outdir=~/out
genome_chr=~/references/built_genomes/star/c.elegans.latest
fastq_dir=~/data/
star_index=~/references/built_genomes/star/c.elegans.latest
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
HOLDINPUT=false
MERGEID=merged

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --threads <number of threads> --input <input path> --id <Run ID> --holdinput --mergeID <merge ID> "
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l threads: -l input: -l id: -l holdinput: -l mergeID -- "$@") || exit_with_bad_args

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
        --holdinput)
           HOLDINPUT="true"
           echo "Holding input fastq files in inputs folder"
           ;;
        --mergeID)
            shift
            MERGEID="$1"
            ;;
         --)
            shift
            break
            ;;
    esac
    shift
done

analysis_out_dir=${outdir}/${RUNID}
mkdir $analysis_out_dir

echo "$analysis_out_dir"

cd $fastq_dir

# Make array to store fastq name
declare -A FILES

#Get all fastq names from input folder
for f in *fastq.gz; do                  # search the files with the suffix
    base=${f%_${MERGEID}_*}                        # remove after "_L001_" To make sample ID the hash key
    if [[ $f == $base* ]] && [[ $f == *"R1"* ]]; then    # if the variable is the current sample ID and is forward
        FILES[$base]=$f                  # then store the filename
    elif [[ $f == $base* ]] && [[ $f == *"R2"* ]]; then # if the variable is the current sample and is reverse
        FILES[$base]+=" $f"
    fi
done

#Loops through the fastq names, make directories for their output and run fastqc
for base in "${!FILES[@]}"; do 
    echo "${base}_${MERGEID}_R1_001.fastq.gz"
    echo "${base}_${MERGEID}_R2_001.fastq.gz"

    mkdir ${analysis_out_dir}/${base} 
    mkdir ${analysis_out_dir}/${base}/fastq
    mkdir ${analysis_out_dir}/${base}/fastqc
    mkdir ${analysis_out_dir}/${base}/star
    mkdir ${analysis_out_dir}/${base}/kalisto
    cd ${analysis_out_dir}/${base}/fastq
    ln -s $fastq_dir/${base}_${MERGEID}_R*_001.fastq.gz .

    #Carry out fastqc
    #fastqc ${analysis_out_dir}/${base}/fastq/${base}_${MERGEID}_R1_001.fastq.gz ${analysis_out_dir}/${base}/fastq/${base}_${MERGEID}_R2_001.fastq.gz \
    #--outdir ${analysis_out_dir}/${base}/fastqc

    STAR --readFilesCommand zcat \
    --runThreadN ${THREADS} \
    --genomeDir $star_index \
    --readFilesIn ${base}_${MERGEID}_R1_001.fastq.gz ${base}_${MERGEID}_R2_001.fastq.gz \
    --outFileNamePrefix ${analysis_out_dir}/${base}/star/ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrIHstart 0 \
    --outWigType wiggle \
    --twopassMode Basic

   #kallisto quant -i $kallisto_idx \
   #-b 100 \
   #-o ${analysis_out_dir}/${base}/kalisto \
   #-t 6 \
   #--rf-stranded ${base}_${MERGEID}_R1_001.fastq.gz ${base}_${MERGEID}_R2_001.fastq.gz

done
