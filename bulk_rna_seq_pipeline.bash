#!/bin/bash
#
#SBATCH --job-name=star_kallisto
#SBATCH --output=slurm_out/star_kallisto.%N.%j.out
#SBATCH --error=slurm_err/star_kallisto.%N.%j.err
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aft36@cam.ac.uk

#############################################################################################################
############################## bulk rna-seq bash pipeline ###################################################
# This code will carryout a basic bulk rna-seq analysis pipeline for all fastq files in the input repository
# The pipeline is for use on a slurm hpc
# Options include:
#      threads = Not currently used
#      path = Change the path of the input files
#      id = Change the name of the output folder, the default is a datestamp
#      holdinput = keep input files where they are
# Author Steve Walsh May 2024
##############################################################################################################

#Set the defaults
outdir=/mnt/home3/ahringer/sw2154/out
genome_chr=/mnt/home3/ahringer/sw2154/references/built_genomes/star/c.elegans.latest
fastq_dir=/mnt/home3/ahringer/sw2154/data/
star_index=/mnt/home3/ahringer/sw2154/references/built_genomes/star/c.elegans.latest
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
HOLDINPUT=false

#Set the possible input options
options=$(getopt -o '' -l threads: -l path: -l id: -l holdinput: -- "$@") || exit_with_bad_args

#Get the inputs
eval set -- "$options"
while true; do
    case "$1" in
        --threads)
            shift
            THREADS="$1"
            ;;
        --path)
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
    base=${f%_L001_*}                        # remove after "_L001_" To make sample ID the hash key
    if [[ $f == $base* ]] && [[ $f == *"R1"* ]]; then    # if the variable is the current sample ID and is forward
        FILES[$base]=$f                  # then store the filename
    elif [[ $f == $base* ]] && [[ $f == *"R2"* ]]; then # if the variable is the current sample and is reverse
        FILES[$base]+=" $f"
    fi
done

#THREADCOUNTER=0


#Loops through the fastq names, make directories for their output and run fastqc
for base in "${!FILES[@]}"; do 
    echo "${base}_L001_R1_001.fastq.gz"
    echo "${base}_L001_R2_001.fastq.gz"

    mkdir ${analysis_out_dir}/${base} 
    mkdir ${analysis_out_dir}/${base}/fastq

    cd ${analysis_out_dir}/${base}/fastq 

    if [ $HOLDINPUT == "false" ]; then
        mv $fastq_dir/${base}_L001_R*_001.fastq.gz .
    else
        cp $fastq_dir/${base}_L001_R*_001.fastq.gz .
    fi

    fastqc ${base}_L001_R1_001.fastq.gz ${base}_L001_R2_001.fastq.gz &

#    $THREADCOUNTER = $(( $THREADCOUNTER + 1 ))
#    if [ "$THREADCOUNTER" -ge "$THREADS" ]; then
#        echo "Maximum number of pipelines are running ($THREADS), waiting for them to finish"
#        wait
#        unset COUNTER
#        echo "Running the next group of pipelines now"
#        COUNTER=0
#    fi

done
