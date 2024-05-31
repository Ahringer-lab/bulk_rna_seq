#!/bin/bash
#SBATCH --job-name=Ahringer_RNA_seq  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 
#SBATCH --mem=8gb
#SBATCH --output=pipeline_%j.log # Standard output and error log 

#Set the defaults
outdir=~/out
kallisto_index=~/references/built_genomes/kallisto/c.elegans_full_transcripts.idx
fastq_dir=~/data/
star_index=~/references/built_genomes/star/c.elegans.latest
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

#Set up stats file
STATSFILE=${analysis_out_dir}/stats.csv
echo \#Run ID,${RUNID} >> $STATSFILE

MERGEID=merged
declare -A FILES

cd ~/data

for f in *fastq.gz; do                  # search the files with the suffix
    base=${f%${MERGEID}*}                        # remove after "_L001_" To make sample ID the hash key
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
          srun bash bulk_rna_seq_pipeline.bash --fastqid ${base} --threads ${THREADS} --input ${fastq_dir} --id ${RUNID} --mergeID ${MERGEID} --star_index ${star_index} --kallisto_index ${kallisto_index} &
done
