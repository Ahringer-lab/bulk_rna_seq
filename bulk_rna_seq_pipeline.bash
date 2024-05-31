#!/bin/bash

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
#      mergeID = If the file names have been merged differently the input can be changed here 'Myfile_<Add the flag here>_R1/R2_001.fastq.gz'
#      star_index = The location of the STAR index
#      kallisto_index = The location of the Kallisto index
# The scripts lacks logs, error handling and parallelisation (Beyond program specific multi-threading)
# Author Steve Walsh May 2024
##############################################################################################################

#Set the defaults
outdir=~/out
kallisto_index=~/references/built_genomes/kallisto/c.elegans_full_transcripts.idx
fastq_dir=~/data/
star_index=~/references/built_genomes/star/c.elegans.latest
CHROM_SIZES=/mnt/home3/ahringer/index_files/genomes/c_elegans.PRJNA13758.WS285.genomic.chrom.sizes
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
MERGEID=merged

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --fastqid <fastq suffix> --threads <number of threads> --input <input path> --id <Run ID>  --mergeID <merge ID> --star_index --kallisto_index"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l fastqid: -l threads: -l input: -l id: -l mergeID: -l star_index: -l kallisto_index: -- "$@") || exit_with_bad_args

#Get the inputs
eval set -- "$options"
while true; do
    case "$1" in
        --fastqid)
            shift
            base="$1"
            ;;
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


#Loops through the fastq names, make directories for each output, ${base} holds the sample name
echo "HERE!!!!"
echo "${base}${MERGEID}_R1_001.fastq.gz"
echo "${base}${MERGEID}_R2_001.fastq.gz"

analysis_out_dir=${outdir}/${RUNID}

mkdir ${analysis_out_dir}/${base}
mkdir ${analysis_out_dir}/${base}/fastq
mkdir ${analysis_out_dir}/${base}/trim_galore
trimmedfastq_dir=${analysis_out_dir}/${base}/trim_galore
mkdir ${analysis_out_dir}/${base}/star
mkdir ${analysis_out_dir}/${base}/fastq_screen
mkdir ${analysis_out_dir}/${base}/kallisto
cd ${analysis_out_dir}/${base}/fastq
cp $fastq_dir/${base}${MERGEID}_R*_001.fastq.gz .

#Carry out trimgalore (includes fastqc)
echo "trim_galore --fastqc ${analysis_out_dir}/${base}/fastq/${base}${MERGEID}_R1_001.fastq.gz ${analysis_out_dir}/${base}/fastq/${base}${MERGEID}_R2_001.fastq.gz \
-o ${analysis_out_dir}/${base}/trim_galore \
-j ${THREADS}"

#Carry out fastq screen
echo "fastq_screen ${trimmedfastq_dir}/*.fq.gz  \
--outdir ${analysis_out_dir}/${base}/fastq_screen \
--threads ${THREADS}"

#Carry out STAR alignment
#***N.B.*** Alignement carried out on un-trimmed reads due to the fussy nature of STAR with regard to it's input
echo "Carrying out STAR alignment"
echo "STAR --readFilesCommand zcat \
--runThreadN ${THREADS} \
--genomeDir $star_index \
--readFilesIn ${analysis_out_dir}/${base}/fastq/${base}${MERGEID}_R1_001.fastq.gz ${analysis_out_dir}/${base}/fastq/${base}${MERGEID}_R2_001.fastq.gz \
--outFileNamePrefix ${analysis_out_dir}/${base}/star/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMattrIHstart 0 \
--outWigType wiggle \
--twopassMode Basic"

#Converts wigs to bigwigs
echo "Converting wigs to bw"
echo "wigToBigWig ${analysis_out_dir}/${base}/star/Signal.UniqueMultiple.str1.out.wig ${CHROM_SIZES} ${analysis_out_dir}/${base}/star/${base}.UniqueMultiple.str1.bw"
#wigToBigWig ${analysis_out_dir}/${base}/star/Signal.UniqueMultiple.str2.out.wig ${CHROM_SIZES} ${analysis_out_dir}/${base}/star/${base}.UniqueMultiple.str2.bw
#wigToBigWig ${analysis_out_dir}/${base}/star/Signal.Unique.str1.out.wig ${CHROM_SIZES} ${analysis_out_dir}/${base}/star/${base}.Unique.str1.bw
#wigToBigWig ${analysis_out_dir}/${base}/star/Signal.Unique.str2.out.wig ${CHROM_SIZES} ${analysis_out_dir}/${base}/star/${base}.Unique.str2.bw

#Carry out Kallisto read quantification
echo "Carrying out quantification with Kallisto"
echo "kallisto quant -i ${kallisto_index} \
-b 100 \
-o ${analysis_out_dir}/${base}/kallisto \
-t 6 \
--rf-stranded \
${trimmedfastq_dir}/${base}${MERGEID}_R*_001_trimmed.fq.gz \
--threads=${THREADS}"

#Re-name Kallisto ouput
##cd ${analysis_out_dir}/${base}/kallisto
mv abundance.tsv ${base}abundance.tsv
#mv abundance.h5 ${base}abundance.h5
#mv run_info.json ${base}run_info.json
