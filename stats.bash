#!/bin/bash

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --input <input path> --id <Run ID>"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

INPUT=~/out
RUNID=RUNID="Stats_summary-$(date '+%Y-%m-%d-%R')"

#Set the possible input options
options=$(getopt -o '' -l input: -l id: -- "$@") || exit_with_bad_args

#Get the inputs
eval set -- "$options"
while true; do
    case "$1" in
        --input)
            shift
            INPUT="$1"
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

cd ${INPUT}/stats
echo \#Run ID\: ${RUNID} > summary_stats.txt
echo -n Sample_name, >> summary_stats.txt
awk 'FNR==3{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Forward_fastq, >> summary_stats.txt
awk 'FNR==4{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Reversqe_fastq, >> summary_stats.txt
awk 'FNR==5{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Aligned_reads, >> summary_stats.txt
awk 'FNR==6{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Aligned_percentage, >> summary_stats.txt
awk 'FNR==7{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q30_Aligned_reads, >> summary_stats.txt
awk 'FNR==8{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q30_Aligned_percentage, >> summary_stats.txt
awk 'FNR==9{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q10_Aligned_reads, >> summary_stats.txt
awk 'FNR==10{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q10_Aligned_percentage, >> summary_stats.txt
awk 'FNR==11{printf $0 >> "summary_stats.txt"}' *
sed 's/.$//' summary_stats.txt >> summary_stats.csv
rm summary_stats.txt