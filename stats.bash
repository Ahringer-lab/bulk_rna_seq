#!/bin/bash

analysis_out_dir = $1
RUNID = $2

cd ${analysis_out_dir}/stats
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