### Install

To use these pipelines first clone the repository into your home directory on the cluster.

Then create two directories in your home directory:

* Data

* Out

Finally clone the conda environment as below:

conda env create -f environment.yml

### Running the pipeline

Drop all files into the data folder, by default the pipeline will take files with the following format:

<my_fastq>_merged_R1/2_001.fastq.gz

This will usually be files that have been merged with the file merger script, if you want to use a differnt file name format the mergeID flag should be used (see below)

From within the directory containing the scripts first edit the sample_sheet.csv file, you should list the files you want to analyse in the first column without _merged_R1/2_001.fastq.gz in the file name and 
in the second column list any updates you want to make to the output file names. 

Once the sample_sheet.csv file is ready fun the pipeline with the following command:

sbatch sbatch_main_RNAseq_pipeline.bash

*Note: You may have to change the Slurm cluster settings at the top of sbatch_pipeline.bash depending on how many files you are analysing and how busy the cluster is.

The following flags can be used to change the behaviour of the pipeline:

--threads : To run the individual programs within the pipeline by this many threads (The input file are run in tandem by default)

--input : Change the input folder from the default

--id : Change the ouput folder name from the default of a date/time stamp

--mergeID : Change the default input from <my_fastq>_merged_R1/2_001.fastq.gz. It will change the _merged_ part in the middle.

--star_index : Change the dedault location of the star index from /mnt/home3/ahringer/index_files/built_indexes/star/c.elegans.full.2024/

--kallisto_index : Change the dedault location of the star index from /mnt/home3/ahringer/index_files/built_indexes/kallisto/c.elegans.full.april2024/c_elegans.PRJNA13758.WS285.canonical_geneset.g.idx

### Output

All files will be ouput to the out directory in your home directory, there will be a folder for each sample that contains the following folders:

* fastq  
* fastq_screen  
* star  
* trim_galore

The fastq folder contains a copy of the fastq files used (as a point of reference to check back) and the other folders contain the output of the respective programs.

### Additional information

There are several additional 'mini' pipelines in the repository including:

* feature_counts.bash : This is simply to use feature counts to count against bam files, it is initiated in a similar way to the main pipeline but with sample_sheet_fc.csv instead of sample_sheet.csv
* htseq.bash : This is simply to use htseq to count against bam files, it also uses sample_sheet_fc.csv and runs like feature_counts.bash
* sbatch_pipeline.bash : This is a variant of the main pipeline that simply collects all the pairs of fastq files from the input directory with the defaul input file path <my_fastq>_merged_R1/2_001.fastq.gz. There is no option to update ouput file names with this pipeline.
