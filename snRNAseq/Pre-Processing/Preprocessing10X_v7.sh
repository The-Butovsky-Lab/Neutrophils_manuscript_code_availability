#! /bin/bash

#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
# Author    : Emma Gerrits
# Date      : 08-nov-2018
# Datasets  : Fastq files from sequenced libraries generated with 10x Genomics Single Cell 3' RNA sequencing
# Purpose   : Automated Cell ranger script for the preproccesing of multiple samples in a dataset. Uses Cellranger 7.0.0 and the corresponding references. Does not run ABACUS.

# Updated on 27-05-2019
# added option to align to pre-mRNA and use 10x version 3

# Updated on 17-05-2020
# Switched to Cellranger v5

# Updated on 16-06-2022 by Mirjam Koster
# Switched to Cellranger 7.0.0
# Note: including introns is the standard now, so I changed the regions variable

# REQUIREMENTS
# add this script to the folder containing your fastq.gz files, the samplenames.txt file and the AdditionalCode_v7 folder

#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#

# ANSWER THE QUESTIONS TO SET THE CORRECT SETTINGS

printf "\n
        ---------------------------------------------------------------------------------------
        ---------------------------------------------------------------------------------------
        ---------------------------------------------------------------------------------------\n
        Pipeline for the preprocessing of 10x Genomics single-cell RNA-sequencing data by Emma Gerrits, updated to Cellranger 7 by Mirjam Koster. Abacus no longer included.\n
        Required files in the current directory:
        - Raw fastq files, unmodified extracted from the ERIBA server
        - samplenames.txt file
        - This Preprocessing10X.sh script
        \n
        Required scripts in the AdditionalCode_v7 directory:
        - OrganizeAndRename.py
        - merge_metrics.sh
        \n
        ---------------------------------------------------------------------------------------
        ---------------------------------------------------------------------------------------
        ---------------------------------------------------------------------------------------
        \n"

# Directories of the fastq file should start with the researchers name (created in Preprocessing.py)
# these names are extracted from the samplenames.txt file
printf "What are the common starting characters of your samples? (for example: 2022, MKO)\n"
read yourname

# Choose your reference genome
printf "\nFrom what species are your samples? (human/mouse) Default: human\n"
read species

if [ "$species" = "mouse" ]; then
    ref="mm10"
elif [ "$species" = "human" ] || [ -z "$species" ]; then
    ref="GRCh38"
fi 

# Choose whether to count exonic or exonic+intronic reads
# 10X recommends always including intronic reads. Emma used to exclude them for cells in CR v3.
printf "\nDo you want to include intronic reads in the alignment? (yes/no) Default: yes\n"
read mRNA

if [ "$mRNA" = "no" ]; then
    regions="--include-introns=false"
elif [ "$mRNA" = "yes" ] || [ -z "$mRNA" ]; then
    regions=""
fi

# Choose whether to perform secondary analysis
printf "\nDo you want to include secondary analysis (i.e. PCA, clustering, differential expression)? (yes/no) Default: yes\n"
read second

if [ "$second" = "no" ]; then
    secondary="--nosecondary"
elif [ "$second" = "yes" ] || [ -z "$second" ]; then
    secondary=""
fi

# organizing and reformating fastq files
printf "\nAre your fastq file names in the right format? (yes/no) Default: no
Make sure there are no underscores in the samplename!\n[samplename]_[Sxx]_L001_[Rx]_001.fastq.gz\n"
read format

dir=$(pwd)
TR="/data/bcn/Pipelines/Cellranger"

#-----------------------------------------------------------------------------------------------------------------------------------------------------------#

# RENAME THE FASTQ FILES IF NECESSARY
if [ "$format" = "no" ] || [ -z "$format" ]; then
   ./AdditionalCode_v7/OrganizeAndRename.py # run python organizing script
   mapfile -t IDs < <(ls | grep $yourname) # make array with the directories
   for j in ${!IDs[*]}
       do
   ( cd ${IDs[j]} || continue
        # Rename the fastq files
        printf '\nRenaming files in ./%s ...\n' "${IDs[j]}"
        # make array with the fastq file names
        i=0
        while read line
        do
            array[ $i ]="$line"        
            (( i++ ))
        done < <(ls) # change the name to the variable
        for i in ${!array[*]}
            do
            A="$(cut -d '_' -f1 <<< ${array[i]})"
            B="$(cut -d '_' -f2 <<< ${array[i]})"
            C="$(cut -d '_' -f3 <<< ${array[i]})"
            D="$(cut -d '_' -f4 <<< ${array[i]})"
            E="$(cut -d '_' -f5 <<< ${array[i]})"
            newname="${A}_${C}_L001_${D}_${E}"
            mv ./${array[i]} ./${newname}
        done
        printf "\nFiles renamed in the right format.\n"
)
done
fi

cd $dir

#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#

# make output directory
output="output_v7"
mkdir $output

printf "\nThe following samples were found:" 
ids=$(ls -d $yourname*)
set -f     
ids=(${ids/ / })
printf '\n%s\n' "${ids[@]}"
printf "\nDo you want to start Cellranger for these samples? (yes/no) Default: yes\n"
read decision

if [ "$decision" = "yes" ] || [ -z "$decision" ]; then
    #---------------------------------------------------------------------------------------------------------------------------------#
    # RUN CELLRANGER
    for i in ${!ids[*]}
        do
        # lines to run
        ID=${ids[$i]}
        FS=${ids[$i]}
        SM=( $(ls $dir/$ID | cut -d _ -f 1 | uniq | tr -d :) )
        # print the lines, will be added to a new .sh script
        printf "/data/bcn/Pipelines/Cellranger/cellranger-7.0.0/cellranger count --id=%s-out --fastqs=$dir/%s --sample=%s --transcriptome=$TR/refdata-gex-$ref-2020-A/ $secondary $regions
        %s\n" "$ID" "$FS" "$SM"
        # rename BAM file (for later use with fastqc)
        printf "mv $dir/$output/%s-out/outs/possorted_genome_bam.bam $dir/$output/%s-out/outs/%s.bam
        %s\n" "$ID" "$ID" "$ID"
    done > cellrangercode.sh # print the cellranger lines to a new .sh script

    sed -i "gunzip ./outs/*/*.gz" cellrangercode.sh # unzip the countfiles
    sed -i "1i cd ./$output" cellrangercode.sh # move to output folder
    sed -i "1i #!/bin/bash" cellrangercode.sh # add shebang
    printf "\n
    ---------------------------------------------------------------------------------------
    Running Cell ranger ...\n"
    chmod u+x ./cellrangercode.sh
    ./cellrangercode.sh # run the script
    rm ./cellrangercode.sh # remove the script
    printf "\n 
    ---------------------------------------------------------------------------------------
    Cell ranger is finished, output can be found in the output directory.\n"
    #---------------------------------------------------------------------------------------------------------------------------------#
    # FASTQC and MULTIQC
    printf "\n
    ---------------------------------------------------------------------------------------
    Running FastQC for each sample ...\n"
    ( cd ./$output || continue
    mkdir -p FASTQC
    find . -name "$yourname*.bam" | xargs -n 1 fastqc -o ./FASTQC
    printf "\n
    ---------------------------------------------------------------------------------------
    Running MultiQC ...\n"
    multiqc ./FASTQC -o ./FASTQC
    )
    #---------------------------------------------------------------------------------------------------------------------------------#
    # COMBINE METRIC FILES
    printf "\n
    ---------------------------------------------------------------------------------------
    Combining metric files ...\n"
    ./AdditionalCode_v7/merge_metrics.sh
    #---------------------------------------------------------------------------------------------------------------------------------#
    # RUN ABACUS
    #printf "\n
    #---------------------------------------------------------------------------------------
    #Running abacus ...\n"
    #./AdditionalCode_v7/abacus.sh # run abacus script
    #./AdditionalCode_v7/abacus_analysis.py # run the abacus analysis script
    #---------------------------------------------------------------------------------------------------------------------------------#     
    printf "\n
    ---------------------------------------------------------------------------------------
    Pipeline is finished. Output can be found in the output folder.\n"
elif [ "$decision" = "no" ]; then
    printf "\nPipeline stopped.\n"
fi

#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
