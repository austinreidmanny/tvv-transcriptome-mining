#!/bin/bash

###################################################################################################
# ENSURE THE PIPELINE IS CALLED CORRECTLY
###################################################################################################

# Set up a usage statement in case this program is called incorrectly
usage() { echo -e "\nERROR: Missing input transcriptome(s) and/or input query and/or query type. \n\n" \
              "If using transcriptomes from the NCBI SRA database, \n" \
              "make sure to provide one (or more) SRA run numbers separated by commas \n" \
              "as well as a virus query (in fasta format), and indicate the query type \n" \
              "as either 'nucl' or 'prot' (do not include the quotes). \n\n" \
	      "If using local RNA-seq libraries as input, please indicate that as seen below... \n\n" \
	      "Proper usage (for using transcriptomes from NCBI SRA): \n" \
                "$0 -s SRR10001,SRR10002,SRR... -q query_sequence.fa -t nucl|prot \n\n" \
	      "Proper usage (for using local RNA-seq library with paired-end reads): \n" \
	            "$0 -1 reads_R1.fq -2 reads_R2.fq -q query_sequence.fa -t nucl|prot \n\n" \
	      "Proper usage (for using local RNA-seq library with unpaired or interleaved reads): \n" \
                "$0 -u reads.fq -q VIRUS_QUERY -t nucl|prot \n\n" \
          "Proper usage (for using  a fasta file for input; e.g.,  de novo assembled contigs or scaffolds): \n" \
                "$0 -f de_novo_assemblies.fasta -q query_sequence.fa -t nucl|prot \n\n" \
          "Optional parameters: \n" \
                "-e (evalue, e.g. 100, 1, or 1e-99; [default = 10]) \n" \
                "-c (custom job name, for naming files [HIGHLY RECOMMENDED TO PROVIDE THIS]; \n" \
                "    [default=none, guess from input filenames]) \n" \
                "-m (maximum amount of memory to use [in GB]; [default=16] ) \n" \
                "-o (output directory for saving results; [default="./results"]) \n" \
                "-p (path to directory for saving SRA files; [default='/tmp/'] ) \n" \
                "-d (sets nucloetide program to discontiguous-megablast; [default=blastn] ) \n" \
                "-g (sets nucleotide program to megablast; [default=blastn] ) \n\n" \
              "Example of a complex run: \n" \
              "$0 -1 sampleA1_reads_R1.fq -2 sampleA1_reads_R2.fq -q my_virus.fasta -t nucl -c sampleA1 -e 1e-3 -m 30 -d -o ./virus-blast-results \n\n" \
              "Exiting program. Please retry with corrected parameters..." >&2; exit 1; }

# Make sure the pipeline is invoked correctly, with project and sample names
while getopts "s:q:t:e:c:m:o:p:dg1:2:u:f:" arg; do
        case ${arg} in
                s ) # Take in the sample name(s)
                  set -f
                  IFS=","
                  ALL_SAMPLES=(${OPTARG}) # call this when you want every individual sample
                        ;;
                q ) # Take in the name of the virus query file (fasta format)
                  VIRUS_QUERY=${OPTARG}
                        ;;
                t ) # Take in the query type (nucleotide or protein)
                  QUERY_TYPE=${OPTARG}
                        ;;
                e ) # set evalue
                  E_VALUE=${OPTARG}
                        ;;
                c ) # set custom name for naming files
                  CUSTOM_NAME=${OPTARG}
                        ;;
                m ) # set max memory to use (in GB; if any letters are entered, discard those)
                  MEMORY_ENTERED=${OPTARG}
                  MEMORY_TO_USE=$(echo $MEMORY_ENTERED | sed 's/[^0-9]*//g')
                        ;;
                o ) # set the output directory
                  OUTPUT_DIRECTORY=${OPTARG}
                        ;;
                p ) #set path to SRA FILES
                  USER_PROVIDED_SRA_DIR=${OPTARG}
                        ;;
                d ) # switch to discontiguous_megablast
                  BLAST_TASK="dc-megablast"
                        ;;
                g ) # switch to megablast
                  BLAST_TASK="megablast"
                        ;;
                1 ) # if providing local transcriptome with paired-end reads, give path to forward reads fastq
                  FORWARD_READS=${OPTARG}
                  LOCAL_FILES="TRUE"
                        ;;
                2 ) # if providing local transcriptome with paired-end reads, give path to reverse reads fastq
                  REVERSE_READS=${OPTARG}
                  LOCAL_FILES="TRUE"
                        ;;
                u ) # if user wants to provide local transcriptome with unpaired reads, give path to reads fastq
                  UNPAIRED_READS=${OPTARG}
                  LOCAL_FILES="TRUE"
                        ;;
                f ) # if user wants to provide a fasta file for input (ie de novo assembled contigs or scaffolds), give path to fasta
                  FASTA_INPUT=${OPTARG}
                  LOCAL_FILES="TRUE"
                        ;;
                * ) # Display help
                  usage
                        ;;
        esac
done
shift $(( OPTIND-1 ))
###################################################################################################

###################################################################################################
# PROCESS THE USER PROVIDED PARAMETERS
###################################################################################################
# If the pipeline is not called correctly, tell that to the user and exit
if [[ -z "${VIRUS_QUERY}" ]] || [[ -z "${QUERY_TYPE}" ]] ; then
  usage
fi

# Check that SRA accessions or local files were provided
if [[ ! -z "${ALL_SAMPLES}" ]] && [[ ! -z "${LOCAL_FILES}" ]] ; then
  usage
fi

# But not both
if [[ -z "${ALL_SAMPLES}" ]] && [[ -z "${LOCAL_FILES}" ]] ; then
  usage
fi

# If local files are provided, check that only forward+reverse reads OR unpaired reads OR fasta input are given
if [[ ! -z "${LOCAL_FILES}" ]]; then
     if [[ ! -z "${FORWARD_READS}" ]] && [[ ! -z "${REVERSE_READS}" ]]; then
         if [[ ! -z "${UNPAIRED_READS}" ]]; then
             if [[ ! -z "${FASTA_INPUT}" ]]; then
                echo "If using local files, must provide only forward+reverse reads OR single file with unpaired reads OR fasta input";
                usage
            fi
         fi
     fi
fi

# If using SRA files, name the samples according to first & last
if [[ -z "${LOCAL_FILES}" ]]; then

  # Retrieve the name of last sample (using older but cross-platform compatible BASH notation)
  LAST_SAMPLE=${ALL_SAMPLES[${#ALL_SAMPLES[@]}-1]}

  # Create a variable that other parts of this pipeline can use mostly for naming
  SAMPLES="${ALL_SAMPLES[0]}-${LAST_SAMPLE}"

# If local files are provided, use the naming scheme of the forward reads files
else
    if [[ ! -z "${FORWARD_READS}" ]] && [[ ! -z "${REVERSE_READS}" ]]; then
        FORWARD_READS_FILE_WITH_NO_PATH=${FORWARD_READS##*/}
        FORWARD_READS_NO_PATH_NO_EXT=${FORWARD_READS_FILE_WITH_NO_PATH%.*}

        REVERSE_READS_FILE_WITH_NO_PATH=${REVERSE_READS##*/}
        REVERSE_READS_NO_PATH_NO_EXT=${REVERSE_READS_FILE_WITH_NO_PATH%.*}

        SAMPLES=${FORWARD_READS_NO_PATH_NO_EXT}
        ALL_SAMPLES=${SAMPLES}

    elif [[ ! -z "${UNPAIRED_READS}" ]]; then
        UNPAIRED_READS_FILE_WITH_NO_PATH=${UNPAIRED_READS##*/}
        UNPAIRED_READS_NO_PATH_NO_EXT=${UNPAIRED_READS_FILE_WITH_NO_PATH%.*}

        SAMPLES=${UNPAIRED_READS_NO_PATH_NO_EXT}
        ALL_SAMPLES=${SAMPLES}

    elif [[ ! -z "${FASTA_INPUT}" ]] ; then
        FASTA_INPUT_WITH_NO_PATH=${FASTA_INPUT##*/}
        FASTA_INPUT_NO_PATH_NO_EXT=${FASTA_INPUT_WITH_NO_PATH%.*}

        SAMPLES=${FASTA_INPUT_NO_PATH_NO_EXT}
        ALL_SAMPLES=${SAMPLES}
    fi
fi

# If user provides custom job name, use that for naming
if [[ ! -z "${CUSTOM_NAME}" ]]; then
    SAMPLES=${CUSTOM_NAME}
fi

# Reset global expansion
set +f

# Handle the query type provided by the user, using that to determine which type of blast to use
## Nucleotide query
if [[ ${QUERY_TYPE} == 'nucl' ]]; then
    BLAST_TYPE='blastn'

    # If -d (dc_megablast) or -g (megablast) flags were not given by user, default to megablast
    if [[ -z "${BLAST_TASK}" ]]; then
        BLAST_TASK='blastn'
    fi

## Protein query
elif [[ ${QUERY_TYPE} == 'prot' ]]; then
    BLAST_TYPE='tblastn'
    BLAST_TASK='tblastn'

## Other/Error
else
        echo "QUERY_TYPE is ${QUERY_TYPE}"
        echo "QUERY_TYPE must be 'nucl' or 'prot' (do not include quotes)"
        usage
fi

# If e-value wasn provided by user, then set it to 10
if [[ -z ${E_VALUE} ]]; then
    E_VALUE="10"
fi
###################################################################################################

###################################################################################################
# ENSURE THAT ALL REQUIRED SOFTWARE IS INSTALLED
###################################################################################################
command -v seqtk > /dev/null || \
{ echo -e "This program requires 'seqtk'. \n" \
        "Please install from github.com/lh3/seqtk and retry. \n" \
        "Exiting..."; exit 2
}

{ command -v makeblastdb > /dev/null && \
  command -v ${BLAST_TYPE} > /dev/null
} \
|| \
{ echo -e "This program requires 'makeblastdb'. \n" \
        "Please install from NCBI website or NCBI github and retry. \n" \
        "Exiting..."; exit 3
}
###################################################################################################

###################################################################################################
# Set up number of CPUs to use and RAM
###################################################################################################
# CPUs (aka threads aka processors aka cores):
#   If 8 CPUs were used, BLAST fails & gives Segmentation Fault. Error stopped if <= 4 CPUs are used
#   Strategy: Use up to 4 CPUs, or maximum available if less than 4 CPUs available

# Use `nproc` if installed (Linux or MacOS with gnu-core-utils); otherwise use `systctl`
{ \
    command -v nproc > /dev/null && \
    MAX_NUM_THREADS=`nproc` && \
    echo "Number of processors available (according to nproc): ${NUM_THREADS}"; \
} || \
{ \
    command -v sysctl > /dev/null && \
    MAX_NUM_THREADS=`sysctl -n hw.ncpu` && \
    echo "Number of processors available (according to sysctl): ${NUM_THREADS}";
}

# If maximum available threads is less than or equal to 4, use all threads; else use 4
if (( ${MAX_NUM_THREADS} > 4 )); then
    NUM_THREADS=4
elif (( ${MAX_NUM_THREADS} <= 4 )); then
    NUM_THREADS=${MAX_NUM_THREADS}
else
    echo "Error. Could not determine number of CPUs to use. Exiting..."
    exit 4
fi

# Set memory usage to 16GB if none given by user
if [[ -z ${MEMORY_TO_USE} ]]; then
    echo "No memory limit set by user. Defaulting to 16GB"
    MEMORY_TO_USE="16"
fi
###################################################################################################

###################################################################################################
# CREATE DIRECTORIES AND PREPARE NAMES FOR BLAST
###################################################################################################

# Create an output directory to run & store the BLAST files
## If user didn't provide one, just use a subdirectory in working directory: "./results"
if [[ -z ${OUTPUT_DIRECTORY} ]]; then
    OUTPUT_DIRECTORY="./results"
    mkdir -p ${OUTPUT_DIRECTORY}

else 
    # If user provided a desired output directory: check to make sure output directory doesn't exist; 
    # then create output directory; if error, just default to a results subdirectory within current dir
    if [[ ! -d ${OUTPUT_DIRECTORY} ]]; then
        mkdir -p ${OUTPUT_DIRECTORY} || \
        {
            echo "Cannot create user-provided output directory. Defaulting to ./results/."
            OUTPUT_DIRECTORY="./results/"
            mkdir ${OUTPUT_DIRECTORY}    
        }
    fi
fi

# Create names for BLAST output file:
## truncates file path, leaving just the filename itself
VIRUS_QUERY_FILE=${VIRUS_QUERY##*/}

## eliminates file extension, giving a cleaner name for blast
BLAST_NAME_VIRUS_QUERY=${VIRUS_QUERY_FILE%.*}

# Create log file
readonly LOG_FILE="${OUTPUT_DIRECTORY}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.stand_alone_blast.log"
touch ${LOG_FILE}

# Copy initial launch command into the log
echo -e "sab was launched with the following command: \n $0 $@ \n" > ${LOG_FILE}

# Set directory to save SRA files and build a blast database
  # If user provided an SRA_DIR on the command line, then set use that
  if [[ ! -z ${USER_PROVIDED_SRA_DIR} ]]; then
    SRA_DIR=${USER_PROVIDED_SRA_DIR}

  # If user did not provide a directory to save SRAs, but has the SRA_DIR in their environment (Nibert lab bioinformatics workstation does), use that
  elif [[ -z ${USER_PROVIDED_SRA_DIR} ]] && [[ ! -z ${SRA_DIR} ]]; then
    SRA_DIR=$SRA_DIR

  # If user did not provide a directory and there's no default, use /tmp/
  elif [[ -z ${USER_PROVIDED_SRA_DIR} ]] && [[ -z ${SRA_DIR} ]]; then
    SRA_DIR="/tmp/"

  # Any other case would be an error
  else
    echo "Error trying to set up SRA directory. Exiting" && exit 11

  fi

# Creat that directory
mkdir -p ${SRA_DIR}

# Make directory to save resulting BLASTDB
BLAST_DB_DIR="${SRA_DIR}/blastdbs"
mkdir -p ${BLAST_DB_DIR}

# Read inputs back to the user and store them in the log
if [[ -z "${LOCAL_FILES}" ]]; then
  echo -e "\n" \
          "SRA Accessions provided: ${ALL_SAMPLES[@]} \n" \
          "Virus query file provided: ${VIRUS_QUERY} \n" \
          "Molecule type (nucl or prot) of input query: ${QUERY_TYPE} \n" \
          "e-value: ${E_VALUE} \n" \
          "Blast program: ${BLAST_TYPE} > ${BLAST_TASK} \n" \
          "Number of processors to use: ${NUM_THREADS} \n" \
          "Memory limit: ${MEMORY_TO_USE}GB \n" \
          "Output directory: ${OUTPUT_DIRECTORY} \n\n" | tee -a ${LOG_FILE}
else
    if [[ ! -z ${FORWARD_READS} ]]; then
        echo -e "\n" \
          "User-provided files for sample: ${FORWARD_READS_NO_PATH_NO_EXT} ${REVERSE_READS_NO_PATH_NO_EXT} \n" \
          "Virus query file provided: ${VIRUS_QUERY} \n" \
          "Molecule type (nucl or prot) of input query: ${QUERY_TYPE} \n" \
          "e-value: ${E_VALUE} \n" \
          "Blast program: ${BLAST_TYPE} > ${BLAST_TASK} \n" \
          "Number of processors to use: ${NUM_THREADS} \n" \
          "Memory limit: ${MEMORY_TO_USE}GB \n" \
          "Output directory: ${OUTPUT_DIRECTORY} \n\n" | tee -a ${LOG_FILE}

    elif [[ ! -z ${UNPAIRED_READS} ]]; then
        echo -e "\n" \
          "User-provided files for sample: ${UNPAIRED_READS} \n" \
          "Virus query file provided: ${VIRUS_QUERY} \n" \
          "Molecule type (nucl or prot) of input query: ${QUERY_TYPE} \n" \
          "e-value: ${E_VALUE} \n" \
          "Blast program: ${BLAST_TYPE} > ${BLAST_TASK} \n" \
          "Number of processors to use: ${NUM_THREADS} \n" \
          "Memory limit: ${MEMORY_TO_USE}GB \n" \
          "Output directory: ${OUTPUT_DIRECTORY} \n\n" | tee -a ${LOG_FILE}
    
    elif [[ ! -z ${FASTA_INPUT} ]]; then
        echo -e "\n" \
          "User-provided files for sample: ${FASTA_INPUT} \n" \
          "Virus query file provided: ${VIRUS_QUERY} \n" \
          "Molecule type (nucl or prot) of input query: ${QUERY_TYPE} \n" \
          "e-value: ${E_VALUE} \n" \
          "Blast program: ${BLAST_TYPE} > ${BLAST_TASK} \n" \
          "Number of processors to use: ${NUM_THREADS} \n" \
          "Memory limit: ${MEMORY_TO_USE}GB \n" \
          "Output directory: ${OUTPUT_DIRECTORY} \n\n" | tee -a ${LOG_FILE}
    fi
fi

# If user provided a custom job name, print that to the user and store in log, as well
if [[ ! -z ${CUSTOM_NAME} ]]; then
    echo -e "Custom run name: ${CUSTOM_NAME} \n\n" |
    tee -a ${LOG_FILE}
fi
###################################################################################################

###################################################################################################
# BLAST DATABASE STEPS
###################################################################################################
# Devise name for blast database from every SRA given (eg: SRA1_SRA2_SRA3_db)
BLAST_DB_NAME="`echo -e ${ALL_SAMPLES[@]} | tr ' ' '_'`_db"

# Concatenate all SRA FASTQs together to make blastdb from & extract hits from
CONCATENATED_FASTQ=${BLAST_DB_DIR}/${BLAST_DB_NAME}.fq
CONCATENATED_FASTA=${BLAST_DB_DIR}/${BLAST_DB_NAME}.fasta

# Check to see if an intact blastdb & concatenated FASTA exist; if so, skip blastdb building step
blastdbcmd -db ${BLAST_DB_DIR}/${BLAST_DB_NAME} -info >> ${LOG_FILE} 2>&1
if [[ $? -eq 0 && -f ${CONCATENATED_FASTA} ]] ; then
    echo -e "BLAST DB for these SRA accessions already exists. \n" \
            "Continuing to BLAST alignment..." | tee ${LOG_FILE}

else

  ##################################################################################################
  # Download SRA files
  ##################################################################################################

  # If local provided files, just concatenate them together (if paired-end) or call the unpaired as the concat fastq
  if [[ ! -z ${LOCAL_FILES} ]]; then

      if [[ ! -z ${FORWARD_READS} ]] ; then
        cat ${FORWARD_READS} ${REVERSE_READS} > ${CONCATENATED_FASTQ}

      # If unpaired reads, call that fasta as the concat fastq
      elif [[ ! -z ${UNPAIRED_READS} ]] ; then
          CONCATENATED_FASTQ=${UNPAIRED_READS}

      # If fasta input, skip the fastq stage, but copy fasta file to the blast directory for later, to call reads from
      elif [[ ! -z ${FASTA_INPUT} ]] ; then
          cp ${FASTA_INPUT} ${CONCATENATED_FASTA}
      fi

  else
      # Add the download from SRA step to the timelog file
      echo -e "Downloading input FASTQs from the SRA at: `date` \n" | tee -a ${LOG_FILE}

      # Download fastq files from the SRA
      for SAMPLE in ${ALL_SAMPLES[@]};
        do \
           echo -e "Downloading FASTQ file(s) for ${SAMPLE} at: \n `date`" | tee -a ${LOG_FILE};
           fasterq-dump \
           --split-3 \
           -t "${SRA_DIR}/tmp/" \
           --progress --verbose \
           --skip-technical --rowid-as-name --print-read-nr \
           --threads=${NUM_THREADS} \
           --mem=${MEMORY_TO_USE}"GB" \
           --bufsize=1000MB \
           --curcache=1000MB \
           --outdir ${SRA_DIR} \
           ${SAMPLE} >> ${LOG_FILE};
        done

      echo -e "\n Finished downloading FASTQs from the SRA at: \n" \
              "`date` \n" | tee -a ${LOG_FILE}

      echo -e "These SRA FASTQ files are located at: ${SRA_DIR} \n\n" | tee -a ${LOG_FILE}
  fi
  ##################################################################################################

  ##################################################################################################
  # MAKE BLASTDB FROM TRANSCRIPTOME FILES
  ##################################################################################################
  # Put starting time of blastdb building in log file
  echo -e "Building BLAST database from SRA files at: \n `date`" | tee -a ${LOG_FILE}

  # If using SRA files, then concatenate the samples together
  if [[ -z ${LOCAL_FILES} ]]; then

      # If there is already a concat fastq, then it could be a partial file, so delete it
      if [[ -f ${CONCATENATED_FASTQ} ]]; then
        rm ${CONCATENATED_FASTQ}
      fi

      # For each SRA accession, add its reads to the master concatendated fastq file
      for SRA_FASTQ in ${ALL_SAMPLES[@]}
        do \
          cat ${SRA_DIR}/${SRA_FASTQ}*fastq >> ${CONCATENATED_FASTQ};
        done
  fi
  
  # If using fasta input, then no need to make a concatenated fasta
  if [[ -z ${FASTA_INPUT} ]]; then
    
     # Transform that FASTQ into FASTA so it's readable by BLAST
     seqtk seq -A ${CONCATENATED_FASTQ} > ${CONCATENATED_FASTA}
  fi
  
  # Create the blast databse
  makeblastdb \
  -dbtype nucl \
  -in  ${CONCATENATED_FASTA} \
  -title ${BLAST_DB_NAME} \
  -out ${BLAST_DB_DIR}/${BLAST_DB_NAME} \
  -logfile ${BLAST_DB_DIR}/${SAMPLES}_makeblastdb.log

  # Delete temporary concatenated fastq object (unless it's just the original unpaired reads file)
  if [[ -z "${UNPAIRED_READS}" ]] ; then 
      if [[ -f ${CONCATENATED_FASTQ} ]]; then
          rm ${CONCATENATED_FASTQ}
      fi
  fi

  # Put finishing time of blastdb building in log file
  echo -e "Completed buiding BLAST database at: `date` \n" \
          "Log file for this step available at: ${BLAST_DB_DIR}/${SAMPLES}_makeblastdb.log" |
          tee -a ${LOG_FILE}
fi
###################################################################################################

###################################################################################################
# RUN BLAST
###################################################################################################
# Print time that blast starts and write to log file
echo -e "Began running ${BLAST_TYPE} with samples ${SAMPLES} at:" | tee -a ${LOG_FILE}
date | tee -a ${LOG_FILE}

# Run blast
${BLAST_TYPE} \
-task ${BLAST_TASK} \
-db ${BLAST_DB_DIR}/${BLAST_DB_NAME} \
-query ${VIRUS_QUERY} \
-out ${OUTPUT_DIRECTORY}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.stand_alone_blast.results.txt \
-outfmt "6 qseqid sseqid evalue" \
-num_threads ${NUM_THREADS} \
-evalue ${E_VALUE} \
-max_target_seqs 100000000 \
-max_hsps 1

# Print time completed and write to log file as well
echo -e "Finished running ${BLAST_TYPE} at:" | tee -a ${LOG_FILE}
date | tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# Create FASTA sequence file of hits 
#   (uses seqtk subseq rather than blastdbcmd since blastbdcmd cannot parse sequence names longer than 50 nt;
#    this limitation could prevent accurate retrieval of sequences that rely on long names)
###################################################################################################
seqtk subseq  ${CONCATENATED_FASTA} \
    <(cut -f 2 ${OUTPUT_DIRECTORY}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.stand_alone_blast.results.txt | sort -u) > \
    ${OUTPUT_DIRECTORY}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.stand_alone_blast.fasta

# Print number of hits and save to log file
echo -e "\nNumber of hits, saved in fasta file:" | tee -a ${LOG_FILE}
grep -c "^>" ${OUTPUT_DIRECTORY}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.stand_alone_blast.fasta | tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# Create a summary of hits, with names of the hits & the number of reads that mapped to each 
###################################################################################################
echo -e "Counts\tName" > ${OUTPUT_DIRECTORY}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.stand_alone_blast.summary.txt

# The sed step removes all leading spaces from the counts column
cut -f 1 ${OUTPUT_DIRECTORY}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.stand_alone_blast.results.txt |
    sort |
    uniq -c |
    sed -e 's/^[ ]*//g' |
    tr " " "\t" >> \
    ${OUTPUT_DIRECTORY}/${BLAST_TASK}.${SAMPLES}.${BLAST_NAME_VIRUS_QUERY}.stand_alone_blast.summary.txt

###################################################################################################



###################################################################################################
# STAND-ALONE-BLAST (sab) FINISHED SUCCESSFULLY
###################################################################################################
echo -e "\n StandAloneBlast (sab) finished successfully at: \n" \
        "`date`" | tee -a ${LOG_FILE}
