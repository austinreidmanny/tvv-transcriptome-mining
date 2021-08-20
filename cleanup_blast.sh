#!/bin/bash

###################################################################################################
# ENSURE THE SCRIPT IS CALLED CORRECTLY
###################################################################################################
set -eo pipefail

usage() { echo -e "The objective of this script is to clean up the results of stand-alone-blast (sab). \n\n" \
                  "This program takes in a FASTA (given as output by sab) of reads that mapped \n" \
                  "to the previous reference and searches those again, this time against a \n" \
                  "reference nucleotide blast database. The output of this script is a \n" \
                  "tab-delimited text file with the following fields: \n" \
                  "(1) name of the read (2) title of best hit in NCBI (3) e-value score of match \n\n" \
                  "Usage: $0 -i input.fasta [options] \n\n" \
                  "Output: input.cleanup.results.txt \n\n" \
                  "Optional parameters: \n" \
                        "-e (evalue, e.g. 100, 1, or 1e-99; [default = 10]) \n" \
                        "-m (maximum amount of memory to use [in GB]; [default=16] ) \n" \
                        "-o (output directory for saving results; [default="./cleanup_results"]) \n" \
                        "-p (path to directory for database; " \
                            "[default='~/Documents/Research/sra/blastdbs/tvv1-5_no-spaces-in-headers_db'] ) \n" \
                        "-n (sets nucleotide program to blastn; [default= dc-megablast] ) \n" \
                        "-g (sets nucleotide program to megablast; [default= dc-megablast] ) \n\n" \
                      "Example of a complex run: \n" \
                      "$0 -i input.fasta -e 1e-3 -m 26 -g \n\n" \
                      "Exiting program. Please retry with corrected parameters..." >&2 && exit 1
        }
###################################################################################################

###################################################################################################
# Note about the RefSeq Viral nucleotide database used
###################################################################################################
# User will need to download this database from ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ ,
# download the two files 'viral.1.1.genomic.fna.gz' and 'viral.2.1.genomic.genomic.fna.gz',
# unzip and concatenate them together into a single viral refseq fasta file, and use that as
# input into the NCBI `makeblastdb` tool to make a blast nucleotide database.
# For posterity, this tool was tested with these steps performed at 10:30AM on 2019-05-03
###################################################################################################

###################################################################################################
# TAKE IN THE USER-PROVIDED PARAMETERS
###################################################################################################
# Store all parameters so I can save them to the log later
ALL_PARAMETERS=$@

# Read in the user-provided parameters
while getopts "i:e:m:o:p:ng*:" arg; do
        case ${arg} in
                i ) # Take in the input fasta
                  SEQDUMP=${OPTARG}
                        ;;
                e ) # set evalue
                  E_VALUE=${OPTARG}
                        ;;
                m ) # set max memory to use (in GB; if any letters are entered, discard those)
                  MEMORY_ENTERED=${OPTARG}
                  MEMORY_TO_USE=$(echo $MEMORY_ENTERED | sed 's/[^0-9]*//g')
                        ;;
                o ) # set the output directory
                  OUTPUT_DIRECTORY=${OPTARG}
                        ;;
                p ) #set path to NCBI nt database
                  PATH_TO_NT_DB=${OPTARG}
                        ;;
                n ) # switch to blastn
                  BLAST_TASK="blastn"
                        ;;
                g ) # switch to megablast
                  BLAST_TASK="megablast"
                        ;;
                * ) # Display help
                  usage
                        ;;
        esac
done
shift $(( OPTIND-1 ))
###################################################################################################

###################################################################################################
# PROCESS THOSE PARAMETERS
###################################################################################################
# If no input is provided, tell that to the user and exit
if [[ -z "${SEQDUMP}" ]] ; then
    echo -e "\nERROR: No input detected. \n"
    usage
fi

# If e-value wasn't provided by user, then set it to 10
if [[ -z ${E_VALUE} ]]; then
    E_VALUE="10"
fi

# If -n (blastn) or -g (megablast) flags were not given by user, default to dc-megablast
if [[ -z "${BLAST_TASK}" ]]; then
    BLAST_TASK="dc-megablast"
fi

# If path to NCBI nt database was not given, give default path
if [[ -z "${PATH_TO_NT_DB}" ]]; then
    PATH_TO_NT_DB="${HOME}/Documents/Research/sra/blastdbs/tvv1-5_no-spaces-in-headers_db"
fi
###################################################################################################

###################################################################################################
# Set up number of CPUs to use (use all available) and RAM (can be set by user, defaults to 16GB)
###################################################################################################
# CPUs (aka threads aka processors aka cores):
#   If 8 CPUs were used, BLAST fails & gives Segmentation Fault. Error stopped if <= 4 CPUs are used
#   Strategy: Use up to 4 CPUs, or maximum available if less than 4 CPUs available

# Use `nproc` if installed (Linux or MacOS with gnu-core-utils); otherwise use `systctl`
{   command -v nproc > /dev/null &&
    MAX_NUM_THREADS=$(nproc)
} ||
{   command -v sysctl > /dev/null &&
    MAX_NUM_THREADS=$(sysctl -n hw.ncpu)
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

# Set memory usage
if [[ -z ${MEMORY_TO_USE} ]]; then
    echo "No memory limit set by user. Defaulting to 16GB"
    MEMORY_TO_USE="16"
fi
###################################################################################################

###################################################################################################
# CREATE DIRECTORIES AND PREPARE NAMES FOR BLAST
###################################################################################################

# Create an output directory to run & store the BLAST files
## If user didn't provide one, just use a subdirectory in working directory: "./cleanup_results"
if [[ -z ${OUTPUT_DIRECTORY} ]]; then
    OUTPUT_DIRECTORY="./cleanup_results"
    mkdir -p ${OUTPUT_DIRECTORY}

else
    # If user provided a desired output directory: check to make sure output directory doesn't exist;
    # then create output directory; if error, just default to a results subdirectory within current dir
    if [[ ! -d ${OUTPUT_DIRECTORY} ]]; then
        mkdir -p ${OUTPUT_DIRECTORY} ||
        {
            echo "Cannot create user-provided output directory. Defaulting to ./cleanup_results/."
            OUTPUT_DIRECTORY="./cleanup_results/"
            mkdir ${OUTPUT_DIRECTORY}
        }
    fi
fi

###################################################################################################

###################################################################################################
# ENSURE THAT ALL REQUIRED SOFTWARE IS INSTALLED
###################################################################################################
command -v blastn > /dev/null ||
{ echo -e "This program requires 'blastn'. \n" \
        "Please install from NCBI website or NCBI github and retry. \n" \
        "Exiting..."; exit 2
    }
###################################################################################################

###################################################################################################
# PROCESS INPUT NAME SO IT'S EASIER TO HANDLE
###################################################################################################
# Create names for BLAST output file: first, truncate file path, leaving just the filename itself
SEQDUMP_FILE=${SEQDUMP##*/}

# Next, eliminate ".stand_alone_blast" and file extension, giving a cleaner name for blast
BLAST_NAME_SEQDUMP=$(basename ${SEQDUMP_FILE} ".stand_alone_blast.fasta")
###################################################################################################

###################################################################################################
# READ ALL INPUTS BACK TO USER
###################################################################################################
# Create log file
LOG_FILE=${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.log
touch ${LOG_FILE} ||
    {
    echo -e "ERROR: Cannot create log file. \n" \
            "Output directory may not be accessible for writing files. \n" \
            "Exiting... \n\n"
    exit 3
    }

# Copy initial launch command into the log
echo -e "\ncleanup_blast was launched with the following command: \n    '${BASH_SOURCE} ${ALL_PARAMETERS}' \n" \
        "at: `date`" | tee ${LOG_FILE}

# Read inputs back to the user and store them in the log
echo -e "\n" \
        "Input fasta (seqdump) file provided: ${SEQDUMP} \n" \
        "Database used: ${PATH_TO_NT_DB} \n" \
        "e-value: ${E_VALUE} \n" \
        "Blast program: nucleotide-blast > ${BLAST_TASK} \n" \
        "Number of processors to use: ${NUM_THREADS} \n" \
        "Memory limit: ${MEMORY_TO_USE}GB \n" \
        "Output directory: ${OUTPUT_DIRECTORY} \n\n" | tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# Run nucleotide blast
###################################################################################################
blastn \
-task ${BLAST_TASK} \
-db ${PATH_TO_NT_DB} \
-query ${SEQDUMP} \
-out ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.txt \
-evalue ${E_VALUE} \
-num_threads ${NUM_THREADS} \
-outfmt "6 qseqid evalue stitle qlen" \
-max_target_seqs 10000000 \
-max_hsps 1

# max_hsps only restricts to best level match PER VIRUS;
# if one read matches to multiple viruses, there will be multiple viruses reported

# Sort the reads, first by name then by evalue with lowest (strongest) on top; then, deduplicate
sort -k1,1 -k2g,2g ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.txt |
    sort -k1,1 -u > \
    ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.sorted.temp

# Create a final results file, with a header
echo -e "Query_name\t" \
        "e-value\t" \
        "Hit_name\t" \
        "Query_length" > \
        ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.txt

# Sort the hits so that the longest contigs will be on top
sort \
    -k4,4nr \
    -t $'\t' \
    ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.sorted.temp >> \
    ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.txt
###################################################################################################

###################################################################################################
# OUTPUT LOGS
###################################################################################################
# Make a token to indicate the job finished correctly
echo -e "Finished nucleotide BLAST (${BLAST_TASK}), using file '${BLAST_NAME_SEQDUMP}' to query against \n" \
        "blast database '${PATH_TO_NT_DB}' at: \n $(date) \n" | tee -a ${LOG_FILE}

# Print number of sequences in the input file:
echo -e "Number of sequences in original input file: "\
        "$(grep -c "^>" ${SEQDUMP})" | tee -a ${LOG_FILE}

# Print number of hits
echo -e "Number of hits in cleaned output hits list: " \
        "$(tail -n +2 ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.txt |
           wc -l |
           awk '{$1=$1};1' |
           cut -d " " -f 1) \
           \n" |
        tee -a ${LOG_FILE}

###################################################################################################

###################################################################################################
# Create a summary of hits, with names of the hits & the number of reads that mapped to each
###################################################################################################
echo -e "Counts\tNames" | tee ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.summary.txt

# The sed step removes all leading spaces from the counts column
tail -n +2 ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.txt |
    cut -f 3 |
    sort |
    uniq -c |
    sed -e 's/^[ ]*//g' |
    tr " " "\t" |
    tee -a ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.summary.txt
###################################################################################################

###################################################################################################
# Create a fasta file of all cleaned hits
###################################################################################################
seqtk \
    subseq \
    ${SEQDUMP} \
    <(tail -n +2 ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.txt | cut -f 1) > \
    ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.fasta
###################################################################################################

###################################################################################################
# Delete temporary files
###################################################################################################
#if using macOS, which has a pretty underpowered `find` utility, cannot use -name flag; so just use `rm`
find ${OUTPUT_DIRECTORY} -name "*temp" -print0 | xargs -0 rm || 
    rm "${OUTPUT_DIRECTORY}/*temp"
###################################################################################################

###################################################################################################
# Indicate to user that this analysis is complete
###################################################################################################
echo -e "\n$0 has finished analysis at: $(date)" | tee -a ${LOG_FILE}
###################################################################################################

###################################################################################################
# FIND READS THAT MAPPED INITIALLY, BUT FELL OUT DURING THIS CLEANUP STEP
###################################################################################################

#=================================================================================================#
#     I don't believe this analysis is critical at this time; 
#     just adds more files that could confuse the user;
#     disabled for now, may restore in the future
#=================================================================================================#


# Create temp files with reads in input vs. reads in output
#INPUT_READS="${OUTPUT_DIRECTORY}/input_reads.${BLAST_NAME_SEQDUMP}.temp"
#OUTPUT_READS="${OUTPUT_DIRECTORY}/output_reads.${BLAST_NAME_SEQDUMP}.temp"

#grep "^>" ${SEQDUMP} |
#    sed 's/>//g' |
#    cut -d " " -f 1 |
#    sort > \
#    ${INPUT_READS}

#cut -f 1 ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.results.txt |
#    cut -d " " -f 1 |
#    sort > \
#    ${OUTPUT_READS}

# Find the reads in one but not the other; extract those reads from the input
#seqtk \
#    subseq \
#    ${SEQDUMP} \
#    <(comm -3 ${INPUT_READS} ${OUTPUT_READS}) > \
#    ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.dropped-reads.fasta

# Print number of dropped out reads
#echo -e "Number of reads that dropped out of analysis during this cleaning step: " \
#        "`grep -c "^>" ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.dropped-reads.fasta | \
#          cut -d " " -f 1` \n" >> \
#        ${LOG_FILE}

# Save the names of those dropped out reads to the log
#echo -e "Names of any reads that dropped out during the analysis: \n" >> ${LOG_FILE}

#grep "^>" ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.dropped-reads.fasta |
#    tr -d ">" >> \
#    ${LOG_FILE}

# Remove the temporary files
#rm ${INPUT_READS}
#rm ${OUTPUT_READS}
#rm ${OUTPUT_DIRECTORY}/${BLAST_NAME_SEQDUMP}.cleanup_blast.dropped-reads.fasta

#find -name "*temp" -print0 | xargs -0 rm || 
#    rm "*temp" #if using macOS, which has a pretty underpowered `find` utility, cannot use -name flag; so just use `rm`
###################################################################################################
