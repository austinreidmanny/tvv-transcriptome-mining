#!/bin/bash

# ------------------------------------------------------------------------------------------------ #
# refineContigs.sh
#     Take TVV reads and map them to TVV contigs for assembly refinement
# ------------------------------------------------------------------------------------------------ #
echo -e "Welcome to 'refineContigs.sh' !"

# If one step fails, stop the script and exit
set -eo pipefail

# ------------------------------------------------------------------------------------------------ #
# Ensure the script is called correctly
# ------------------------------------------------------------------------------------------------ #
# Set up a usage statement in case this program is called incorrectly
usage() { echo -e "\nERROR: Missing sample name (needed for naming files) or input files. \n" \
                  "Proper usage for mapping tvv-reads to tvv-contigs: \n\n" \
                  "$0 -s mysample -r sample-tvv-mapped-reads.fq -c tvv-contigs.fasta \n\n" \
                  "Optional parameters: \n" \
                  "-o (output directory for saving trimmed files; [default = './' (current directory)]) \n" \
                  "-t (number of threads/CPUs); [default = 1] \n\n"
                  "Example of a complex run: \n" \
                  "$0 -s my_sample-r sample-tvv-mapped-reads.fq -c tvv-contigs.fasta -t 4 -o './' \n\n" \
                  "Exiting program. Please retry with corrected parameters..." >&2; exit 1;
        }

# Make sure the pipeline is invoked correctly, with project and sample names
while getopts "s:f:r:c:o:t:" arg; do
        case ${arg} in
                s ) # Take in the sample name for naming
                  sample=${OPTARG}
                        ;;
                f ) # path to forward reads fastq
                  forward_reads=${OPTARG}
                        ;;
                r ) # path to forward reads fastq
                  reverse_reads=${OPTARG}
                        ;;
                c ) # path to reverse reads fastq
                  contigs=${OPTARG}
                        ;;
                o ) # set the output directory
                  output_directory="${OPTARG%/}"
                        ;;
                t ) # set the number of threads/processors/CPUs
                  threads=${OPTARG}
                        ;;
                * ) # Display help
                  usage
                        ;;
        esac
done
shift $(( OPTIND-1 ))

# Check that required parameters are provided
if [[ -z "${sample}" ]] || [[ -z "${forward_reads}" ]] || [[ -z "${reverse_reads}" ]] || [[ -z "${contigs}" ]]; then
    usage
fi

# Set up an empty log file
cat /dev/null > ${sample}.refineContigs.log

## If user didn't provide number of threads, set to 1
if [[ -z ${threads} ]]; then
    threads="1"
fi


# Create an output directory to store the output files
## If user didn't provide one, just use a subdirectory in working directory
if [[ -z ${output_directory} ]]; then
    output_directory="./"
    mkdir -p ${output_directory}

else
    # If user provided a desired output directory: check to make sure output directory doesn't exist;
    # then create output directory; if error, just default to a results subdirectory within current dir
    if [[ ! -d ${output_directory} ]]; then
        mkdir -p ${output_directory} || \
        {
            echo "Cannot create user-provided output directory. Defaulting to current working directory './'" | \
            tee ${sample}.refineContigs.log

            output_directory="./"
            mkdir ${output_directory}
        }
    fi
fi
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
function refinement () {

    # --FUNCTION---------------------------------------------------------------------------------- #
    # Name:        refinement
    # Description: Map TVV species-specific reads to the TVV species-specific contigs in order to
    #              refine the contigs and resolve any SNPs or other errors made during contig assembly
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Ensure that the necessary software is installed
    command -v bwa > /dev/null && \
    command -v samtools > /dev/null && \
    command -v bcftools > /dev/null || {
        echo -e "ERROR: This script requires 'bwa' 'samtools' and 'bcftools' but could not found. \n" \
                "Please install these applications. \n" \
                "Exiting with error code 6..." >&2; exit 2
        }
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Log info
    echo "Began mapping to TVVs at:    $(date)" | \
    tee -a ${sample}.refineContigs.log
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Refinement
    # -------------------------------------------------------------------------------------------- #
    # Create a BWA index of the contigs
    bwa index \
    -p "${output_directory}/${sample}_contigs_index" \
    $contigs

    # Map the reads to the contigs
     bwa mem \
     -t $threads \
     "${output_directory}/${sample}_contigs_index" \
     $forward_reads $reverse_reads > \
     "${output_directory}/${sample}.reads_mapped_to_contigs.sam"

    # Get summary stats of the mapping
    samtools flagstat \
    --threads $threads \
    "${output_directory}/${sample}.reads_mapped_to_contigs.sam" > \
    "${output_directory}/${sample}.reads_mapped_to_contigs.stats"

    # Remove unmapped reads and sort output (save only contig-mapped reads
    samtools view --threads $threads -F 4 -bh "${output_directory}/${sample}.reads_mapped_to_contigs.sam" | \
    samtools sort --threads $threads - > "${output_directory}/${sample}.reads_mapped_to_contigs.sorted.bam"

    # Remove (very large) uncompressed sam file
    rm "${output_directory}/${sample}.reads_mapped_to_contigs.sam"

    # Add variable for mapped bam for easier reading
    mapped_bam="${output_directory}/${sample}.reads_mapped_to_contigs.sorted.bam"

    # Convert the TVV-aligned-reads (BAM) into a pileup (VCF)
    bcftools mpileup \
        --threads $threads \
        -d 1000000 \
        -f $contigs \
        $mapped_bam > \
        "${output_directory}/${sample}.pileup.vcf"

    # Call the variants (BCF file)
    bcftools call \
        --threads $threads \
        -m -Ob \
        -o "${output_directory}/${sample}.variants_called.bcf" \
        "${output_directory}/${sample}.pileup.vcf"

    # Index the calls.bcf file
    bcftools index --threads $threads "${output_directory}/${sample}.variants_called.bcf"

    # Combine the reference fasta and the called-variants into a consensus FASTA
    bcftools consensus \
        -f $contigs \
        "${output_directory}/${sample}.variants_called.bcf" > \
        "${output_directory}/${sample}.refined_contigs.fasta"

    # Compress the consensus fasta and move it to the main folder
    gzip "${output_directory}/${sample}.refined_contigs.fasta"
    # -------------------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------------------- #
    # Adapter trimming log info
    echo "Finished refinement at:    $(date)" | \
    tee -a "${sample}.refineContigs.log"
   # -------------------------------------------------------------------------------------------- #

}
# ------------------------------------------------------------------------------------------------ #

refinement
