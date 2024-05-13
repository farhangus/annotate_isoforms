#!/bin/bash

# Function to display usage information
display_aanotate_pipeline_usage() {
    echo "Options:"
    echo "  -h,--help                          Display this help message"
    echo "  -i,--input-file      <isoforms_list>     Input isoforms_list file"
    echo "  -b,--bed       <provided_bedfile>  provided_bedfile "
    echo "  -f,--fraction       Similarity fraction  default=0.9"
    echo "  -r,--reference <REF>               reference fasta file"
    # echo "  -m,--margin    base_pairs_margin   default=10"
    echo "  -o,--output    output_fle_path   output file path"
}

# Initialize variables
ISOFORMS_LIST=""
PROVIDED_BEDFILE=""
REFERENCE=""
OUTPUT_FILE="annotated_isoforms.bed"
MARGIN=0
TOTAL_NUMBER_OF_ISOFORMS=0
FOUND_ISOFORMS=0
FRACTION=.9
NUMBER_OF_ISOFORMS=""
ISOFORMS_NOT_FOUND=""
ISOFORMS_TRINITY=""

# Parse command-line options
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-file)
            ISOFORMS_LIST=$2
            shift 2
            ;;
        -b|--bed)
            PROVIDED_BEDFILE=$2
            shift 2
            ;;
        -r|--reference)
            REFERENCE=$2
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE=$2
            shift 2
            ;;
        -m|--margin)
            MARGIN=$2
            shift 2
            ;;
        -f|--fraction)
            FRACTION=$2
            shift 2
            ;;
        -h|--help)
            display_aanotate_pipeline_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            display_aanotate_pipeline_usage
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z $ISOFORMS_LIST || -z $PROVIDED_BEDFILE || -z $REFERENCE ]]
then
    echo "Missing required arguments."
    display_aanotate_pipeline_usage
    exit 1
fi
isoforms_locations=$(mktemp)
isoforms_locations='isoforms_locations.txt'
grep -v '^\s*$' "$ISOFORMS_LIST" >tmp_${isoforms_locations}
rm -f $isoforms_locations
# Iterate through each line in the ISOFORMS_LIST file
while IFS= read -r isoform; do
    # Check if the isoform exists in the provided BED file
    if grep -q "$isoform" "$PROVIDED_BEDFILE"; then
        # Get the lines from the BED file matching the isoform
        grep "$isoform" "$PROVIDED_BEDFILE" >> $isoforms_locations
    else 
        echo -e "_\t _\t _\t ${isoform}" >> $isoforms_locations
    fi
done < tmp_${isoforms_locations}

ISOFORMS_NOT_FOUND=$(grep -c '^_*_*_' $isoforms_locations)


echo -e "\e[31mTrinity Isoforms found  but not found in RefSeq \e[34m Fraction: ${FRACTION}\e[0m" > 'log.txt'
bedtools intersect -a "$isoforms_locations" -b mm39_genes.bed -f $FRACTION  -wa -wb  -v >> 'log.txt'
echo -e "\e[31mIsoforms intersection between Trinity and RefSeq \e[34m Fraction: ${FRACTION}\e[0m" >> 'log.txt'
bedtools intersect -a "$isoforms_locations" -b mm39_genes.bed -f $FRACTION  -wa -wb  >> 'log.txt'
awk '!seen[$1,$2,$3]++' log.txt  > "$OUTPUT_FILE"
echo -e "\e[31m***************************************************************************************\e[0m" > STATS.log
NUMBER_OF_ISOFORMS="$(grep -v '^\s*$' ${ISOFORMS_LIST} | wc -l  | grep -o '[0-9]*' )"
awk '!seen[$4]++' "${OUTPUT_FILE}" > "uniq_${OUTPUT_FILE}"
echo "NUMBER_OF_ISOFORMS: ${NUMBER_OF_ISOFORMS}" >> STATS.log
echo "ISOFORMS_NOT_FOUND: ${ISOFORMS_NOT_FOUND}" >> STATS.log
ISOFORMS_TRINITY=$(awk '/Trinity Isoforms found /{flag=1; next} /Isoforms intersection between Trinity and RefSeq/{ exit} {if(flag) {count++}} END {print count}' "${OUTPUT_FILE}")
ISOFORMS_TRINITY_UNIQ=$(awk '/Trinity Isoforms found /{flag=1; next} /Isoforms intersection between Trinity and RefSeq/{ exit} {if(flag) {count++}} END {print count}' "uniq_${OUTPUT_FILE}")
echo "NUMBER OF ISOFORMS ONLY REPORTED IN TRINITY: ${ISOFORMS_TRINITY}" >> STATS.log
total_isoform_found=$(wc -l "${OUTPUT_FILE}"| grep -o '[0-9]*')
total_isoform_found_uniq=$(wc -l "uniq_${OUTPUT_FILE}"| grep -o '[0-9]*')
echo "NUMBER OF SHARED ISOFORMS: $((total_isoform_found -ISOFORMS_TRINITY -2 ))" >> STATS.log
echo "NUMBER OF UNIQ ISOFORMS ONLY REPORTED IN TRINITY: ${ISOFORMS_TRINITY_UNIQ}" >> STATS.log
echo "NUMBER OF SHARED UNIQ ISOFORMS: $((total_isoform_found_uniq - ISOFORMS_TRINITY_UNIQ -2 ))" >> STATS.log

# echo "Tr_chr,Tr_start_,Tr_end,Tr_isoform,Tr_gene,Ref_chrom,Ref_chromStart,Ref_chromEnd,Ref_name,Ref_geneName,Ref_geneName2,Ref_geneType"
# while IFS= read -r line; do

#     bedfile_chromosome=$(awk '{print $1}' <<< "$line")
#     gene_start_position=$(awk '{print $2}' <<< "$line")   
#     gene_end_position=$(awk '{print $3}' <<< "$line")   
#     trinity_isoform_name=$(awk '{print $4}' <<< "$line")   
#     trinity_gene_name=$(sed 's/_[^_]*$//'<<<$trinity_isoform_name)


#     if [[ ! -z "$bedfile_chromosome" && ! -z "$gene_start_position" ]]
#     then
#     gene_start_position_start=$((gene_start_position - MARGIN))
#     gene_start_position_end=$((gene_start_position + MARGIN))
#     gene_end_position_start=$((gene_end_position - MARGIN))
#     gene_end_position_end=$((gene_end_position + MARGIN))
#     found_match=false
# #   Iterate over the range
#     for ((i = gene_start_position_start; i <= gene_start_position_end; i++)); do
#         if grep -q  -E -m 1 "^$bedfile_chromosome\s+$i\s+" "$REFERENCE"
#         then
#             line=$(grep -E -m 1 "^$bedfile_chromosome\s+$i\s+" "$REFERENCE")
#             new_string=$(sed 's/[[:space:]]\+/,/g' <<<"$line")
#             echo "$bedfile_chromosome","$gene_start_position","$gene_end_position" ,"$trinity_isoform_name" ,"$trinity_gene_name","$new_string"
#             found_match=true
#             ((TOTAL_NUMBER_OF_ISOFORMS++))
#             ((FOUND_ISOFORMS++))
#             break
#         else
#             line=""
#         fi
#     done
#     if ! $found_match
#     then
# for ((i = gene_end_position_start; i <= gene_end_position_end; i++)); do
#     if grep -E -q -m 1 "${bedfile_chromosome}\s.*\s${i}" "$REFERENCE"; then
#         line=$(grep -E -m 1 "${bedfile_chromosome}\s.*\s${i}" "$REFERENCE")
#         new_string=$(sed 's/[[:space:]]\+/,/g' <<<"$line")
#         echo "$bedfile_chromosome","$gene_start_position","$gene_end_position" ,"$trinity_isoform_name" ,"$trinity_gene_name","$new_string"
#         found_match=true
#         ((TOTAL_NUMBER_OF_ISOFORMS++))
#         ((FOUND_ISOFORMS++))
#         break
#     else
#         line=""
#     fi
# done
# fi

#     if ! $found_match 
#     then
#         ((TOTAL_NUMBER_OF_ISOFORMS++))
#         echo "$bedfile_chromosome","$gene_start_position","$gene_end_position" ,"$trinity_isoform_name" ,"$trinity_gene_name","$line"
#     fi
#     fi  
        
# done < "$isoforms_locations"
# # Print the first line in red
# echo -e "\e[31m***********************************************************************************\e[0m"

# # Print the second line in blue
# echo -e "\e[34mTOTAL_NUMBER_OF_ISOFORMS: $TOTAL_NUMBER_OF_ISOFORMS\e[0m"

# # Print the third line in blue
# echo -e "\e[34mFOUND_ISOFORMS: $FOUND_ISOFORMS\e[0m"
