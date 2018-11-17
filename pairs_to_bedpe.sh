#!/usr/bin/env bash
PROGNAME=$0

usage() {
    cat << EOF >&2
    Usage: $PROGNAME [-i <input_pairs.gz_file>] [-o <output_bedpe.gz_file>]
    
    Dependency: zcat, awk, gzip, 

    This script will convert the pairs.gz format file to bedpe.gz format file. The pairs.gz format file must
    include cigar columns of both ends.

    -i : Input pairs.gz file
    -o : Output file name with directoy
    -h : show usage help
EOF
    exit 1
}

while getopts :i:o:h opt; do
    case $opt in
        i) input_file=${OPTARG};;
        o) output_file=${OPTARG};;
        h) usage;;
    esac
done

zcat $input_file | \
    awk 'BEGIN{FS="\t"; OFS="\t"}/^#/{next}
        {
            m1 = gensub(/[0-9]+S/, "", "g", $11); gsub(/M/, "", m1);
            m2 = gensub(/[0-9]+S/, "", "g", $12); gsub(/M/, "", m2); 
            if($6=="+"){
                start1=$3-1; end1=start1+m1;
            }else{
                end1=$3;start1=end1-m1;
            }; 
            if($7=="+"){
                start2=$5-1; end2=start2+m2;
            }else{
                end2=$5;start2=end2-m2;
            }; print $2, start1, end1, $4, start2, end2, $1, 1, $6, $7}' | gzip > $output_file