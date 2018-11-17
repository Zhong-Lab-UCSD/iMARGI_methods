#!/usr/bin/env bash
PROGNAME=$0

usage() {
    cat << EOF >&2
    Usage: $PROGNAME [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>] [-o <output_dir>] [-d <drop>] 
                [-t <threads>] [-b <block_size>]
    
    Dependency: seqtk, zcat, awk, parallel

    This script will clean the DNA end reads (R2) of iMARGI sequencing Fastq data. According to the iMARGI design,
    DNA end reads (R2) must start with "CT". So it cleans the data by filtering out those R2 reads not starting with "CT".
    After filtering R2, it fixes the paired reads in R1. If "-d" was set as "true", it will drop all the non "CT" started R2
    reads and paired R1 reads, which outputs two fastq files with prefix "clean_". If "-d" was "false", the filtered read
    pairs would also be outputed in a pair of fastq files with prefix "drop_". The default setting of "-d" is "false".
    The input fastq files must be gzip files, i.e., fastq.gz or fq.gz. The output files are also fastq.gz.

    -1 : R1 fastq.gz file, if there are multiple files, just separated with space, such as -1 l1_R1.fq l2_R1.fq
    -2 : R2 fastq.gz file, if there are multiple files, just separated with space, such as -1 l1_R2.fq l2_R2.fq
    -o : output directoy
    -d : Flag of dropping. Default is false. Option "true".
    -t : CPU threads for parallelized processing based on GNU parallel, Default 2
    -b : Fastq data block size (number of reads) for each thread. Default 2000000.
    -h : show usage help
EOF
    exit 1
}

while getopts :1:2:o:d:t:b:h opt; do
    case $opt in
        1) R1=${OPTARG};;
        2) R2=${OPTARG};;
        o) output_dir=${OPTARG};;
        d) dflag=${OPTARG};;
        t) threads=${OPTARG};;
        b) block=${OPTARG};;
        h) usage;;
    esac
done

[ ! -f "$R1" ] && echo "Error!! Fastq R1 not exist: "$R1 && usage
[ ! -f "$R2" ] && echo "Error!! Fastq R2 not exist: "$R2 && usage
[ ! -d "$output_dir" ] && echo "Error!! Output directory not exist: "$output_dir && usage
[  -z "$dflag" ] && echo "Use default setting '-d false'." && dflag=false
if [ "$dflag" != "false" ] && [ "$dflag" != "true" ]; then
    echo "Error!! Only true or false is acceptable for -d." && usage
fi
[  -z "$threads" ] && echo "Use default thread number 2'." && threads=2
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && usage 
fi
[  -z "$block" ] && echo "Use default block size of reads for each thread 2000000'." && block=2000000
if ! [[ "$block" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -b" && usage 
fi

clean_R1=$output_dir"/clean_"$(basename $R1)
clean_R2=$output_dir"/clean_"$(basename $R2)
drop_R1=$output_dir"/drop_"$(basename $R1)
drop_R2=$output_dir"/drop_"$(basename $R2)
[ -f "$clean_R1" ] && echo "Error!! Output clean fastq file exists: "$clean_R1 && usage
[ -f "$clean_R2" ] && echo "Error!! Output clean fastq file exists: "$clean_R2 && usage
[ -f "$drop_R1" ] && echo "Error!! Output drop fastq file exists: "$drop_R1 && usage
[ -f "$drop_R2" ] && echo "Error!! Output drop fastq file exists: "$drop_R2 && usage

tmpdir=$output_dir"/tmp_ct_clean_"$RANDOM""$RANDOM
mkdir $tmpdir

echo "Start processing reads $R1 and $R2 using $threads CPU cores with $block single-thread reads..."

seqtk mergepe $R1 $R2 | parallel -j $threads --pipe -L8 -N$block \
    "awk -v dflag=\"$dflag\" 'BEGIN{keep=0; 
                srand(systime()\"\"{%});
                randname=int(1000000*rand());
                OFS=\"\n\";
             }{
                idx=NR%8;
                if(idx==6){
                    if(\$0~/^[Cc][Tt]/){
                        keep=1;
                    }else{
                        keep=0;
                    }
                }
                if(idx==2 || idx==4){
                    tmp=substr(\$0, 3);
                    seqinfo[idx]=tmp;
                }else{
                    seqinfo[idx]=\$0;
                }
                if(idx==0){
                    if(keep==1){
                        print seqinfo[1],seqinfo[2],seqinfo[3],seqinfo[4] | \
                            \"gzip > $tmpdir/clean_\"randname\"_R1.fastq.gz\";
                        print seqinfo[5],seqinfo[6],seqinfo[7],seqinfo[0] | \
                            \"gzip > $tmpdir/clean_\"randname\"_R2.fastq.gz\";
                    }else{
                        if(dflag==\"false\"){
                            print seqinfo[1],seqinfo[2],seqinfo[3],seqinfo[4] | \
                                \"gzip > $tmpdir/drop_\"randname\"_R1.fastq.gz\";
                            print seqinfo[5],seqinfo[6],seqinfo[7],seqinfo[0] | \
                                \"gzip > $tmpdir/drop_\"randname\"_R2.fastq.gz\";
                        }
                    }            
                }
             }'"  
   
cat $tmpdir/clean_*_R1.fastq.gz >$clean_R1
cat $tmpdir/clean_*_R2.fastq.gz >$clean_R2
if [ "$dflag" == "false" ]; then   
    cat $tmpdir/drop_*_R1.fastq.gz >$drop_R1
    cat $tmpdir/drop_*_R2.fastq.gz >$drop_R2
fi

rm -r $tmpdir