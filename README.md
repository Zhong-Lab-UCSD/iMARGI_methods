# iMARGI Data Processing Methods

In this repository, we introduce the data processing methods (tools and usages) of iMARGI, from fastq to RNA-DNA interaction map.

- [iMARGI Data Processing Methods](#imargi-data-processing-methods)
    - [1. Tools and data dependencies](#1-tools-and-data-dependencies)
        - [1.1. Required tools](#11-required-tools)
        - [1.2. Reference data](#12-reference-data)
    - [2. Data processing](#2-data-processing)
        - [2.1. Input data and directory](#21-input-data-and-directory)
        - [2.2. Fastq file cleaning](#22-fastq-file-cleaning)
        - [2.3. Mapping](#23-mapping)
        - [2.4. Pair parsing, de-duplication and filtering](#24-pair-parsing-de-duplication-and-filtering)
        - [2.5. Multiple output formats](#25-multiple-output-formats)
            - [2.5.1. `.pairs` format](#251-pairs-format)
            - [2.5.2. BEDPE format](#252-bedpe-format)
    - [3. Further analysis tips](#3-further-analysis-tips)

## 1. Tools and data dependencies

To incorporate with standard pipelines of 4D Nucleome project (4DN), we use the same tools and reference data adopted
by 4DN DCIC (4D Nucleome Data Coordination and Integration Center) to process iMARGI data.

### 1.1. Required tools

Make sure that all the required tools are installed on your linux system with correct exec path configurations.

- [bwa](https://github.com/lh3/bwa) (version 0.7.17)
- [pairtools](https://github.com/mirnylab/pairtools) (version 0.2.0)
- [samtools](http://www.htslib.org/download/) (version 1.9)
- [seqtk](https://github.com/lh3/seqtk) (version 1.3)
- [pbgzip](https://github.com/nh13/pbgzip)
- [GNU parallel](https://www.gnu.org/software/parallel/)
- [ct_clean_pefastq.sh](./ct_clean_pefastq.sh): provided in this repository. Put it in your exec path and `chmod +x`,
   or specify directory when you use it.
- [pairs_to_bedpe.sh](./pairs_to_bedpe.sh): provided in this repository. Put it in your exec path and `chmod +x`, or specify directory
   when you use it.
- [GNU awk](https://www.gnu.org/software/gawk/manual/html_node/Quick-Installation.html): most of Linux distributions
  have GNU awk default installed.
- `bash`, `sort`, `gunzip`, `zcat`: Almost all of Linux distributions have these tools default installed.

### 1.2. Reference data

The reference genome and gene annotation GTF file used in 4DN are the same as ENCODE project. They can be downloaded
from [ENCODE project website](https://www.encodeproject.org/data-standards/reference-sequences/).

- Reference Genome: GRCh38_no_alt_analysis_set_GCA_000001405.15. [download link](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)
- Gene Annotations:  GRCh38 GENCODE V24 Primary (GTF format). [download link](https://www.encodeproject.org/files/ENCFF824ZKD/@@download/ENCFF824ZKD.gtf.gz)
- chromSize.txt. You can generate it or download from UCSC genome browser. We also provided it [here](./chromSize.txt).

## 2. Data processing

As the iMARGI data we generated are big high-throughput sequencing data (hundreds of millions read pairs), so a high
performance computing machine with Linux system is necessary, such as 8 CPU cores with 32 GB memory. Almost all of the
tools used here support parallel computing.

### 2.1. Input data and directory

The raw data generated by iMARGI experimental technique are paired-end sequencing read pairs in FASTQ format. That's
the only input data we needed. Besides, we also need reference genome and gene annotation file.
Here, we assume the input files and directory like this:

``` bash
|-- |
    |-- ref
    |   |-- GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    |   |-- gencode.v24.primary_assembly.annotation.gtf
    |
    |-- raw_data
    |   |-- HEK_iMARGI_R1.fastq.gz
    |   |-- HEK_iMARGI_R2.fastq.gz
    |   |-- chromSize.txt
    |   
    |-- clean_fastq
    |   |--
    |
    |-- mapping_output
    |   |--
    |
    |-- pairs_output
        |--
```

### 2.2. Fastq file cleaning

- input: raw paired FASTQ files `./raw_data/HEK_iMARGI_R1.fastq.gz` and `./raw_data/HEK_iMARGI_R2.fastq.gz`
- output: clean paired FASTQ files `./clean_fastq/clean_HEK_iMARGI_R1.fastq.gz` and `./clean_fastq/clean_HEK_iMARGI_R2.
  fastq.gz`

According to the design of iMARGI sequencing library construction, we need to clean the fastq files before mapping.
The first 2 bases (5' end) of RNA end (Read 1) are random nucleotides, which need to be removed before mapping.
The first 2 bases (5' end) of DNA end (Read 2) must be "CT" for correct ligation of RNA and DNA fragments. So we will
filter our the read pairs without "CT" as the first two bases of DNA end (Read 2).

We use the bash script `ct_clean_pefastq.sh` to clean the paired fastq files. We can use one command line to do that.
The `-t` and `-b` parameters are used for parallel computing, which means use 16 CPU cores and read in 4 million read
pairs each time in each CPU core. You need to adjust the parameters based on you machine configurations.

``` bash
bash ct_clean_pair.sh \
    -1 ./raw_data/HEK_iMARGI_R1.fastq.gz \
    -2 ./raw_data/HEK_iMARGI_R2.fastq.gz \
    -d true
    -o ./clean_fastq/ \
    -t 16 -b 4000000
```

Here is the usage and help info of the `ct_clean_pefastq.sh` script.

    Usage: ct_clean_pair.sh [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>] [-o <output_dir>] [-d <drop>]
                [-t <threads>] [-b <block_size>]

    Dependency: seqtk, zcat, awk, parallel

    This script will clean the DNA end reads (R2) of iMARGI sequencing Fastq data. According to the iMARGI design,
    DNA end reads (R2) must start with "CT". So it cleans the data by filtering out those R2 reads not starting with
    "CT". After filtering R2, it fixes the paired reads in R1. If "-d" was set as "true", it will drop all the non "CT"
    started R2 reads and paired R1 reads, which outputs two fastq files with prefix "clean_". If "-d" was "false", the
    filtered read pairs would also be outputed in a pair of fastq files with prefix "drop_". The default setting of "-d"
    is "false". The input fastq files must be gzip files, i.e., fastq.gz or fq.gz. The output files are also fastq.gz.

    -1 : R1 fastq.gz file, if there are multiple files, just separated with space, such as -1 l1_R1.fq l2_R1.fq
    -2 : R2 fastq.gz file, if there are multiple files, just separated with space, such as -1 l1_R2.fq l2_R2.fq
    -o : output directoy
    -d : Flag of dropping. Default is false. Option "true".
    -t : CPU threads for parallelized processing based on GNU parallel, Default 2
    -b : Fastq data block size (number of reads) for each thread. Default 2000000.
    -h : show usage help

### 2.3. Mapping

- input: clean paired FASTQ files `./clean_fastq/clean_HEK_iMARGI_R1.fastq.gz` and `./clean_fastq/clean_HEK_iMARGI_R2.
  fastq.gz`
- output: BAM file `./mapping_output/clean_HEK_iMARGI.bam`

We use `bwa mem` with `-SP5M` parameters to map cleaned read pairs to Human genome (hg38/GRch38). Before mapping, we
need to create the bwa index files from human reference genome sequence.

``` bash
# build bwa index
bwa index -p ./ref/bwa_index/bwa_index_GRCh38 ./ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fa

# mapping
bwa mem -t 16 -SP5M ./ref/bwa_index/bwa_index_GRCh38 \
    ./clean_fastq/clean_HEK_iMARGI_R1.fastq.gz \
    ./clean_fastq/clean_HEK_iMARGI_R2.fastq.gz |\
    samtools view -Shb - >./mapping_output/clean_HEK_iMARGI.bam
```

### 2.4. Pair parsing, de-duplication and filtering

- input: BAM file `./mapping_output/clean_HEK_iMARGI.bam`
- final output: `pairs` format (4DN) file `./pairs_output/final_HEK_iMARGI.pairs.gz`

We use pairtools in four steps. Here are some brief descriptions of the usages of pairtools in our work. If you want
to know more, please read the [Github repo](https://github.com/mirnylab/pairtools) and [documentation](https://pairtools.readthedocs.io/en/latest/) of pairtools. 

- `pairtools parse`: Parse the types of read pairs based on their mapping results. Output a `.pairsam` file. The 
  parameters `--add-columns cigar`, `--no-flip`, `--walks-policy 5any` and `--max-inter-align-gap 100` are important,
  don't change it unless you know what you want to do.
- `pairtools dedup`: Mark duplications in the `.pairsam` file.
- `pairtools split`: Split the `.pairsam` file into a SAM file and a `.pairs` file without sam1/sam2 fields. Use
  the `.pairs` file in further analysis as its size is much smaller than the `.pairsam` file.
- `pairtools select`: Filter out those read pairs which are marked as duplications, multiple mappings and 5' most end
  was not unique mapped. We used a long filter string designed for iMARGI, don't change it.

Here are the command lines we used.
``` bash
# parse
pairtools parse \
    -c ref/chromSize.txt \
    --assembly hg38 \
    --add-columns cigar \
    --no-flip \
    --output-stats stats_clean_HEK_iMARGI.pairs.txt \
    --walks-policy 5any \
    --max-inter-align-gap 100
    --nproc-in 16 \
    --nproc-out 8 \
    --output ./pairs_output/clean_HEK_iMARGI.pairsam.gz
    ./mapping_output/clean_HEK_iMARGI.bam

# mark duplications
pairtools dedup \
    --output-dups \
    --output-unmapped \
    --nproc-in 16 \
    --nproc-out 8 \
    --output ./pairs_output/markdup_clean_HEK_iMARGI.pairsam.gz \
    ./pairs_output/clean_HEK_iMARGI.pairsam.gz

# split pairsam to pairs and sam
pairtools split \
    --output-pairs ./pairs_output/markdup_clean_HEK_iMARGI.pairs.gz \
    --output-sam ./pairs_output/markdup_clean_HEK_iMARGI.bam \
    --nproc-in 8 \
    --nproc-out 16 \
    ./pairs_output/markdup_clean_HEK_iMARGI.pairsam.gz

# filter
pairtools select \
    '((strand1=="+" and regex_match(cigar1, "^\d+M.*")) or (strand1=="-" and regex_match(cigar1, ".*\d+M$"))) and ((strand2=="+" and regex_match(cigar2, "^\d+M.*")) or (strand2=="-" and regex_match(cigar2, ".*\d+M$"))) and regex_match(pair_type, "[UuR][UuR]")' \
    --chrom-subset ./ref/chromSize.txt \
    --nproc-in 8 \
    --nproc-out 16 \
    --output  ./pairs_output/final_HEK_iMARGI.pairs.gz \
    ./pairs_output/markdup_clean_HEK_iMARGI.pairs.gz
```  

### 2.5. Multiple output formats

The final output file is `final_HEK_iMARGI.pairs.gz` which is in `.pairs` format.

#### 2.5.1. `.pairs` format

`.pairs` file format is designed by 4DN.
You can read its [specification document](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md).

#### 2.5.2. BEDPE format

We also provide a script tool [`pairs_to_bedpe.sh`](./pairs_to_bedpe.sh) to convert the output the gzip compressed
`.pairs` format file to a gzip compressed BEDPE format file. **The pairs format file must include cigar columns**.
If you want to know more about BEDPE format, please read the
[docs of bedtools](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format).

``` bash
bash pairs_to_bedpe.sh \
    -i ./pairs_output/final_HEK_iMARGI.pairs.gz \
    -o ./pairs_output/final_HEK_iMARGI.bedpe.gz
```

## 3. Further analysis tips

According to iMARGI library construction design, the RNA end (Read 1) is reverse strand specific. It means that when
you annotate the RNA end with gene annotations, you need to reverse the strand, i.e., "+" -> "-" and "-" -> "+".
The DNA end (Read 2) is not strand specific.

Many tools can be used to annotate, such as bedtools and BEDOPS. But you need to write some scripts to annotate RNA
end and DNA end separately.

If you are using R, you can try to use Bioconductor packages
[GenomicRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html),
[GenomicInteractions](http://bioconductor.org/packages/release/bioc/html/GenomicInteractions.html), and
[InteractionSet](http://bioconductor.org/packages/release/bioc/html/InteractionSet.html), which will help you do more
analysis not only annotation.