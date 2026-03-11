#!/bin/bash

# map_reads.sh

VERSION='1.2.2'
NTHREAD=1
REF=''
SAMPLE_ID=''
SAMPLE_LIB=''
NAME_MOD=''
FQDIR=''
OUTDIR=''
OUTPREFIX=''
PLATFORM=''
DATADS=''
BWAOPT=''
SINGLE=0
EXCLUDE_UNMAPPED=0
ONLY_PAIRED=0

if [ $# -lt 8 ]
then
#	>&2 echo "map_reads.sh <reference fasta> <sample name, e.g. 'T_101'> <fastq directory> <output directory> <sequencing platform> <data description> [number threads]"
	>&2 printf "map_reads.sh <input>\n\n"
	>&2 printf "Required inputs:\n"
	>&2 printf "%s\n" "--sample_id   STRING   Sample name, e.g. FSJ_101"
	>&2 printf "%s\n" "--infq_dir    STRING   Directory that contains fastq files to map"
	>&2 printf "%s\n" "--ref         STRING   reference fasta"
	>&2 printf "%s\n\n" "--outdir      STRING   Output directory name"
	>&2 printf "Optional inputs:\n"
	>&2 printf "%s\n" "--sample_lib  STRING   Sequencing library name, e.g. TL09"
	>&2 printf "%s\n" "--name_mod    STRING   Extra text in fastq file names"
	>&2 printf "%s\n" "--platform    STRING   Sequencing platform, e.g. Illumina, Aviti"
	>&2 printf "%s\n" "--ds          STRING   Data description contained in single or double quotes, e.g. 'Florida Scrub-Jay whole genome sequencing data'"
	>&2 printf "%s\n" "--single               Evoke to indicate that data are single-end."
	>&2 printf "%s\n" "--outprefix   STRING   Output bam file will be outdir/outprefix_untrimmed.bam. Default bam name is sample_id along with sample_lib and name_mod information if provided."
	>&2 printf "%s [%i]\n" "--threads     INT      Number of threads" $NTHREAD
	>&2 printf "%s\n" "--bwaopt      STRING   bwa mem options other than -t and -R (entire string should be within single or double quotes)"
	>&2 printf "%s\n" "--exclude_unmapped     Exclude unmapped reads."
	>&2 printf "%s\n\n" "--only_paired          Only process and output reads that map in a proper pair."
	>&2 printf "Usage notes:\n"
	>&2 printf "\n%s\n" "Paired-end input fastq files must be gziped and named as <sample id>[_sample lib][_name mod][_batch number]_R<1|2>.fastq.gz"
	>&2 printf "%s\n" "Single-end input fastq files must be gziped and named as <sample id>[_sample lib][_name mod][_batch number]*.fastq.gz"
	>&2 printf "%s\n" "<> = required, [] = optional (notice that underscores separate optional strings)"
	>&2 printf "%s\n" "R1 = forward reads, R2 = reverse reads"
	>&2 printf "\n%s\n" "All batches for a sample contained in the input fastq directory will be mapped."
	>&2 printf "\n%s\n" "Currently this script only supports reads with headers that are in the Casava 1.8 format."
	>&2 printf "\nversion %s\n" $VERSION
	exit
fi

# parse inputs
while [[ $# -gt 1 ]]; do
	key="$1"
	case $key in
		--sample_id)
		  SAMPLE_ID="$2"
		  shift;;
		--infq_dir)
		  FQDIR="$2"
		  shift;;
		--ref)
		  REF="$2"
		  shift;;
		--platform)
		  PLATFORM="$2"
		  shift;;
		--ds)
		  DATADS="$2"
		  shift;;
		--outdir)
		  OUTDIR="$2"
		  shift;;
		--outprefix)
		  OUTPREFIX="$2"
		  shift;;
		--sample_lib)
		  SAMPLE_LIB="$2"
		  shift;;
		--name_mod)
		  NAME_MOD="$2"
		  shift;;
		--bwaopt)
		  BWAOPT="$2"
		  shift;;
		--threads)
		  NTHREAD="$2"
		  shift;;
		--single)
		  SINGLE=1
		  ;;
		--exclude_unmapped)
		  EXCLUDE_UNMAPPED=1
		  ;;
		--only_paired)
		  ONLY_PAIRED=1
		  ;;
		*)
		>&2 echo "Error: Unknown argument $1"
		exit 1;;
	esac
	shift
done

SAMPLE_NAME=''
if [ ! -z "$SAMPLE_ID" ]; then
	SAMPLE_NAME="$SAMPLE_ID"
else
	>&2 echo "Error: --sample_id argument is required"
	exit 1
fi

if [ ! -z "$SAMPLE_LIB" ]; then SAMPLE_NAME+="_${SAMPLE_LIB}"; fi

if [ ! -z "$NAME_MOD" ]; then SAMPLE_NAME+="_${NAME_MOD}"; fi

if [ ! -z "$FQDIR" ]; then
	if [ ! -d "$FQDIR" ]; then echo >&2 "Error: Input fastq file directory $FQDIR does not exist"; exit 1; fi
else
	>&2 echo "Error: --infq_dir argument is required"
	exit 1
fi

if [ ! -z "$REF" ]; then
	if [ ! -f "$REF" ]; then
		>&2 echo "Error: Unable to locate reference FASTA $REF";
		exit 1
	fi
else
	>&2 echo "Error: --sample_id argument is required"
	exit 1
fi

if [ ! -z "$OUTDIR" ]; then
	if [ ! -d "$OUTDIR" ]; then
		>&2 echo "Error: Output directory $OUTDIR does not exist"
		exit 1
	fi
else
	>&2 echo "Error: --outdir argument is required"
	exit 1
fi

if [ $NTHREAD -lt 1 ]; then >&2 echo "Error: Cannot request fewer than 1 threads"; exit 1; fi

subset_arg="" # optional arguments to samtools to only keep subsets of reads based on how they map

if [[ "$EXCLUDE_UNMAPPED" -eq 1 ]]; then
	subset_arg="-F 0x4"
fi

if [[ "$ONLY_PAIRED" -eq 1 ]]; then
	if [[ "$SINGLE" -eq 1 ]]; then
		>&2 echo "Error: specifying --only_paired and --single is incompatible"
		exit 1
	fi
	if [[ "$EXCLUDE_UNMAPPED" -eq 1 ]]; then >&2 echo "Warning: both --exclude_unmapped and --only_paired specified. Only properly paired reads will be processed."
	subset_arg="-f 0x3"
fi

if [[ "$FQDIR" == */ ]]; then FQDIR=$(echo "$FQDIR" | sed 's/\/$//'); fi

if [[ "$OUTDIR" == */ ]]; then OUTDIR=$(echo "$OUTDIR" | sed 's/\/$//'); fi
OUTFULL="${OUTDIR}/"
if [ ! -z "$OUTPREFIX" ]; then
	OUTFULL+="${OUTPREFIX}"
else
	OUTFULL+="${SAMPLE_NAME}"
fi

# collect all fastq files for the sample

printf "Processing all fastq files for $SAMPLE_NAME\n"

declare -a fqlist;
if [[ $SINGLE == 1 ]]; then
	fqlist=($(find $FQDIR -name "${SAMPLE_NAME}*.fastq.gz" | tr '\n' ' '))
else
	fqlist=($(find $FQDIR -name "${SAMPLE_NAME}*_R1.fastq.gz" | tr '\n' ' '))
fi
numlanes=${#fqlist[@]} # this is the number of fastq batches to process

# create an array of batch numbers
lane_n=()
for fq in ${fqlist[@]}; do
	lnn=$(echo "$fq" | perl -s -ne 'print $1 if $_ =~ /${samp}_(\d+)/' -- -samp="$SAMPLE_NAME")
	if [[ ! " ${lane_n[*]} " =~ [[:space:]]${lnn}[[:space:]] ]]; then lane_n+=("$lnn"); fi
done

# LOOP OVER LANES TO MAP EACH SET OF FASTQ FILES
declare -a readheader
lanecounter=1
for FWDFQ in ${fqlist[@]}
do
	REVFQ=''
	if [[ $SINGLE == 0 ]]; then REVFQ=$( echo "$FWDFQ" | sed 's/_R1\.fastq\.gz$/_R2\.fastq\.gz/'); fi

	prevrg='N'
	readgroup=''
	nmatch=0
	lstart=-3

	# need to check multiple headers to make sure sample barcodes converge (i.e. no sequencing error)
	while [ $nmatch -lt 10 ]
	do
		lstart=$((lstart + 4))
		lend=$((lstart + 1))
		readheader=$(zcat "$FWDFQ" | sed -n "${lstart}p;${lend}q" | sed 's/ [[:digit:]]//')
		readinfo=($(echo "$readheader" | perl -e 'chomp($read = <>); @arr = split(/:/,$read); print "@arr[1,2,3,$#arr]";'))
		barcode=$(sed 's/+/-/' <<< "${readinfo[3]}")
		readgroup="@RG\tID:${readinfo[0]}.${readinfo[1]}.${readinfo[2]}\tBC:${barcode}\tPU:${readinfo[1]}.${readinfo[2]}\tSM:${SAMPLE_ID}"

		if [ "$readgroup" = "$prevrg" ] && ! [[ "$barcode" =~ 'N' ]]; then ((nmatch++)); fi
		prevrg="$readgroup"
	done
	if [ ! -z "$SAMPLE_LIB" ]; then readgroup+="\tLB:${SAMPLE_LIB}"; fi
	if [ ! -z "$PLATFORM" ]; then readgroup+="\tPL:${PLATFORM}"; fi
	if [ ! -z "DATADS" ]; then readgroup+="\tDS:${DATADS}"; fi

	printf "\n--MAPPING READS (batch %i/%i)--\n" "$lanecounter" "$numlanes"
	#echo "$readgroup"; exit # debug

	# map reads

	OUTBAM="${OUTFULL}_raw"
	if [ $numlanes -gt 1 ]; then
		lanenum=${lane_n[$(($lanecounter-1))]}
		OUTBAM+="_${lanenum}.bam"
	else
		OUTBAM+=".bam"
	fi

	MAPCMD="bwa mem -t $NTHREAD -R '$readgroup'"
	if [ ! -z "$BWAOPT" ]; then MAPCMD+=" $BWAOPT"; fi
	if [[ $SINGLE == 1 ]]; then
		MAPCMD+=" $REF $FWDFQ"
	else
		MAPCMD+=" $REF $FWDFQ $REVFQ"
	fi
	MAPCMD+=" | samtools view -b"
	if [ ! -z "$subset_arg" ]; then MAPCMD+=" $subset_arg"; fi
	MAPCMD+=" > $OUTBAM"
	printf "\n%s\n\n" "$MAPCMD"
	eval $MAPCMD
	wait

	((lanecounter++))

done

printf "\n--MERGING SEPARATE BATCH BAMS--\n"

if [ $numlanes -gt 1 ]; then
	declare -a mergelist
	for i in ${lane_n[@]}
	do
		mergelist+=("${OUTFULL}_raw_${i}.bam")
	done
	MERGE_CMD="samtools merge -n -c -@ $NTHREAD -o ${OUTFULL}_raw.bam ${mergelist[@]}"
	printf "\n%s\n\n" "$MERGE_CMD"
	eval $MERGE_CMD
else
	printf "\nSingle batch, no merging necessary.\n"
fi

printf "\n--SORTING AND MARKING DUPLICATES--\n"

# The following pipe
# 1) Adds ms (mate score) tag for duplicate marking (bwa output should already be grouped by read name)
# 2) sorts the bam
# 3) marks duplicates

DIVTHREAD=$(perl -e '$thread = <>; print int($thread/3)' <<< "$NTHREAD")
OUTBAM2="${OUTFULL}_untrimmed.bam"

POSTMAP_CMD="samtools fixmate -m -@ $DIVTHREAD -u ${OUTFULL}_raw.bam - | samtools sort -u -@ $DIVTHREAD | samtools markdup -s -S -O BAM -@ $DIVTHREAD - ${OUTBAM2}"

printf "\n%s\n\n" "$POSTMAP_CMD"

eval $POSTMAP_CMD

printf "\n--INDEXING BAM--\n"

INDEX_CMD="samtools index -c -@ $NTHREAD ${OUTBAM2}"

printf "\n%s\n" "$INDEX_CMD"

eval $INDEX_CMD

printf "\n--FINISHED--\n"

exit
