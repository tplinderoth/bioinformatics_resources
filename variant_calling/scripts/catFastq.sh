#!/bin/bash

DIR=$1
DIR=$(echo $DIR | sed 's/\/$//')

OUTDIR=$2
OUTDIR=$(echo $OUTDIR | sed 's/\/$//')

FWDFQ=$(ls $DIR | grep "_1.fq.gz$")
REVFQ=$(ls $DIR | grep "_2.fq.gz$")

OUTPREFIX=$(echo $FWDFQ | cut -f1 -d ' ' | perl -e '$samp = <>; print $1 if $samp =~ /^([^_]+_\d+_[^_]+)_/')

# cat forward reads
FWDFQ=$( echo "$FWDFQ" | sed s:^:$DIR/:g )
printf "\nCAT FORWARD READS:\n%s\n" "$FWDFQ"
if (( $(echo "$FWDFQ" | wc -l) > 1 )); then
	FWDFQ=$(echo "$FWDFQ" | tr '\n' ' ')
fi
FWDFQ=$(echo "$FWDFQ" | sed 's/ $//')

zcat $FWDFQ | gzip > "${OUTDIR}/${OUTPREFIX}_1.fq.gz"

# cat reverse reads
REVFQ=$( echo "$REVFQ" | sed s:^:$DIR/:g )
printf "\nCAT REVERSE READS:\n%s\n" "$REVFQ"
if (( NREV=$(echo "$REVFQ" | wc -l) > 1 )); then
        REVFQ=$(echo "$REVFQ" | tr '\n' ' ')
fi
REVFQ=$( echo "$REVFQ" | sed 's/ $//')

zcat $REVFQ | gzip > "${OUTDIR}/${OUTPREFIX}_2.fq.gz"

exit
