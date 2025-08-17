#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J trim_adapters_atelopus_aviti1
#! Account name for group, use SL2 for paying queue:
#SBATCH -A bioinformaticscore
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=trim_adapters_atelopus_aviti1_%A_%a.out
#! Errors filename:
#SBATCH --error=trim_adapters_atelopus_aviti1_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=8
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=24:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=24000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=2-22

#! This is the partition name.
#! #SBATCH -p cclake

#! mail alert at start, end and abortion of execution
#! emails will default to going to your email address
#! you can specify a different email address manually if needed.
#SBATCH --mail-type=FAIL

#! Don't put any #SBATCH directives below this line

#! Modify the environment seen by the application. For this example we need the default modules.
#! . /etc/profile.d/modules.sh                # This line enables the module command
#! module purge                               # Removes all modules still loaded
#! module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! export OMP_NUM_THREADS=1

#! The variable $SLURM_ARRAY_TASK_ID contains the array index for each job.
#! In this example, each job will be passed its index, so each output file will contain a different value
echo "This is job" $SLURM_ARRAY_TASK_ID

#! Command line that we want to run:
#! jobDir=Job_$SLURM_ARRAY_TASK_ID
#! mkdir $jobDir
#! cd $jobDir

workdir="$SLURM_SUBMIT_DIR" # The value of SLURM_SUBMIT_DIR sets workdir to the directory
cd $workdir

EXEC='cutadapt'
EXEC2='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/clean_data/scripts/batch2_aviti/rescueUnpaired.pl'
METADATA='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/metadata/atelopus_wgs_metadata.tsv'
FQDIR='/mnt/gs21/scratch/lindero1/atelopus/wgs/deduplicate_fastq'
FQPREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$METADATA" | cut -f2,4 | sed "s/[[:space:]]/_/")
BATCH='aviti1'
FQFWD="${FQDIR}/${FQPREFIX}_${BATCH}_R1.fastq.gz"
FQREV="${FQDIR}/${FQPREFIX}_${BATCH}_R2.fastq.gz"
OUTDIR='/mnt/gs21/scratch/lindero1/atelopus/wgs/trim_adapters'
OUTFWD="${OUTDIR}/${FQPREFIX}_${BATCH}_qcpass_R1.fastq.gz"
OUTREV="${OUTDIR}/${FQPREFIX}_${BATCH}_qcpass_R2.fastq.gz"
SHORTR1="${OUTDIR}/${FQPREFIX}_${BATCH}_short_R1.fastq.gz"
SHORTR2="${OUTDIR}/${FQPREFIX}_${BATCH}_short_R2.fastq.gz"
OUTUNPAIR="${OUTDIR}/${FQPREFIX}_${BATCH}_qcpass_U.fastq.gz"
LOGFILE="/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/clean_data/trim_adapters/stats/${FQPREFIX}_${BATCH}.cutadapt.json"

CMD="$EXEC -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --overlap 3 --trim-n --minimum-length 70 -o $OUTFWD -p $OUTREV \
--too-short-output $SHORTR1 --too-short-paired-output $SHORTR2 --json $LOGFILE --discard-casava $FQFWD $FQREV"

printf "\n%s\n\n" "$CMD"

eval $CMD

wait

CMD2="$EXEC2 --fq1 $SHORTR1 --fq2 $SHORTR2 --outname $OUTUNPAIR --minlen 70"
printf "\n%s\n\n" "$CMD2"
eval $CMD2

wait

rm $SHORTR1
rm $SHORTR2
