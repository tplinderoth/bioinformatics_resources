#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J split_msb_fastq
#! Account name for group, use SL2 for paying queue:
#SBATCH -A bioinformaticscore
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=split_msb_fastq_%A_%a.out
#! Errors filename:
#SBATCH --error=split_msb_fastq_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=12
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=36:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=36000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=1-253

#! This is the partition name.
#! #SBATCH -p cclake

#! mail alert at start, end and abortion of execution
#! emails will default to going to your email address
#! you can specify a different email address manually if needed.
#SBATCH --mail-type=FAIL

#! Don't put any #SBATCH directives below this line

#! Modify the environment seen by the application. For this example we need the default modules.
#! module purge                               # Removes all modules still loaded

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! export OMP_NUM_THREADS=1

#! The variable $SLURM_ARRAY_TASK_ID contains the array index for each job.
#! In this example, each job will be passed its index, so each output file will contain a different value
echo "This is job" $SLURM_ARRAY_TASK_ID

workdir="$SLURM_SUBMIT_DIR" # The value of SLURM_SUBMIT_DIR sets workdir to the directory
cd $workdir

EXEC='/mnt/research/Fitz_Lab/projects/msb/scripts/split_fastq_batches.pl'
SAMP_LIST="/mnt/research/Fitz_Lab/projects/msb/metadata/msb_libraries_unique.tsv"
SAMPID=$(tail -n+2 $SAMP_LIST | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f1)
TISSUE=$(tail -n+2 $SAMP_LIST | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f2)
FQFWD="/mnt/gs21/scratch/lindero1/msb/fastq/trim_adapters/${SAMPID}_qcpass_R1.fastq.gz"
FQREV="/mnt/gs21/scratch/lindero1/msb/fastq/trim_adapters/${SAMPID}_qcpass_R2.fastq.gz"
FQU="/mnt/gs21/scratch/lindero1/msb/fastq/trim_adapters/${SAMPID}_qcpass_U.fastq.gz"
LIB="CSM.${SAMPID}.${TISSUE}"
OUT="/mnt/gs21/scratch/lindero1/msb/fastq/batch_split/${SAMPID}_${LIB}_qcpass"

CMD="$EXEC --pe $FQFWD $FQREV --single $FQU --outprefix $OUT"

printf "\n%s\n\n" "$CMD"

eval $CMD
