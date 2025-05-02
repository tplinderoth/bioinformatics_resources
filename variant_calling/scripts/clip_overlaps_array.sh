#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J clip_atelopus_overlap_reads
#! Account name for group, use SL2 for paying queue:
#SBATCH -A bioinformaticscore
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=clip_overlap_atelopus_reads_%A_%a.out
#! Errors filename:
#SBATCH --error=clip_overlap_atelopus_reads_%A_%a.err

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
#SBATCH --array=1-21

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

EXEC='bam clipOverlap'
IDLIST='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/metadata/atelopus_wgs_specimen_id.txt'
SAMPID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDLIST")
INBAM="/mnt/gs21/scratch/lindero1/atelopus/wgs/untrimmed_bam/${SAMPID}_untrimmed.bam"
OUTBAM="/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/map/filtered_bam/${SAMPID}_filtered.bam"
NTHREAD=10

CMD="$EXEC --in $INBAM --out $OUTBAM --stats --excludeFlags 0x20C --params --poolSize 1000000000 --noPhoneHome"
printf "\n%s\n" "$CMD"
eval $CMD

CMD2="samtools index -c -@ $NTHREAD $OUTBAM"
printf "\n%s\n" "$CMD2"
eval $CMD2
