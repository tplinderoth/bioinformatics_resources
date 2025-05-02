#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J map_atelopus_reads
#! Account name for group, use SL2 for paying queue:
#SBATCH -A general
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=atelopus_map_reads_%A_%a.out
#! Errors filename:
#SBATCH --error=atelopus_map_reads_%A_%a.err

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

EXEC='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/map/scripts/map_reads.sh'
REF='/mnt/research/Fitz_Lab/ref/frog/atelopus_ignescens/AteIgnes.fasta'
IDLIST='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/metadata/atelopus_wgs_specimen_id.txt'
FQDIR='/mnt/gs21/scratch/lindero1/atelopus/wgs/resplit'
OUTDIR='/mnt/gs21/scratch/lindero1/atelopus/wgs/untrimmed_bam'
SAMPID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDLIST")
NTHREAD=10

CMD="$EXEC $REF $SAMPID $FQDIR $OUTDIR 'Atelopus toad whole genome sequencing data produced by the Michigan State University RTSF Genomics Core, project JAY14502 KEJ_Atelopus_WG_RR_2023' $NTHREAD"
printf "\n%s\n" "$CMD"
eval $CMD
