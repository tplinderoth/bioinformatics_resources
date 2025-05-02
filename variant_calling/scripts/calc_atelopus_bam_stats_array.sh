#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J atelopus_bam_stats
#! Account name for group, use SL2 for paying queue:
#SBATCH -A general
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=atelopus_bam_stats_%A_%a.out
#! Errors filename:
#SBATCH --error=atelopus_bam_stats_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=32
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=36:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=96000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=1-12

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

EXEC='/mnt/research/Fitz_Lab/software/bin/bamstats'
SCAFLIST="/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/scaffold_sets/AteIgnes_scaffolds_${SLURM_ARRAY_TASK_ID}.txt"
REF='/mnt/research/Fitz_Lab/ref/frog/atelopus_ignescens/AteIgnes.fasta'
BAM='/mnt/gs21/scratch/lindero1/atelopus/wgs/mask/atelopus_all_qc_merged.bam'
OUTFILE="/mnt/gs21/scratch/lindero1/atelopus/wgs/mask/atelopus_all_qc_${SLURM_ARRAY_TASK_ID}.bamstats"

printf "chr\tpos\tdepth\taverage_baseq\trms_baseq\tfraction_baseq0\taverage_mapq\trms_mapq\tfraction_mapq0\n" > $OUTFILE

while read -r SCAF
do
	CMD="$EXEC -A -d 1000000 -f $REF -q 0 -Q 0 -r $SCAF --ff UNMAP,SECONDARY,QCFAIL,DUP -s -aa $BAM | tail -n+2 >> $OUTFILE"
	printf "\n%s\n\n" "$CMD"
	eval $CMD
	wait
done < $SCAFLIST
wait
