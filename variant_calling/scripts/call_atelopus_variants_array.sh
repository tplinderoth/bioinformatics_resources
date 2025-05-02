#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J call_atelopus_variants
#! Account name for group, use SL2 for paying queue:
#! #SBATCH -A general
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=call_atelopus_variants_%A_%a.out
#! Errors filename:
#SBATCH --error=call_atelopus_variants_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=40
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=48:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=120000mb
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

REF='/mnt/research/Fitz_Lab/ref/frog/atelopus_ignescens/AteIgnes.fasta'
BAMS='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/map/atelopus_filtered_bam_list.txt'
SCAFFLIST="/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/scaffold_sets/AteIgnes_${SLURM_ARRAY_TASK_ID}.rf"
GROUPFILE='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/metadata/atelopus_wgs_qcpass_specimen_groups.tsv'
OUTBCF="/mnt/gs21/scratch/lindero1/atelopus/wgs/variants/atelopus_${SLURM_ARRAY_TASK_ID}.bcf.gz"

# Use call P = 0.0028, which is the average nucleotide diversity for contemporary Atelopus varius (pi = 0.0031)
# and captive Atelopus zeteki (pi = 0.0025) from Byrne etal (2020).

CMD="bcftools mpileup \
-f $REF \
-b $BAMS \
-R $SCAFFLIST \
-C 0 \
-d 1000000 \
-L 1000000 \
-q 20 \
-Q 20 \
--ns UNMAP,SECONDARY,QCFAIL,DUP \
-a FORMAT/AD,FORMAT/DP,FORMAT/QS,FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/SCR \
-p \
-m 2 \
-F 0.05 \
-O u \
--threads 8 \
| bcftools call \
--ploidy 2 \
-a PV4,GQ,GP \
-m \
-P 0.0028 \
-G $GROUPFILE \
-O u \
--threads 8 \
| bcftools +fill-tags \
-O u \
--threads 8 \
-- -t 'AF,ExcHet,NS' \
| bcftools norm \
-f $REF \
-O b \
-o $OUTBCF \
--threads 8"

printf "\n%s\n\n" "$CMD"

eval $CMD

wait

tabix -p bcf $OUTBCF

printf "\nfinished\n"
