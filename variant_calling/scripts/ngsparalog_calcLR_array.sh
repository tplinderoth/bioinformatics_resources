#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J atelopus_ngsparalog
#! Account name for group, use SL2 for paying queue:
#SBATCH -A bioinformaticscore
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=atelopus_ngsparalog_lr_%A_%a.out
#! Errors filename:
#SBATCH --error=atelopus_ngsparalog_lr_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=8
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=36:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=24000mb
#! Submit a job array with index values between 0 and 31
#! NOTE: This must be a range, not a single number
#SBATCH --array=1-1000

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

EXEC='/mnt/research/Fitz_Lab/software/ngsParalog/ngsParalog'
REF='/mnt/research/Fitz_Lab/ref/frog/atelopus_ignescens/AteIgnes.fasta'
BAM='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/map/ignescens_filtered_bam_list.txt'
POSFILE="/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/mask/ngsparalog_analysis/sites/pos_split2/atelopus_ignescens_putative_biallelic_snps_${SLURM_ARRAY_TASK_ID}.pos"
SCAFF=$(head -n 1 $POSFILE | cut -f 1)
START=$(head -n 1 $POSFILE | cut -f 2)
END=$(tail -n 1 $POSFILE | cut -f 2)
REGION="${SCAFF}:${START}-${END}"
OUTFILE="/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/mask/ngsparalog_analysis/dup_lr/ignescens_putative_biallelic_snps_duplr_${SLURM_ARRAY_TASK_ID}.tsv"

CMD="samtools mpileup -b $BAM -f $REF -r $REGION -l $POSFILE -A -a -q 0 -Q 0 --ff UNMAP,DUP,QCFAIL -d 110000 | $EXEC calcLR -infile - -outfile $OUTFILE -minQ 20 -minind 8 -mincov 1"

printf "\n%s\n\n" "$CMD"

eval $CMD
