#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J annotate_atelopus_vcf
#! Account name for group, use SL2 for paying queue:
#SBATCH -A bioinformaticscore
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=annotate_atelopus_vcf_%A_%a.out
#! Errors filename:
#SBATCH --error=annotate_atelopus_vcf_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=28
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=48:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=84000mb
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

EXEC='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/scripts/annotate_vcf.pl'
BCF="/mnt/gs21/scratch/lindero1/atelopus/wgs/variants/atelopus_${SLURM_ARRAY_TASK_ID}.bcf.gz"
BEDMASK='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/mask/bam_stats/atelopus_all_qc_bamstats_mask_genome_fail.bed'
DEVFILE='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/mask/ngsparalog_analysis/dup_lr/ignescens_putative_biallelic_snps_duplr_genome.tsv'
GROUPFILE='/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/scripts/atelopus_group_file.txt'
OUTVCF="/mnt/gs21/scratch/lindero1/atelopus/wgs/variants/annotated/atelopus_annotated_${SLURM_ARRAY_TASK_ID}.vcf.gz"

CMD="bcftools view --no-version $BCF | $EXEC --dpbounds 50,314 --bed $BEDMASK --devlr_file $DEVFILE --max_devlr 37.19075 --overwrite --genorep $GROUPFILE | bgzip > $OUTVCF"

printf "\n%s\n\n" "$CMD"
eval $CMD
wait
tabix -p vcf -C $OUTVCF

printf "\nfinished\n"
