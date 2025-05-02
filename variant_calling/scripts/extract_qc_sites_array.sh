#!/bin/bash
#! This line is a comment
#! Make sure you only have comments and #SBATCH directives between here and the end of the #SBATCH directives, or things will break
#! Name of the job:
#SBATCH -J extract_qc_sites
#! Account name for group, use SL2 for paying queue:
#SBATCH -A bioinformaticscore
#! Output filename:
#! %A means slurm job ID and %a means array index
#SBATCH --output=atelopus_qc_sites_%A_%a.out
#! Errors filename:
#SBATCH --error=atelopus_qc_sites_%A_%a.err

#! Number of nodes to be allocated for the job (for single core jobs always leave this at 1)
#SBATCH --nodes=1
#! Number of tasks. By default SLURM assumes 1 task per node and 1 CPU per task.
#SBATCH --ntasks=1
#! How many many cores will be allocated per task?
#SBATCH --cpus-per-task=12
#! Estimated runtime: hh:mm:ss (job is force-stopped after if exceeded):
#SBATCH --time=36:00:00
#! Estimated maximum memory needed (job is force-stopped if exceeded):
#SBATCH --mem=360000mb
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

VCF="/mnt/gs21/scratch/lindero1/atelopus/wgs/variants/annotated/atelopus_annotated_${SLURM_ARRAY_TASK_ID}.vcf.gz"
SITES="/mnt/research/Fitz_Lab/projects/atelopus/atelopus_wgs/variants/vcf/qc_sites/sp_comparisons/atelopus_qc_sites_${SLURM_ARRAY_TASK_ID}.pos"

bcftools view -f "PASS" -i 'N_PASS(FMT/DP[0-10] > 1) > 9 && N_PASS(FMT/DP[11] > 1) > 0 && N_PASS(FMT/DP[12-13] > 1) > 1 && N_PASS(FMT/DP[14] > 1) > 0 && N_PASS(FMT/DP[15] > 1) > 0 && N_PASS(FMT/DP[16] > 1) > 0 && N_PASS(FMT/DP[17-18] > 1) > 1 && N_PASS(FMT/DP[19] > 1) > 0 && (VT="." || VT="SNP")' -M 2 $VCF \
| bcftools query -f '%CHROM\t%POS\n' | uniq > $SITES

wait

angsd sites index $SITES

printf "\nfinished\n"
