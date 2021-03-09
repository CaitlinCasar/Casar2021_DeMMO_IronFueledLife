#!/bin/bash
#SBATCH -A p30777               # Allocation
#SBATCH -p long                # Queue
#SBATCH -t 10:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 8                    # Number of threads
#SBATCH --mem=48G               # Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --mail-user=casar@u.northwestern.edu  # Designate email address for job communications
#SBATCH --mail-type=ALL     # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --job-name="checkM"       # Name of job

module load python #python/anaconda is the default version, but 2.7 and 3 are available if needed
#module load idba/2016_12_longread
module load checkm/1.0.7

module load samtools/1.6

#get completeness, contamination, and putative lineage with CheckM

#run on metagenomes - note that the file extensions MUST be '.fasta', '.fa' is not acceptable.
#loop over directories 
for path in /projects/p30777/metagenome_data/genomes/*; do
    [ -d "${path}" ] || continue # if not a directory, skip
    dirname="$(basename "${path}")"
    echo calculating $dirname genome stats...
    checkm lineage_wf --tab_table -t 8 -x fasta /projects/p30777/metagenome_data/genomes/$dirname /projects/p30777/metagenome_data/checkm/$dirname -f /projects/p30777/metagenome_data/checkm/$dirname/checkm_stats.txt
done

