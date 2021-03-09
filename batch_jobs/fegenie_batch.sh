#!/bin/bash
#SBATCH -A p30777               # Allocation
#SBATCH -p short                # Queue
#SBATCH -t 4:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 1                    # Number of cores
#SBATCH --mem=18G               # Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --mail-user=casar@u.northwestern.edu  # Designate email address for job communications
#SBATCH --mail-type=ALL     # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --job-name="fegenie batch"       # Name of job

#https://github.com/Arkadiy-Garber/FeGenie/wiki/Tutorial

#run FeGenie

#srun --account=p30777 --time=04:00:00 --partition=short --mem=18G --pty bash

module load python 
source activate fegenie

echo running FeGenie on genomes...

#to run on .faa files from binned genomes 
#/home/cpc7770/FeGenie/FeGenie.py -bin_dir DeMMO6_ORFs -bin_ext faa -out DeMMO6_fegenie_output --orfs

#loop over directories 
# for path in /projects/p30777/metagenome_data/genomes/*; do
#     [ -d "${path}" ] || continue # if not a directory, skip
#     dirname="$(basename "${path}")"
#     echo annotating $dirname genomes...
#     ~/FeGenie/FeGenie.py -bin_dir $path  -bin_ext faa -out /projects/p30777/metagenome_data/fegenie_genome_annotations/$dirname --orfs
# done

#annotate Momper 2017
~/FeGenie/FeGenie.py -bin_dir /projects/p30777/metagenome_data/Momper2017/Momper2017_MAGs  -bin_ext fa -out /projects/p30777/metagenome_data/fegenie_genome_annotations/Momper2017_Mags --meta

echo done!


echo running FeGenie on metagenomes...

# to run FeGenie on the metagenome assemblies 
# gene calls worked for all but D6 assemblies
# /home/cpc7770/FeGenie/FeGenie.py -bin_dir /projects/p30777/metagenome_data/contigs -bin_ext fa -out DeMMO_metagenomes_fegenie_output --meta


 #~/FeGenie/FeGenie.py -bin_dir /projects/p30777/metagenome_data/contigs/  -bin_ext faa -out /projects/p30777/metagenome_data/fegenie_metagenome_annotations/ --orfs

 ~/FeGenie/FeGenie.py -bin_dir /projects/p30777/metagenome_data/Momper2017/  -bin_ext fna -out /projects/p30777/metagenome_data/fegenie_metagenome_annotations/Momper2017 --meta


echo done!