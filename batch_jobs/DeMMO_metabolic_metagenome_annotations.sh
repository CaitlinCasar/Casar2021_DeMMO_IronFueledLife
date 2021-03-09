#!/bin/bash
#SBATCH -A p30777               # Allocation
#SBATCH -p long                # Queue
#SBATCH -t 120:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 28                    # Number of threads
#SBATCH --mem=100G               # Memory in GB needed for a job. 
#SBATCH --mail-user=casar@u.northwestern.edu  # Designate email address for job communications
#SBATCH --mail-type=ALL     # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --job-name="DeMMO metabolic metagenome annotations"       # Name of job

module load anaconda3
source activate metabolic


#copy files from RDSS to Quest
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/contigs;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO1;  mget DeMMO1_simplified_contigs.fa' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/contigs;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO2;  mget DeMMO2_simplified_contigs.fa' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/contigs;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO3;  mget DeMMO3_simplified_contigs.fa' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/contigs;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO4;  mget DeMMO4_simplified_contigs.fa' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/contigs;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO5;  mget DeMMO5_simplified_contigs.fa' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/contigs;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO6;  mget DeMMO6_simplified_contigs.fa' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/contigs;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/White_creek;  mget WC_simplified_contigs.fa' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/contigs;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/Service_water;  mget SW_simplified_contigs.fa' -U "ads\cpc7770"

#rename from .fa to .fasta
#cd /projects/p30777/metagenome_data/genomes
#find . -name "*.fa" -exec rename .fa .fasta {} +

#run on metagenomes - note that the file extensions MUST be '.fasta', '.fa' is not acceptable.
perl /home/cpc7770/METABOLIC/METABOLIC-G.pl -in-gn /projects/p30777/metagenome_data/contigs -o /projects/p30777/metagenome_data/metabolic_metagenome_annotations -p meta -t 28
