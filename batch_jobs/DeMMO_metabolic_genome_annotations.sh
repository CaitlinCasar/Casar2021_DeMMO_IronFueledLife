#!/bin/bash
#SBATCH -A p30777               # Allocation
#SBATCH -p long                # Queue
#SBATCH -t 120:00:00             # Walltime/duration of the job
#SBATCH -N 1                    # Number of Nodes
#SBATCH -n 28                    # Number of threads
#SBATCH --mem=100G               # Memory in GB needed for a job. 
#SBATCH --mail-user=casar@u.northwestern.edu  # Designate email address for job communications
#SBATCH --mail-type=ALL     # Events options are job BEGIN, END, NONE, FAIL, REQUEUE
#SBATCH --job-name="DeMMO metabolic genome annotations"       # Name of job

module load anaconda3
source activate metabolic

#run on metagenomes - note that the file extensions MUST be '.fasta', '.fa' is not acceptable.
# perl /home/cpc7770/METABOLIC/METABOLIC-G.pl -in-gn /projects/p30777/metagenome_data/contigs -o /projects/p30777/metagenome_data/metabolic_metagenome_annotations


#copy files from RDSS to Quest
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/reads;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO1;  mget D1_S1_L003_R*' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/reads;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO2;  mget D2_S2_L003_R*' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/reads;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO3;  mget D3_S3_L003_R*' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/reads;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO4;  mget D4_S4_L003_R*' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/reads;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO5;  mget D5_S5_L003_R*' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/reads;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/DeMMO6;  mget D6_S6_L003_R*' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/reads;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/White_creek;  mget WC_S8_L003_R*' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/reads;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/Service_water;  mget SW_S7_L003_R*' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/genomes/WC;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/White_creek/White_creek_genomes;  mget *' -U "ads\cpc7770"
# smbclient '//resfiles.northwestern.edu/OSBURN_LAB' -c 'prompt;recursive;lcd  /projects/p30777/metagenome_data/genomes/SW;  cd DeMMO/DNA_Sequence_Data/Environmental_Metagenomes/Momper_Data/Service_water/Service_water_genomes;  mget *' -U "ads\cpc7770"

#run on binned genomes

#rename from .fa to .fasta
#cd /projects/p30777/metagenome_data/genomes
#find . -name "*.fa" -exec rename .fa .fasta {} +


#run on metagenomes - note that the file extensions MUST be '.fasta', '.fa' is not acceptable.
#loop over directories 
for path in /projects/p30777/metagenome_data/genomes/*; do
    [ -d "${path}" ] || continue # if not a directory, skip
    dirname="$(basename "${path}")"
    echo annotating $dirname genomes...
    perl ~/METABOLIC/METABOLIC-C.pl -in-gn $path  -r /projects/p30777/metagenome_data/reads/$dirname -o /projects/p30777/metagenome_data/metabolic_genome_annotations2/$dirname -p meta -t 28
done


echo done!