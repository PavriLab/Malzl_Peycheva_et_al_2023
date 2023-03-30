#!/bin/bash
#
#SBATCH --job-name=harp
#SBATCH --cpus-per-task=16
#SBATCH --mem=80GB
#SBATCH --qos=medium
#SBATCH --output=harp.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=daniel.malzl@imp.ac.at
module load bowtie2/2.2.9-foss-2017a
module load samtools/1.9-foss-2017a

for fq in `ls fastqs/*`;
do
	python3 utils/pyharp.py -f ${fq} -bti1 /groups/pavri/bioinfo/daniel/analysisMEF/repliseq/btindex/129_sv/129_sv -bti2 /groups/pavri/bioinfo/daniel/analysisMEF/repliseq/btindex/cast_ei/cast_ei -t 16 -o processed
done
