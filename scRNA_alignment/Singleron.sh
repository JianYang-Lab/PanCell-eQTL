# ==============================================================================
# Shell Script for Singleron Single-Cell RNA-Seq Analysis
# ==============================================================================


# Step 1: Create a my.mapfile as below
#SRR14141993	/path/to/tissue/PRJNA719728	SRR14141993
#SRR14141994	/path/to/tissue/PRJNA719728	SRR14141994
#SRR14141995	/path/to/tissue/PRJNA719728	SRR14141995
#SRR14141996	/path/to/tissue/PRJNA719728	SRR14141996
#SRR14141997	/path/to/tissue/PRJNA719728	SRR14141997

# Note: use tab to separate, first column = prefix for fastq; second column = fastq file path;
# third column = sample name

# Step 2: activate conda environment and create bash script
tissue=$1
project=$2
sample=$3

cd /path/to/${tissue}/${project}
mkdir log
conda activate celescope

multi_rna \
	--mapfile my.mapfile \ # mapfile
	--genomeDir /path/to/ref/celescope/hs_ensembl_99 \ # Reference path
	--thread 8 \ # threads
	--mod shell \ # output file format
	--outdir ${sample} # output path

# Step 3: add sbatch parameters to bash scripts, then submit

for file in shell/*.sh; do
    {
        echo '#!/bin/bash' 
        echo "#SBATCH --qos=normal"
        echo "#SBATCH -p intel-sc3,amd-ep2"
        echo "#SBATCH --mem=80G"
        echo "#SBATCH -o log/%A-%a_out.txt"
        echo "#SBATCH -e log/%A-%a_error.txt"
        cat "$file"
    } > "${file}.tmp" && mv "${file}.tmp" "$file"
done

for file in shell/*.sh; do sbatch $file;echo $file;done
