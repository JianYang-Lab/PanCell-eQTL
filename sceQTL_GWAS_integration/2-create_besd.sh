#!/bin/bash
#SBATCH --mem=100G 
#SBATCH -c 1
#SBATCH -p intel-sc3,amd-ep2
#SBATCH --qos=normal
#SBATCH --job-name=sceQTL_besd
#SBATCH --output=sceQTL_besd_%j.out
#SBATCH --error=sceQTL_besd_%j.err

##############################################################################
# Create eQTL besd file from query for SMR
##############################################################################

tissue=$1
ancestry=$2
celltype=$3

SMR=/path/to/smr-1.3.1-linux-x86_64/smr
query_path=/path/to/${tissue}/sc-eQTL/results/query
besd_path=/path/to/${tissue}/sc-eQTL/results/besd
celltype_list=/path/to/${tissue}/script_eqtl/cell_type.list

dir_local=/data/$$
mkdir -p $dir_local

cd $dir_local

if [ -f ${query_path}/${tissue}_${celltype}_${ancestry}_sceQTL_query.txt ]; then	
		cp ${query_path}/${tissue}_${celltype}_${ancestry}_sceQTL_query.txt $dir_local
		query_file=${tissue}_${celltype}_${ancestry}_sceQTL_query.txt
		besd_file=${tissue}_${celltype}_${ancestry}_sceQTL
		if [ -f $query_file ]; then
#			if [ ! -f ${besd_path}/${besd_file}.besd ]; then
				$SMR --qfile ${query_file} --make-besd --out ${besd_file}
#			fi
		fi	
fi


cp *.besd ${besd_path}
cp *.epi ${besd_path}
cp *.esi ${besd_path}

cd /data
rm -rf $dir_local

echo "All jobs completed for ${ancestry}"
