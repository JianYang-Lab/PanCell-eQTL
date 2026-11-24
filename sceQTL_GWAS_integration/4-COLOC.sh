#!/bin/bash
#SBATCH -p intel-sc3,amd-ep2
#SBATCH -c 8
#SBATCH --mem=60G
#SBATCH --qos=huge
#SBATCH --job-name=sceQTL_COLOC
#SBATCH --output=sceQTL_COLOC_%j.out
#SBATCH --error=sceQTL_COLOC_%j.err

set -euo pipefail
module load R/4.3.1

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 {default|diff|fdr} GWAS_ROW GWAS_ANCESTRY TISSUE EQTL_ANCESTRY CELL_TYPE" >&2
    exit 1
fi
mode=$1
GWAS_id=$2
GWAS_ancestry=$3
tissue=$4
ancestry=$5
celltype=$6

case "$mode" in
    default) result_subdir=raw ;;
    diff)    result_subdir=raw_diff ;;
    fdr)     result_subdir=raw_fdr ;;
    *) echo "Unknown COLOC mode: $mode" >&2; exit 1 ;;
esac

# sbatch can execute a spooled copy of this file. Find both R files before changing directory.
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
submit_dir=${SLURM_SUBMIT_DIR:-$PWD}
if [[ -n ${COLOC_SCRIPT_DIR:-} ]]; then
    candidates=("$COLOC_SCRIPT_DIR")
else
    candidates=("$script_dir" "$submit_dir" "$submit_dir/sceQTL_GWAS_integration")
fi
r_dir=""
for candidate in "${candidates[@]}"; do
    if [[ -r "$candidate/COLOC_perm.R" && -r "$candidate/CommonFunc.r" ]]; then
        r_dir=$(cd -- "$candidate" && pwd -P)
        break
    fi
done
if [[ -z "$r_dir" ]]; then
    echo "Set COLOC_SCRIPT_DIR to the folder containing COLOC_perm.R and CommonFunc.r." >&2
    exit 1
fi

SMR=/path/to/smr-1.3.1-linux-x86_64/smr
GWAS_list=/path/to/GWAS_${GWAS_ancestry}_list.txt
GWAS_trait=$(awk -F '\t' -v n="$GWAS_id" 'NR == n {print $1; exit}' "$GWAS_list")
GWAS_path=$(awk -F '\t' -v n="$GWAS_id" 'NR == n {print $4; exit}' "$GWAS_list")
GWAS_N=$(awk -F '\t' -v n="$GWAS_id" 'NR == n {print $5; exit}' "$GWAS_list")
if [[ -z "$GWAS_trait" || -z "$GWAS_path" ]]; then
    echo "Missing trait or GWAS input at row $GWAS_id of $GWAS_list" >&2
    exit 1
fi

qtl_name="${tissue}_${celltype}_${ancestry}"
besd_prefix="${qtl_name}_sceQTL"
eQTL_path="/path/to/${tissue}/sc-eQTL/results/besd/$besd_prefix"
bed_dir="/path/to/${tissue}/sc-eQTL/celltype_matrix/bed"
probe_path="/path/to/${tissue}/sc-eQTL/results/sig"
bed_name="${celltype}_${ancestry}_expr.bed"
probe_name="${qtl_name}_snp_probe_fdr05.list"
meta=/path/to/crossTissue_analysis/tissue_celltype_ancestry_ss_new.txt
QTLN=$(awk -v t="$tissue" -v a="$ancestry" -v ct="$celltype" \
    '($1 == t || $2 == t) && ($1 == a || $2 == a) && $3 == ct {print $4}' "$meta")
if ! [[ "$QTLN" =~ ^[1-9][0-9]*$ ]]; then
    echo "Expected one sample size for $qtl_name in $meta" >&2
    exit 1
fi

out_dir="/path/to/${tissue}/sceqtl_gwas_integration/COLOC/$result_subdir"
mkdir -p "$out_dir"
dir_local=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/data}}/coloc.XXXXXXXX")
publish_dir=""
cleanup() {
    rm -rf -- "$dir_local"
    if [[ -n "$publish_dir" ]]; then rm -rf -- "$publish_dir"; fi
}
trap cleanup EXIT
publish_dir=$(mktemp -d "$out_dir/.coloc.XXXXXXXX")

cp -- "$GWAS_path" "$bed_dir/$bed_name" "$dir_local/"
for ext in besd esi epi; do
    cp -- "$eQTL_path.$ext" "$dir_local/"
done
cd -- "$dir_local"

args=(--SMR "$SMR" --sumstats "$(basename -- "$GWAS_path")"
      --besd "$besd_prefix" --QTLName "$qtl_name" --QTLN "$QTLN"
      --expr_bed "$bed_name" --R_functions "$r_dir/CommonFunc.r" --tmp_dirt "$dir_local")
if [[ -n "$GWAS_N" && "$GWAS_N" != "NA" ]]; then
    args+=(--GWASN "$GWAS_N")
fi
if [[ "$mode" == fdr ]]; then
    cp -- "$probe_path/$probe_name" .
    args+=(--perm TRUE --probe_path "$probe_name")
else
    args+=(--perm FALSE)
fi

has_coloc_header() {
    [[ -s "$1" ]] && awk -F '\t' 'NR == 1 {
        exit !($1 == "probe_name" && $2 == "gene_name" && $3 == "pp0" &&
               $4 == "pp1" && $5 == "pp2" && $6 == "pp3" && $7 == "pp4")
    }' "$1"
}

for chr in {1..22}; do
    name="COLOC_${GWAS_trait}_${qtl_name}_${GWAS_ancestry}_chr${chr}"
    if has_coloc_header "$out_dir/$name.coloc"; then continue; fi
    Rscript "$r_dir/COLOC_perm.R" "${args[@]}" --chr "$chr" --out "$publish_dir/$name"
    if ! has_coloc_header "$publish_dir/$name.coloc"; then
        echo "Missing or invalid COLOC output for chromosome $chr" >&2
        exit 1
    fi
    mv -f -- "$publish_dir/$name.coloc" "$out_dir/$name.coloc"
done

echo "Completed COLOC for GWAS ${GWAS_trait}_${GWAS_ancestry}, eQTL $qtl_name"
