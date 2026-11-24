#!/bin/bash
# Run Cell Ranger for the donor on line i of sample.list.
set -euo pipefail

module load cellranger

ref="/path/to/genome_reference/GRCh38_10X/refdata-gex-GRCh38-2020-A"
tissue=$1
project=$2
i=$3
CR_dir="/path/to/scRNA/${tissue}/raw_data/${project}"

if ! [[ "$i" =~ ^[1-9][0-9]*$ ]]; then
    echo "Expected a positive sample.list line number, got: $i" >&2
    exit 1
fi
sample=$(sed -n "${i}p" "${CR_dir}/sample.list")
sample=${sample%$'\r'}
if [[ -z "$sample" || "$sample" == */* || "$sample" == "." || "$sample" == ".." ]]; then
    echo "Missing or invalid donor on line $i of ${CR_dir}/sample.list" >&2
    exit 1
fi

# By default the FASTQ library prefix is the donor ID. For multiple libraries
# or different library IDs, provide sample_fastq_prefixes.tsv in the project
# directory: donor ID<TAB>prefix1,prefix2. List only donors needing overrides.
prefix_csv="$sample"
prefix_map="${CR_dir}/sample_fastq_prefixes.tsv"
if [[ -f "$prefix_map" ]]; then
    if ! prefix_csv=$(awk -F '\t' -v donor="$sample" '
        ("x" $1) == ("x" donor) { print $2; n++ }
        END { if (n > 1) exit 1; if (n == 0) print donor }
    ' "$prefix_map"); then
        echo "Could not read a unique library mapping for donor $sample from $prefix_map" >&2
        exit 1
    fi
    prefix_csv=${prefix_csv%$'\r'}
fi
if [[ -z "$prefix_csv" || "$prefix_csv" == ,* || "$prefix_csv" == *, || "$prefix_csv" == *,,* || "$prefix_csv" == */* ]]; then
    echo "Invalid FASTQ library prefix list for donor $sample: $prefix_csv" >&2
    exit 1
fi
IFS=, read -r -a prefixes <<< "$prefix_csv"

dir_local=$(mktemp -d "/data/cellranger.XXXXXXXX")
trap 'rm -rf -- "$dir_local"' EXIT
shopt -s nullglob
fastq_pattern='^(.+)_S[0-9]+(_L[0-9]{3})?_[RI][0-9]+_[0-9]{3}\.fastq\.gz$'
for prefix in "${prefixes[@]}"; do
    selected=()
    for fastq in "$CR_dir/${prefix}"_S*.fastq.gz; do
        name=${fastq##*/}
        if [[ "$name" =~ $fastq_pattern && "${BASH_REMATCH[1]}" == "$prefix" ]]; then
            selected+=("$fastq")
        fi
    done
    if (( ${#selected[@]} == 0 )); then
        echo "No matching FASTQs for donor $sample, library $prefix in $CR_dir" >&2
        exit 1
    fi
    cp -- "${selected[@]}" "$dir_local/"
done

cd "$dir_local"
cellranger count --id="$sample" \
    --create-bam=true \
    --transcriptome="$ref" \
    --fastqs=./ \
    --sample="$prefix_csv" \
    --localcores=12 \
    --localmem=56

cp -r -- "$sample" "$CR_dir/"
