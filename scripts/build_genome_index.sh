#!/usr/bin/env bash
#PBS -l mem=96gb,nodes=1:ppn=16,walltime=48:00:00
#PBS -q fat
#PBS -N genome-index
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR

source activate wta

GTF_FILE="$1"
shift
FASTA_FILES="$@"

filename=$(basename -- "$GTF_FILE")
extension="${GTF_FILE##*.}"
filename="${filename%.*}"

mkdir -p data/reference_indexes

# Set `--genomeSAsparseD 2` to reduce RAM requirement when mapping.
# You can't add double quotes to the value of `--genomeFastaFiles`,
#   otherwise the parameter won't be parsed properly.
STAR \
    --runThreadN 16 \
    --limitGenomeGenerateRAM 103079215104 \
    --genomeSAsparseD 2 \
    --runMode genomeGenerate \
    --genomeDir data/reference_indexes \
    --genomeFastaFiles $FASTA_FILES \
    --sjdbGTFfile $GTF_FILE

tar --transform "s#^\.#${filename}#" \
    -cvf - -C data/reference_indexes . | pigz -p 16 > data/${filename}.tar.gz