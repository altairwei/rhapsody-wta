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

GTF_FILES="$1"
shift
FASTA_FILES="$@"

filename=$(basename -- "$GTF_FILES")
extension="${GTF_FILES##*.}"
filename="${filename%.*}"

mkdir -p data/reference_indexes

# Set `--genomeSAsparseD 2` to reduce RAM requirement when mapping.
STAR \
    --runThreadN 16 \
    --limitGenomeGenerateRAM 103079215104 \
    --genomeSAsparseD 2 \
    --runMode genomeGenerate \
    --genomeDir data/reference_indexes \
    --genomeFastaFiles "$FASTA_FILES" \
    --sjdbGTFfile "$GTF_FILES"

tar --transform "s#^\.#${filename}#" \
    -cvf - -C data/reference_indexes . | pigz -p 16 > data/reference_sequences/${filename}.tar.gz