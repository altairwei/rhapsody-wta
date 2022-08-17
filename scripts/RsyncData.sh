#!/usr/bin/env bash

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

HOST=$1
DEST=$2
SOURCES=$(echo -n 0DPI-MOCK-{1,2} {1,2,3}DPI-{MOCK,PNR2,TR4}-{1,2} 3DPI-PNR2-3 | xargs -d' ' -I{} echo -n "$HOST{} ")

rsync -avr --exclude "*.BAM" --exclude "*/SeuratAnalysis/" --exclude "*/QualityControl/" ${SOURCES} ${DEST}
