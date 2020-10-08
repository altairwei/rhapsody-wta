#!/usr/bin/env bash

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

IMAGE=$1
LAYER_HASH=$2
TARGET=$3

tar --to-command="tar -xf - $TARGET" -xf $IMAGE $LAYER_HASH/layer.tar