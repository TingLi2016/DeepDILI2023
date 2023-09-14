#!/bin/sh
set -e

die() { echo "$*" 1>&2 ; exit 1; }

[ $# -eq 1 ] || die "Expected arguments: image"

IMAGE=$1

docker run --rm -v $(pwd)/data:/data "$IMAGE" \
  /deep-dili/main_training.py \
  /data \
  /data/training \
  /data/training/online-objs.pkl
