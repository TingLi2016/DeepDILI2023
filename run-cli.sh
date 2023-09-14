#!/bin/sh
set -e

die() { echo "$*" 1>&2 ; exit 1; }

[ $# -eq 1 ] || die "Expected arguments: image"

IMAGE=$1

docker run -i --rm "$IMAGE" /deep-dili/main_cli.py /data/training/online-objs.pkl /dev/stdin
