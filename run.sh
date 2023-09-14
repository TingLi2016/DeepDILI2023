#!/bin/sh
set -e

# docker pull scharris/deep-dili

IMAGE=${1:-"scharris/deep-dili:latest"}

(docker rm -vf deep-dili 2> /dev/null) || echo "no running container found"
docker run --name deep-dili -d -v $(pwd)/data:/data --network=host "$IMAGE"
