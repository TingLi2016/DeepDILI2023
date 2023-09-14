#!/bin/bash
SCRIPTDIR="$(cd "$(dirname "$0")"; pwd)";

rm "$SCRIPTDIR"/*/*/*.csv
rm "$SCRIPTDIR"/probabilities_output/*.csv
