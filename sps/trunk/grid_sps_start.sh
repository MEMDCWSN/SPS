#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

$(dirname $0)/run_specnets.sh $1 sps.params 0 sps
