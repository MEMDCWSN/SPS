#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

$(dirname $0)/run_csps.sh $1 csps.params 0
