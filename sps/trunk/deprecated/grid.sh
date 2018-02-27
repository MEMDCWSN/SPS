#!/bin/bash
#
# author: Phil Bouchard <pbouchard@cs.ucsd.edu>
#         May 2010

# arguments
while getopts 'h' OPTION
do
	case $OPTION in
	h)	echo "Usage: $(basename $0): [PROJECT DIR] [MATLAB DIR] [MAX NO PROCESSES]" >&2
		exit 2
		;;
	esac
done

# default values
MAX=${3:-128}
DIR=${2:-~/matlab}

EXE=$(dirname $0)
ROOT=${1:-$(pwd)}

# definitions
function halt_over
{
	while (( $(qstat -u $USER | wc -l) >= $1 + 2 )); do sleep 1; done
}

function launch_qsub
{
	for i in $SCRIPT; do
		halt_over $MAX
		cd $(dirname $i)
		qsub -v LD_LIBRARY_PATH $EXE/grid_$1_$2.sh $DIR
	done

	wait
}

function launch_bash
{
	for i in $SCRIPT; do
		halt_over $MAX
		cd $(dirname $i)
		bash $i
	done

	wait
}

# execution
echo
echo "* 1 of 4:"
SCRIPT=$(find $ROOT -name sps.params)
launch_qsub sps start
halt_over 0

echo
echo "* 2 of 4:"
SCRIPT=$(find $ROOT -name run_jobs.sh)
launch_bash
halt_over 0

echo
echo "* 3 of 4:"
SCRIPT=$(find $ROOT -name sps.params)
launch_qsub sps resume
halt_over 0

echo
echo "* 4 of 4:"
SCRIPT=$(find $ROOT -name csps.params)
launch_qsub csps start
halt_over 0
