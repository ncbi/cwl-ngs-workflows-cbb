#!/usr/bin/env bash

usage(){
cmdName=`basename $0`
cat << EOF
usage: ${cmdName} options

This script send jobs to the SGE

OPTIONS:
   -p      Project base directory
   -e      Conda environment to activate
EOF
}

PROJECT=""
CONDAENV=""
while getopts 'p:e:' OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        p)
            PROJECT="$OPTARG"
            ;;
        e)
            CONDAENV="$OPTARG"
            ;;
    esac
done

if [[ -z ${PROJECT} ]]
then
   echo "Project base directory is required"
   exit -1
fi

if [[ -z ${CONDAENV} ]]
then
   echo "Conda environment name is requiered"
   exit -1
fi

source $PROJECT/bin/conda/etc/profile.d/conda.sh
conda activate $CONDAENV

TMPFOLDER=`hexdump -n 16 -v -e '/1 "%02X"' /dev/urandom`

# For slurm temp folder
if [[ ! -z $SLURM_JOB_ID ]] & [[ -e "/lscratch/$SLURM_JOB_ID" ]]
then
    export TMPDIR=/lscratch/$SLURM_JOB_ID/$TMPFOLDER
else
    export TMPDIR=/tmp/$TMPFOLDER
fi
mkdir ${TMPDIR}
echo "Using TPM=${TMPDIR}"
echo "Setting env variables"
source $PROJECT/config/init.sh

echo "Executing command: ${@:5}"

CMD="$PROJECT/bin/conda/envs/$CONDAENV/bin/cwl-runner --on-error continue --no-container  --rm-tmpdir --tmp-outdir-prefix=$TMPDIR/ --tmpdir-prefix=$TMPDIR/ ${@:5}"
eval $CMD
rm -rf ${TMPDIR}

exit 0
