#!/usr/bin/env bash

usage(){
cmdName=`basename $0`
cat << EOF
usage: $cmdName options

This script send jobs to the Torque queue system

OPTIONS:
   -w      Working dir
   -n      Number of processors
   -c      Number of Threads
   -p      Partition of queue name
   -t      Wall time (format: DD-HH-MM-SS)
   -r      Command to run
   -d      Dependency job ID
   -J      Job name
EOF
}


# HEAD INSTANCE

w=`pwd`
n=1
c=1
p=""
t="04:00:00"
r=""
d=""
J=""

while getopts 'hn:c:p:t:r:w:d:' OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        w)
            w="$OPTARG"
            ;;
        d)
            d="$OPTARG"
            ;;
        p)
            p="$OPTARG"
            ;;
        n)
            n="$OPTARG"
            ;;
        c)
            c="$OPTARG"
            ;;
        t)
            t="$OPTARG"
            ;;
        r)
            r="$OPTARG"
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

cd $w || { echo "ERROR: Could not change to working dir" ; exit 1; }

sbatch -c $c
