#!/bin/bash

BETAWARE_ROOT=$PWD
export BETAWARE_ROOT
PATH=$PATH:$BETAWARE_ROOT/bin

if [ $# -lt 2 ]
then
   echo $0 fastafile profilefile
   exit
fi

fasta=$1
profile=$2

echo "======= running  ==============================="
echo "betaware.py -t -f" $fasta "-p " $profile
betaware.py -t -f $fasta -p $profile
