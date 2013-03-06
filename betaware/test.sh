#!/bin/bash

BETAWARE_ROOT=$PWD
export BETAWARE_ROOT
PATH=$PATH:$BETAWARE_ROOT/bin

echo "======= running  ==============================="
echo "betaware.py -f example/1qj8.fasta -p example/1qj8.prof"
betaware.py -f example/1qj8.fasta -p example/1qj8.prof

echo "======= running  ==============================="
echo "betaware.py -t -f example/12ca.fasta -p example/12ca.prof"
betaware.py -t -f example/12ca.fasta -p example/12ca.prof

echo "======= running  ==============================="
echo "betaware.py -o example/1qj8.out -f example/1qj8.fasta -p example/1qj8.prof"
betaware.py -o example/1qj8.out -f example/1qj8.fasta -p example/1qj8.prof

echo "======= running  ==============================="
echo "betaware.py -a ACDEFGHIKLMNPQRSTVWY -f example/1qj8.fasta -p example/1qj8ALPH.prof"
betaware.py -a ACDEFGHIKLMNPQRSTVWY -f example/1qj8.fasta -p example/1qj8ALPH.prof

echo "======= running  ==============================="
echo "betaware.py -s 0.1 -f example/1qj8.fasta -p example/1qj8.prof"
betaware.py -s 0.1 -f example/1qj8.fasta -p example/1qj8.prof

echo "======= running  ==============================="
echo "betaware.py -s 0.95 -f example/12ca.fasta -p example/12ca.prof"
betaware.py -s 0.95 -f example/12ca.fasta -p example/12ca.prof