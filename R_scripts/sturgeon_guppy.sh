#! /bin/bash

echo "sturgeon inputtobed --margin 50 -i $1 -o $1 -s guppy --probes-file $2"
sturgeon inputtobed --margin 50 -i $1 -o $1 -s guppy --probes-file $2 ;
echo "sturgeon predict -p --i $1 -o $1 -m $3"
sturgeon predict -p --i $1 -o $1 -m $3 ;
