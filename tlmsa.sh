#!/bin/bash


echo "Enter the path for data:"
read n

echo "Enter project name:"

read m

Rscript --vanilla retrieveData_bash.R $n $m
python tutorial_bash.py $n $m

