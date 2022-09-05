#!/bin/bash


echo "Enter the path of the folder where you will download the data:"
read n

echo "Enter the project name:"

read m

Rscript --vanilla retrieveData_bash.R $n $m
python tutorial_bash.py $n $m

