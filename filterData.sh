#!/bin/bash

# Go to data folder 
cd Data 

# Use awk to filter and write down 
var=$@
echo "$var"
awk -v FS=';' '$19 ~ "HA-" {print $19 "," $3 "," $12 }' AFLSA208A_PAN_TRYPSIN.csv > filt_AFLSA.csv