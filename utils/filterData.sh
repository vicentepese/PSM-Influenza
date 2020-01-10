#!/bin/bash

# Go to data folder 
cd ../DataMASS 

# Use awk to filter and write down each file
var="$(ls)"

for i in $var ; do
    echo "$i"
    echo "Filtering $i "
    if [[ $i == *"ARP"* ]]; then 
        awk -v FS=';' '$19 ~ "HA-" {print $19 "," $3 "," $12 "," "ARP"}' $i > ../FiltDataMASS/$i
    elif [[ $i == *"PAN"* ]]; then 
        awk -v FS=';' '$19 ~ "HA-" {print $19 "," $3 "," $12 "," "PAN"}' $i > ../FiltDataMASS/$i
    else
         awk -v FS=';' '$19 ~ "hemagglutinin" {print "HA" "," $3 "," $12 "," "FOC"}' $i > ../FiltDataMASS/$i
    fi
done
