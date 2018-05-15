#!/bin/bash

file_names=(plaq_occ topo_chrg avg_plaq_occ)

for item in ${file_names[*]}
do
	echo $item
	g++-7 stats.cpp -o stats -std=c++11
	./stats $item
done

echo All done
