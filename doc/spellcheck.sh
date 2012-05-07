#!/bin/bash
for f in $( ls *.tex)
do
	echo $f;
	aspell -l en_US -c $f
done;
	
