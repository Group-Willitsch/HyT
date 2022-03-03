#!/bin/bash

ls ../*/*/T2jump_test.out
for filename in ../*/*/*/*/T2jump.out; do
	echo "${filename%.*}"
	echo "${filename%.*}_c.out"  # new filename is same with _c added
#	grep -v '^#' $filename | grep  -v '^\[' 	| sed -e 's/^[ \t]*//' > "${filename%.*}_c.out" # remove comments and first screwed line
	sed -e 's/^[ \t]*//' $filename | tr -s ' ' | sed "s/^\[1\]/#\[1\]/g" | sed "s/^#.*/#/g" > "${filename%.*}_c.out" # remove comments and first screwed line
	
done


