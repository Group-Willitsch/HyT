#!/bin/bash

ls ./*/output/out*.dat

for filename in ./*/output/out*.dat; do
	echo "$filename"			# bare filename
#	echo "${filename%.*}.zip"		# filename with zip extension	
#	zip "${filename%.*}.zip" "$filename" 	# zip all
done

for filename in ./*/output/out*.zip; do
	unzip "${filename%.*}.zip" 		# unzip all
done

for filename in ./Stark_deceleration_with_MATLAB_parfor/sequences/*/*.zip; do
#	echo "$filename"
	unzip "${filename}"
done
# find . -name 'out*.dat' -delete 		# find, then delete all the appropriate .dat files

