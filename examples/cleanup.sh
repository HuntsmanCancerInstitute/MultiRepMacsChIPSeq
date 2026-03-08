#!/bin/bash

echo
echo "########## Cleaning up example run directories ###########"

for DIR in se pe tup1_indiv pe_indep pe_nodup pe_gapped pe_nobam pe_nocon pe_nogenome pe_normsize cutsite rpd3 tup1 Rpd3_Tup1_merge
do
	if [[ -e $DIR ]]
	then
		echo "Removing $DIR"
		rm -rf $DIR
	fi
done

echo "Complete"
