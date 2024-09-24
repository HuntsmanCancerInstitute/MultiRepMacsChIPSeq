#!/bin/bash

echo
echo "########## Cleaning up example run directories ###########"

if [[ -e se ]]
then
	echo "Removing se"
	rm -rf se
fi

if [[ -e pe ]]
then
	echo "Removing pe"
	rm -rf pe
fi

if [[ -e tup1_indiv ]]
then
	echo "Removing tup1_indiv"
	rm -rf tup1_indiv
fi

if [[ -e pe_indep ]]
then
	echo "Removing pe_indep"
	rm -rf pe_indep
fi

if [[ -e pe_nodup ]]
then
	echo "Removing pe_nodup"
	rm -rf pe_nodup
fi

if [[ -e pe_gapped ]]
then
	echo "Removing pe_gapped"
	rm -rf pe_gapped
fi

if [[ -e pe_nocon ]]
then
	echo "Removing pe_nocon"
	rm -rf pe_nocon
fi

if [[ -e pe_nogenome ]]
then
	echo "Removing pe_nogenome"
	rm -rf pe_nogenome
fi

if [[ -e cutsite ]]
then
	echo "Removing cutsite"
	rm -rf cutsite
fi

if [[ -e rpd3 ]]
then
	echo "Removing rpd3"
	rm -rf rpd3
fi

if [[ -e tup1 ]]
then
	echo "Removing tup1"
	rm -rf tup1
fi

if [[ -e Rpd3_Tup1_merge ]]
then
	echo "Removing Rpd3_Tup1_merge"
	rm -rf Rpd3_Tup1_merge
fi

echo "Complete"
