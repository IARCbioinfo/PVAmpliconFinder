#!/bin/bash

##############################################################
##	PVAmplicon_install.sh									##
##	Creation 27/03/2019										##
##	Alexis Robitaille : robitaillea@students.iarc.fr		##
##	IARC, LYON												##
##	Last modification = 27/03/2019							##
##	Version 1.0												##
##############################################################
##	This script is used to install configure the lineage file and PaVE database used by PVAmpliconFinder
##	For more information, please visit : https://github.com/SixEl27/PVAmpliconFinder

##	Lineages file from https://github.com/zyxue/ncbitax2lin

#~ Check if already present

#~ Unzip the file (already present - given in the GitHub repository)
gunzip $PWD/databases/lineages/lineages*.csv.gz
#~ Get only viruses related lineages
grep -i "virus"  $PWD/databases/lineages/lineages*.csv >  $PWD/databases/lineages/lineagesVirus.csv
#~ Remove uncompress initial file as it takes a lot of disk space
rm  $PWD/databases/lineages/lineages-*.csv
