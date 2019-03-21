#!/bin/bash

##############################################################
##	PVAmplicon_install.sh									##
##	Creation 21/03/2019										##
##	Alexis Robitaille : robitaillea@students.iarc.fr		##
##	IARC, LYON												##
##	Last modification = 21/03/2019							##
##	Version 1.0												##
##############################################################
##	This script is used to install all the requirement to run PVAmpliconFinder tool on a Linux machine
##	For more information, please visit : https://github.com/SixEl27/PVAmpliconFinder

NC='\033[0m'
BLUE='\033[0;34m'
RED='\033[0;31m'
systime=`date`
echo -e "${RED}[$systime]${NC}"
echo -e "This script will install all the requirement to run PVAmpliconFinder tool on a Linux machine."
echo -e "It may take some time, please be patient..."


##	Get the machine type
MACHINE_TYPE=`uname -m`

##	Check if Python is installed and get the Python version
PYTHON_VERSION=$(python -V 2>&1 | grep -Po '(?<=Python )(.+)')
if [[ -z "$PYTHON_VERSION" ]]
then
    echo -e "Python is not yet installed. Please install at least Python version 2.7.0 on your machine.\n"
    echo -e "For more information, please visit : https://github.com/SixEl27/PVAmpliconFinder\n"
    exit
fi

##	Parse Python version
echo -e "Proceed installation for ${MACHINE_TYPE} machine type, using Python version ${PYTHON_VERSION}"

##########################
##	Conda installation	##
##########################
##	To get latest version of Miniconda, please visit https://conda.io/en/latest/miniconda.html

cd $PWD/program

##	Install the Miciconda version corresponding to the Phyton version installed on the system
if [[ "$parsedVersion" -lt "270" ]]										##	Python version < 2.7.0
then 
    echo -e "Python is not yet installed, or the version is too old. Please install or update your current python version."
    echo -e "For more information, please visit : https://github.com/SixEl27/PVAmpliconFinder"
    exit
elif [[ "$parsedVersion" -lt "300" && "$parsedVersion" -ge "270" ]]		##	Python version => 2.7.0 && Python version < 3.0.0
then
	echo -e "Downloading and installation of CONDA."
	if [ ${MACHINE_TYPE} == 'x86_64' ]
	then
		wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
		chmod +x $PWD/program/Miniconda2-latest-Linux-x86_64.sh
		bash $PWD/program/Miniconda2-latest-Linux-x86_64.sh -b -p $PWD/program/miniconda
	else
		wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86.sh
		chmod +x $PWD/program/Miniconda2-latest-Linux-x86.sh
		bash $PWD/program/Miniconda2-latest-Linux-x86.sh -b -p $PWD/program/miniconda
		echo "${RED}Please considere that under 32bits version, you need to manually install PaPaRa"
		echo "For more information, please visit https://cme.h-its.org/exelixis/web/software/papara/index.html${NC}"
	fi
else																	##	Python version => 3.0.0
then
	if [ ${MACHINE_TYPE} == 'x86_64' ]
	then
		wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
		chmod +x $PWD/program/Miniconda3-latest-Linux-x86_64.sh
		bash $PWD/program/Miniconda3-latest-Linux-x86_64.sh -b -p $PWD/program/miniconda
	else
		wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh
		chmod +x $PWD/program/Miniconda3-latest-Linux-x86.sh
		bash $PWD/program/Miniconda3-latest-Linux-x86.sh -b -p $PWD/program/miniconda
		echo "${RED}Please considere that under 32bits version, you need to manually install PaPaRa"
		echo "For more information, please visit https://cme.h-its.org/exelixis/web/software/papara/index.html${NC}"
	fi
fi

##	Add miniconda into $PATH environment variable
export PATH="$PWD/miniconda/bin:$PATH"

##	Update conda
conda update -y -n base -c defaults conda

##	Add bioconda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

##########################
##	Install the tools	##
##########################
##	All is installed with conda
conda install -y fastqc multiqc trim-galore vsearch blast raxml cap3 krona

##	PaPaRa is not available on Conda
wget https://cme.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz
tar -xvzf papara_nt-2.5-static_x86_64.tar.gz
mv $PWD/program/papara_nt-2.5-static_x86_64 $PWD/program/papara
chmod +x $PWD/program/papara
##	If any issue with PaPaRa, need to be build manually : https://cme.h-its.org/exelixis/web/software/papara/index.html

cd ..

systime=`date`
echo -e "${RED}[$systime]${NC}"
