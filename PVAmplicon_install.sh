#!/bin/bash

##############################################################
##	PVAmplicon_install.sh									##
##	Creation 21/03/2019										##
##	Alexis Robitaille : robitaillea@students.iarc.fr		##
##	IARC, LYON												##
##	Last modification = 26/03/2019							##
##	Version 1.0												##
##############################################################
##	This script is used to install all the requirement to run PVAmpliconFinder tool on a Linux machine
##	For more information, please visit : https://github.com/SixEl27/PVAmpliconFinder

usage="$(basename "$0") [-h]
This script is used to install all the requirement to run PVAmpliconFinder tool on a Linux machine.
It does not take any argument. An internet connexion is required.
Please note that tool installation is managed by conda, that will be installed in your system.
All executable file will be accessible in your PATH.
If you wish to manually install the tools, please visit https://github.com/SixEl27/PVAmpliconFinder to get the complete required tool list."

while getopts ':h' option; do
  case "$option" in
    h) echo -e "$usage"
       exit
       ;;
   esac
done
shift $((OPTIND - 1))
       
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
parsedVersion=$(echo "${PYTHON_VERSION//./}" | grep -P '^\d+' -o)


##########################
##	Conda installation	##
##########################
##	To get latest version of Miniconda, please visit https://conda.io/en/latest/miniconda.html

cd $PWD/program

##	Install the Miciconda version corresponding to the Phyton version installed on the system
if [[ "$parsedVersion" -lt "2700" ]]										##	Python version < 2.7.0
then 
    echo -e "Python is not yet installed, or the version is too old. Please install or update your current python version."
    echo -e "For more information, please visit : https://github.com/SixEl27/PVAmpliconFinder"
    exit
elif [[ "$parsedVersion" -lt "3000" && "$parsedVersion" -ge "2700" ]]		##	Python version => 2.7.0 && Python version < 3.0.0
then
	echo -e "Downloading and installation of CONDA."
	if [ ${MACHINE_TYPE} == 'x86_64' ]
	then
		wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
		chmod +x $PWD/Miniconda2-latest-Linux-x86_64.sh
		bash $PWD/Miniconda2-latest-Linux-x86_64.sh -b -p $PWD/miniconda
	else
		wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86.sh
		chmod +x $PWD/Miniconda2-latest-Linux-x86.sh
		bash $PWD/Miniconda2-latest-Linux-x86.sh -b -p $PWD/miniconda
		echo "${RED}Please considere that under 32bits version, you need to manually install PaPaRa"
		echo "For more information, please visit https://cme.h-its.org/exelixis/web/software/papara/index.html${NC}"
	fi
else																	##	Python version => 3.0.0
	if [ ${MACHINE_TYPE} == 'x86_64' ]
	then
		wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
		chmod +x $PWD/Miniconda3-latest-Linux-x86_64.sh
		bash $PWD/Miniconda3-latest-Linux-x86_64.sh -b -p $PWD/miniconda
	else
		wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh
		chmod +x $PWD/Miniconda3-latest-Linux-x86.sh
		bash $PWD/Miniconda3-latest-Linux-x86.sh -b -p $PWD/miniconda
		echo "${RED}Please considere that under 32bits version, you need to manually install PaPaRa"
		echo "For more information, please visit https://cme.h-its.org/exelixis/web/software/papara/index.html${NC}"
	fi
fi

##	Add miniconda into $PATH environment variable
export PATH="$PWD/miniconda/bin:$PATH"

##	Update conda
#~ conda update -y -c defaults conda

##	Add bioconda channels
conda config --add channels defaults/label/cf201901
conda config --add channels bioconda/label/cf201901
conda config --add channels conda-forge/label/cf201901

##########################
##	Install the tools	##
##########################
##	All is installed with conda 

##	Tools and their dependencies
#~ conda install -y fastqc multiqc trim-galore vsearch blast raxml cap3 krona libxml2 gcc_linux-64 gxx_linux-64 gfortran_linux-64 
##	Perl packages and their dependencies
#~ conda install -y perl-padwalker perl-xml-libxml perl-libxml-perl perl-bioperl  perl-getopt-long perl-math-round perl-statistics-basic perl-list-moreutils perl-module-build perl-bioperl-run perl-text-csv

##	Evrything in once
conda install -y fastqc multiqc trim-galore vsearch blast raxml cap3 krona libxml2 gcc_linux-64 gxx_linux-64 gfortran_linux-64 perl-padwalker perl-xml-libxml perl-libxml-perl perl-bioperl perl-getopt-long perl-math-round perl-statistics-basic perl-list-moreutils perl-module-build perl-bioperl-run perl-text-csv

##	PaPaRa is not available on Conda
#~ wget https://cme.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz
#~ tar -xvzf papara_nt-2.5-static_x86_64.tar.gz
#~ mv $PWD/papara_static_x86_64 $PWD/papara
#~ chmod +x $PWD/papara

##	Go back to root of PVAmpliconFinder
cd $PWD/..

##	If any issue with PaPaRa, need to be build manually : https://cme.h-its.org/exelixis/web/software/papara/index.html
export PATH="$PWD/program:$PATH"

systime=`date`
echo -e "${BLUE}Done!${NC}"
echo -e "${RED}[$systime]${NC}"

##	How to remove everything correctly
##	Everything installed with conda
#~ conda uninstall -y fastqc multiqc trim-galore vsearch blast raxml cap3 krona libxml2 gcc_linux-64 gxx_linux-64 gfortran_linux-64 perl-padwalker perl-xml-libxml perl-libxml-perl perl-bioperl  perl-getopt-long perl-math-round perl-statistics-basic perl-list-moreutils perl-module-build perl-bioperl-run perl-text-csv
##	Conda itself
#~ rm -rf $PWD/program/*

#~ Also remove the hidden .condarc file and .conda and .continuum directories which may have been created in the home directory with
