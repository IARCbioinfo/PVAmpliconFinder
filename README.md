# PVAmpliconFinder

**PVampliconFinder** is a data analysis workflow designed to rapidly identify and classify known and potentially new *papilliomaviridae* sequences from amplicon deep-sequencing with degenerated papillomavirus (PV) primers.

Robitaille, A., Brancaccio, R.N., Dutta, S. et al. PVAmpliconFinder: a workflow for the identification of human papillomaviruses from high-throughput amplicon sequencing. BMC Bioinformatics 21, 233 (2020). https://doi.org/10.1186/s12859-020-03573-8

## Description

PVampliconFinder is based on alignment similarity metrics, but also consider molecular evolution time for an improved identification and taxonomic classification of novel PVs. The final output of the tool includes a list of fully characterized putatively new papillomaviriade sequences, as well as graphical representations of relative abundance of the virome sequence diversity in the tested samples. 

## Prerequisites

The PVampliconFinder workflow is designed for the analysis of sequencing reads generated from **paired-end sequencing** of DNA amplified using degenerated primers targeting specifically the L1 sequence of papillomaviruses ([Chouhy *et al.*, 2010](https://www.ncbi.nlm.nih.gov/pubmed/19948351),[Forslund *et al.*, 1999](https://www.ncbi.nlm.nih.gov/pubmed/10501499),[Forslund *et al.*, 2003](https://www.ncbi.nlm.nih.gov/pubmed/12798239)).

## Installation

[Python2.7](https://www.python.org/download/releases/2.7/) or higher and [Perl v5.22.1](https://www.perl.org/get.html) or higher are required.

The tool has been created under UNIX environment, but installing [clang_osx-64](https://anaconda.org/conda-forge/clang_osx-64), [clangxx_osx-64](https://anaconda.org/conda-forge/clangxx_osx-64) and [gfortran_osx-64](https://anaconda.org/conda-forge/gfortran_osx-64) with conda should provide a functional environment on Mac.

### Automatic installation

PVAmpliconFinder come with a SHELL script [PVAmplicon_install.sh](PVAmplicon_install.sh) that will proceed with the downloading and the installation of all the software required to run PVAmpliconFinder.

```
source PVAmplicon_install.sh [-h] [-p conda installation path]

Description :

    -h  display help
    -p	directory installation path (will be created if not already existing) - default value : PVAmpliconFinder/program

```

**PVAmpliconFinder rely on [Bioconda](https://bioconda.github.io/) to install the software and associated dependencies**

> For 32bits system, [PaPaRa](https://cme.h-its.org/exelixis/web/software/papara/index.html) available binary file is not functionnal, as specified on the webpage of the tool. You need to install manually PaPara following the instruction, and put the binary file in PVAmpliconFinder/program. Note that the binary file must be named "papara".

### Manual installation

The list of tools used by PVAmpliconFinder can be manually downloaded and installed, and corresponding **executable must be present in the [PATH environment variable](http://www.linfo.org/path_env_var.html)**.

> Please note that [PaPaRa](https://cme.h-its.org/exelixis/web/software/papara/index.html) binary file must be named "papara".

### List of software

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [VSEARCH](https://github.com/torognes/vsearch)
- [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [RaxML-EPA](https://cme.h-its.org/exelixis/web/software/epa/index.html)
- [PaPaRa](https://cme.h-its.org/exelixis/web/software/papara/index.html)
- [CAP3](http://seq.cs.iastate.edu/cap3.html)
- [KRONA](https://github.com/marbl/Krona/wiki)

## Databases

### NCBI databases

PVAmpliconFinder need the **nt** and **taxdb** NCBI databases to work properly. You can find thoses databases at the following ftp adress : ftp://ftp.ncbi.nlm.nih.gov/blast/db/. Note that the taxonomy file must be correctly located.

It is advised to use the NCBI script [update_blastdb.pl](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl) to facilitate the installation of the databases.

**Once downloaded and installed, please check that the ```~/.ncbirc file``` is present and point to the correct NCBI nt database location.**

### List of other databases

- [ncbitax2lin](https://github.com/zyxue/ncbitax2lin) 
- [PaVE](https://pave.niaid.nih.gov/)

## Input
  | Type      | Description     |
  |-----------|---------------|
  | -d        | PATH to input fastq directory|

> tests files can be found [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA555194)

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| -s    | pool | suffix of fastq filename |
| -o    | PV_Amplicon_output | PATH to output directory |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| -f   | NA | Tabular file containing information about the samples |
| -b    | nt | Name of the local "nt" blast database |
| -i	| 98	| Threshold of percentage of identity used for the de-novo centroid-based clustering	|
| -t    | 2 | Number of threads |


  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| -h   | Display help |

## Usage
```
sh PVAmpliconFinder.sh [-h] [-t threads] [-b "nt" database] [-f info_file] [-i identity thershold] -s fastq_files_suffix -d input_dir -o output_dir
```
## Output

  | Type      | Description     |
  |-----------|---------------|
  | QC report    | Report on FastQ file quality, before and after trimming |
  | Diversity by tissu    | Excel table of taxonomically classified PV species identified in the samples |
  | Table summary    | Excel table of reads metics |
  | Table putative Known viruses    | Excel table of putative known viruses identified in the samples |
  | Table putative New viruses    | Excel table of putative new viruses identified in the samples | 
  | Putative Known viruses    | Fasta files of putative known viruses ssequences identified in the samples | 
  | Putative New viruses    | Fasta files of putative new viruses ssequences identified in the samples | 
  | KRONA Megablast    | Directory of KRONA graphical representations of the unormalized abundance of viruses identified by Megablast in the samples |
  | KRONA BlastN    | Directory of KRONA graphical representations of the unormalized abundance of viruses identified by BlastN in the samples |
  | KRONA RaxML    | Directory of KRONA graphical representations of the unormalized abundance of viruses identified by RaxML-EPA in the samples |
  | Log file    | File of the logs | 
  
## Detailed description of the output

[Detailed description of the output](docs/output_description.md)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Alexis Robitaille    | alexis.robitaille@orange.fr | Developer to contact for support |
  | Magali Olivier    | olivierm@iarc.fr |  |
  | Massimo Tommasino    | tommasinom@iarc.fr |  |

## Versioning

Version 1.0

## Authors

* **Alexis Robitaille** - [IARC bioinformatic platform](https://github.com/IARCbioinfo)

## License

This project is licensed under GPL-3.0.

## Acknowledgments

## References

[References](docs/references.md)

## FAQ




