# PVAmpliconFinder

**PVampliconFinder** is a data analysis workflow designed to rapidly identify and classify known and potentially new *papilliomaviridae* sequences from amplicon deep-sequencing with degenerated papillomavirus (PV) primers.

## Description

PVampliconFinder is based on alignment similarity metrics, but also consider molecular evolution time for an improved identification and taxonomic classification of novel PVs. The final output of the tool includes a list of fully characterized putatively new papillomaviriade sequences, as well as graphical representations of relative abundance of the virome sequence diversity in the tested samples. 

## Prerequisites

The PVampliconFinder workflow is designed for the analysis of sequencing reads generated from **paired-end sequencing** of DNA amplified using degenerated primers targeting specifically the L1 sequence of papillomaviruses ([Chouhy *et al.*, 2010](https://www.ncbi.nlm.nih.gov/pubmed/19948351),[Forslund *et al.*, 1999](https://www.ncbi.nlm.nih.gov/pubmed/10501499),[Forslund *et al.*, 2003](https://www.ncbi.nlm.nih.gov/pubmed/12798239)).

### Dependencies

#### Programming Language

- Bash/Shell
- Perl
- Python

#### External Software

> For now the following tools need to be manually downloaded and installed, and corresponding executable must be present in the [PATH environment variable](http://www.linfo.org/path_env_var.html). Dependencies are expected to be soon available in this GitHub repository.

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [VSEARCH](https://github.com/torognes/vsearch)
- [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [RaxML-EPA](https://cme.h-its.org/exelixis/web/software/epa/index.html)
- [PaPaRa](https://cme.h-its.org/exelixis/web/software/papara/index.html)
- [CAP3](http://seq.cs.iastate.edu/cap3.html)
- [KRONA](https://github.com/marbl/Krona/wiki)

### Databases

> For now the following database must be manually downloaded. Please settle the ```~/.ncbirc file``` to specify the location of the nt database to the system.

- [blastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) : It is advised to use the NCBI script [update_blastdb.pl](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl). Required databases are :
  - nt
  - taxdb
- [ncbitax2lin](https://github.com/zyxue/ncbitax2lin) : A lineage file is available on the GitHub webpage, but a more recent one can be created manually following the instruction
- [PaVE](https://pave.niaid.nih.gov/#search/search_database/kw?dbNamespace=Genomes&includeNR=true&refCloneOnly=false&sort=Locus_ID&sortType=true&page=600&start=1&showTable=1&) : Select All > Download Fasta

Once the lineage file downloaded, please launch the following command :
```
grep -i "virus" lineages* > lineages-virus.csv
```
Once the PaVE database download, please launch the following command :
```
gettaxidByAcc.sh
```

## Input
  | Type      | Description     |
  |-----------|---------------|
  | -d        | PATH to input fastq directory|

> tests files can be found .....

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
| -p   | NA | Adapter seuqence of Read1 (in case the adapter have been sequenced) |
| -q    | NA | Adapter seuqence of Read2 (in case the adapter have been sequenced)|
| -t    | 2 | Number of threads |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| -h   | Display help |

## Usage
  ```
sh amplicon_processing_HPV_Vlast.sh [-h] [-t threads] [-b "nt" database] [-f info_file] [-p adapter R1] [-q adapter R2] -s fastq_files_suffix -d input_dir -o output_dir
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
  | Alexis Robitaille    | robitaillea@students.iarc.fr | Developer to contact for support (link to specific gitter chatroom) |
  | Magali Olivier    | olivierm@iarc.fr |  |
  | Massimo Tommasino    | tommasinom@iarc.fr |  |

## Versioning

Version 1.0

## Authors

* **Alexis Robitaille** - *Initial work* - [IARC bioinformatic platform](https://github.com/IARCbioinfo)

## License

This project is licensed under ...

## Acknowledgments

## References

[References](docs/references.md)

## FAQ




