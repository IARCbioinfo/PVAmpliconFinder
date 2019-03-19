# PVAmpliconFinder

**PVampliconFinder** is a data analysis workflow designed to rapidly identify and classify known and potentially new *papilliomaviridae* sequences from amplicon deep-sequencing with degenerated papillomavirus (PV) primers.

## Description

PVampliconFinder is based on alignment similarity metrics, but also consider molecular evolution time for an improved identification and taxonomic classification of novel PVs. The final output of the tool includes a list of fully characterized putatively new papillomaviriade sequences, as well as graphical representations of relative abundance of the virome sequence diversity in the tested samples. 

## Prerequisites

The PVampliconFinder workflow is designed for the analysis of sequencing reads generated from **paired-end sequencing** of DNA amplified using degenerated primers targeting specifically the L1 sequence of papillomaviruses ([Chouhy *et al.*, 2010](https://www.ncbi.nlm.nih.gov/pubmed/19948351),[Forslund *et al.*, 1999](https://www.ncbi.nlm.nih.gov/pubmed/10501499),[Forslund *et al.*, 2003](https://www.ncbi.nlm.nih.gov/pubmed/12798239)).

### Dependencies & External Software

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


  Specify the test files location

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --param1    |            xx | ...... |
| --param2    |            xx | ...... |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --param3   |            xx | ...... |
| --param4    |            xx | ...... |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |
| --flag2    |      .... |


## Usage
  ```
  ...
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | output1    | ...... |
  | output2    | ...... |


## Detailed description (optional section)
...

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | contrib1*    |            xx | Developer to contact for support (link to specific gitter chatroom) |
  | contrib2    |            xx | Developer |
  | contrib3    |            xx | Tester |

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

## References (optional)

## FAQ (optional)




