# PVAmpliconFinder

**PVampliconFinder** is a data analysis workflow designed to rapidly identify and classify known and potentially new *papilliomaviridae* sequences from amplicon deep-sequencing with degenerated papillomavirus (PV) primers.

## Description

PVampliconFinder is based on alignment similarity metrics, but also consider molecular evolution time for an improved identification and taxonomic classification of novel PVs. The final output of the tool includes a list of fully characterized putatively new papillomaviriade sequences, as well as graphical representations of relative abundance of the virome sequence diversity in the tested samples. 

## Prerequisites

The PVampliconFinder workflow is designed for the analysis of sequencing reads generated from **paired-end sequencing** of DNA amplified using degenerated primers targeting specifically the L1 sequence of papillomaviruses ([Chouhy *et al.*, 2010](https://www.ncbi.nlm.nih.gov/pubmed/19948351),[Forslund *et al.*, 1999](https://www.ncbi.nlm.nih.gov/pubmed/10501499),[Forslund *et al.*, 2003](https://www.ncbi.nlm.nih.gov/pubmed/12798239)).

### Dependencies & External Software

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [VSEARCH](https://github.com/torognes/vsearch)
- [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [RaxML-EPA](https://cme.h-its.org/exelixis/web/software/epa/index.html)
- [PaPaRa](https://cme.h-its.org/exelixis/web/software/papara/index.html)
- [CAP3](http://seq.cs.iastate.edu/cap3.html)

### Databases

- [blastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/db/) : It is advised to use the script [update_blastddb.pl](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl). the required databases are :
  - nt
  - taxdb
- [ncbitax2lin](https://github.com/zyxue/ncbitax2lin) : A lineage file is available on the GitHub webpage, but a more recent one can be created manually following the instruction
- [PaVE](https://pave.niaid.nih.gov/#search/search_database/kw?dbNamespace=Genomes&includeNR=true&refCloneOnly=false&sort=Locus_ID&sortType=true&page=600&start=1&showTable=1&) : Select All > Download Fasta


```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

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

