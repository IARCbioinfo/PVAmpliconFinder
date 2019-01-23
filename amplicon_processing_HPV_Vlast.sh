#!/bin/bash

##############################################################
##	amplicon_processing_HPV_Vlast.sh						##
##	01/02/2017												##
##	Amplicon sequencing Illumina MiSeq						##
##	Alexis Robitaille : robitaillea@students.iarc.fr		##
##	IARC, LYON												##
##	Last modification = 08/10/2018							##
##	Version 1.1												##
##############################################################

##############
##	USAGE	##
##############
usage="$(basename "$0") [-h] [-t threads] [-b \"nt\" database] [-f info_file] [-p adapter R1] [-q adapter R2] -s fastq_files_suffix -d input_dir -o output_dir -- program to process amplicon-based NGS data
Version 1.0
The fastq filename to process must with the same suffix (option \"-s\").
The Read 1 filename must contain \"R1\" and the Read 2 filename must contain \"R2\" (the pair must otherwise have the same name).
See README for more information about $(basename "$0") usage.

where:
    -h  show this help text
    -s	suffix of fastq filename (ex : \"pool\" or \"sample\")
    -d  PATH to input fastq directory (.fastq | .fq | .zip | .tar.gz | .gz)
    -o	PATH to output directory
    -f	file containing pool information
    -b	\"nt\" blast database name (default \"nt\")
    -p	adapter sequence of read 1
    -q	adapter sequence of read 2
    -t	number of threads (default 2)
    "
    
##	Example
#time ./amplicon_processing_HPV_Vlast.sh -s pool -d /data/robitaillea/NGS1/fastq_files -o /data/robitaillea/NGS1/output -f /data/robitaillea/NGS1/infofile.txt -t 8

##	Get the parameters
while getopts ':hs:d:t:o:b:f:p:q:' option; do
  case "$option" in
    h) echo -e "$usage"
       exit
       ;;
    d) fastq_dir=$OPTARG
       ;;
    o) working_dir=$OPTARG
       ;;
    t) threads=$OPTARG
       ;; 
    s) suffix=$OPTARG
       ;;
    b) dbnt=$OPTARG
       ;; 
    f) info=$OPTARG
       ;; 
    p) primer=$OPTARG
       ;; 
    q) primer2=$OPTARG
       ;; 
    :) echo -e "$usage" >&2
	   printf "missing argument for -%s\n" "$OPTARG" >&2
       exit 1
       ;;
   \?) echo -e "$usage" >&2
	   printf "illegal option: -%s\n" "$OPTARG" >&2
       exit 1
       ;;
   (*) echo -e "$usage"
       exit
       ;;
  esac
done
shift $((OPTIND - 1))

primer="${primer%\\n}"
primer2="${primer2%\\n}"

##	Check if parameters are correct
if [ -z "$working_dir" ] || [ -z "$fastq_dir" ] || [ -z "$suffix" ]
then
   echo -e "$usage"
   exit
fi

if [ -z "$threads" ]
then
	threads=2
fi


if [ -z "$dbnt" ]
then
	dbnt="nt"
fi

##	Define absolute path if needed
if [[ ! "$working_dir" = /* ]] || [[ "$working_dir" = ^.+ ]]; then
	working_dir=${PWD}/${working_dir};
fi

if [[ ! "$fastq_dir" = /* ]] || [[ "$fastq_dir" = ^.+ ]]; then
	fastq_dir=${PWD}/${fastq_dir};
fi

if [ ! -z "$info" ]; then
	if [[ ! "$info" = /* ]] || [[ "$info" = ^.+ ]]; then
		info=${PWD}/${info};
	fi
fi

mkdir -p $working_dir;

##############################
##	Location of software	##
##############################
#~ fastqc="~/app/FastQC/fastqc";
#~ usearch="~/app/usearch9.2.64_i86linux32";
#~ trim_galore="~/app/trim_galore";
#~ blast="~/app/ncbi-blast-2.6.0+/bin/";
#~ multiQC="~/app/miniconda2/bin/multiqc";

##################
##	Log file	##
##################
##	Log file
logfile=$working_dir"/logfile.txt";
touch ${logfile} 2> /dev/null;

######################################
##	Reads Quality Control : FastQC	##
######################################
##	FastQC & MultiQC
fastqc_output_dir=${working_dir}"/fastqc_raw";

echo -e "##########################################\n##\tFastQC control of the raw reads\t\t##\n##########################################" >> ${logfile};
echo -e "##########################################\n##\tFastQC control of the raw reads\t##\n##########################################";

##	Check if QC on raw data already done
if [ ! -d ${working_dir}"/fastqc_raw" ]
then

	cd ${fastq_dir};
	
	##	Decompression if needed
	if [[ ! -z `find . -name "${suffix}*.zip"` ]]
	then
		echo -e "decompression zip";
		find . -name "${suffix}*.zip" | xargs --max-args=1 --max-procs=${threads} unzip 2>${logfile};
	elif [[ ! -z `find . -name "${suffix}*.tar.gz"` ]]
	then
		echo -e "decompression tar.gz";
		find . -name "${suffix}*.tar.gz" | xargs --max-args=1 --max-procs=${threads} tar xzvf 2>${logfile};
	elif [[ ! -z `find . -name "${suffix}*.gz"` ]]
	then
		echo -e "decompression gz";
		find . -name "${suffix}*.gz" | xargs --max-args=1 --max-procs=${threads} gunzip 2>${logfile};
	fi
	
	##	Output directory of fastqc result
	mkdir ${fastqc_output_dir} 2>> ${logfile};

	##	FastQC
	find . -name "${suffix}*.fastq" -o -name "${suffix}*R1*.fq" | xargs --max-args=1 --max-procs=${threads} fastqc -o ${fastqc_output_dir} -t ${threads} -q 2> ${logfile};
	
	cd ${fastqc_output_dir};

	##	MultiQC
	multiqc . 2>> ${logfile};
	
	ln -s ${fastqc_output_dir}"/multiqc_report.html" ${working_dir}"/multiQC_report_on_raw.html";
	
	echo "Done";
else
	echo -e "FastQC control of raw fastq files already done";
	echo -e "FastQC control of raw fastq files already done" >> ${logfile};
fi

##############################
##	Remove adapter sequence	##
##############################
##	Trim Galore & FastQC
trim_galore=${working_dir}"/trim_galore";

fastq_filtered=${working_dir}"/fastq_filtered";

cd ${working_dir};

echo -e "##########################################\n##\tRemove adapter sequence\t\t##\n##########################################" >> ${logfile};
echo -e "##########################################\n##\tRemove adapter sequence\t\t##\n##########################################";

if [ ! -d "fastq_filtered" ]
then

	mkdir ${trim_galore} 2>> ${logfile};

	cd ${fastq_dir};
	
	##	Trimming TrimGalore!
	if [ -z "$primer" ]
	then
		find . -name "${suffix}*R1*.fastq" -o -name "${suffix}*R1*.fq" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'r2="${0/R1/R2}"; echo Pair : ${0} $r2; trim_galore --fastqc --paired --retain_unpaired -o '${trim_galore}' ${0} $r2 &>>'${logfile}';';
	else
		find . -name "${suffix}*R1*.fastq" -o -name "${suffix}*R1*.fq" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'r2="${0/R1/R2}"; echo Pair : ${0} $r2; trim_galore --fastqc --paired --retain_unpaired -o '${trim_galore}' -a $primer -a2 $primer2 ${0} $r2 &>>'${logfile}';';
		#~ find . -name "${suffix}*R1*.fastq" -o -name "${suffix}*R1*.fq" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'r2="${0/R1/R2}"; echo Pair : ${0} $r2; trim_galore --fastqc --paired --retain_unpaired -o '${trim_galore}' -a '${primer}' -a2 '${primer2}' --clip_R1 9 --clip_R2 9 --three_prime_clip_R1 2 --three_prime_clip_R2 2 ${0} $r2 &>>'${logfile}';';
	fi
	
	cd ${trim_galore};
	
	##	MultiQC
	multiqc . 2>> $logfile
	
	ln -s ${trim_galore}"/multiqc_report.html" ${working_dir}"/multiQC_report_on_filtered.html";
	
	mkdir ${fastq_filtered} 2>> ${logfile};
	
	find . -name "${suffix}*val_1.fq" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'r2="${0/val_/unpaired_}"; cat $r2 >> ${0}; mv ${0} '${fastq_filtered}';'
	
	find . -name "${suffix}*val_2.fq" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'r2="${0/val_/unpaired_}"; cat $r2 >> ${0}; mv ${0} '${fastq_filtered}';'
	
	cd ${fastq_filtered};
	
	find . -name "${suffix}*R1*.fq" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'mv ${0} $(basename "${0/_val_1.fq/}").fq'
	find . -name "${suffix}*R2*.fq" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'mv ${0} $(basename "${0/_val_2.fq/}").fq'
	
	echo "Done";
else
	echo -e "Trim_galore already done";
	echo -e "Trim_galore already done" >> $logfile;
fi

##################################
##	Metagenomic steps : VSEARCH	##
##################################
##	Merging Pairs
##	Dereplication
##	Chimerics removal
##	Clustering

echo -e "##########################################\n##\tClustering step : VSEARCH\t##\n##########################################" >> $logfile;
echo -e "##########################################\n##\tClustering step : VSEARCH\t##\n##########################################";

cd ${working_dir};

outputdir=${working_dir}"/vsearch";

if [ ! -d ${working_dir}"/vsearch" ]
then
	
	mkdir ${outputdir} 2>> $logfile
	
	tmpdir=${working_dir}"/tmp";
	
	mkdir ${tmpdir} 2>> $logfile
	
	cd ${fastq_filtered};

	##	Merging Pairs
	echo -e "~~	MergePair	~~";
	find . -name "${suffix}*R1*.fq" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'r2="${0/R1/R2}"; [[ ${0} =~ \.*/*(${suffix}.*)[\._-]R1.* ]]; echo Pair : ${0} $r2 - Label : ${BASH_REMATCH[1]}; vsearch --quiet --fastq_mergepairs ${0} --reverse $r2 --minseqlength 30 --threads '${threads}' --fastq_allowmergestagger --label_suffix ${BASH_REMATCH[1]} --fastaout '${tmpdir}/'${BASH_REMATCH[1]}.fasta &>> '${logfile}';'
	
	cd ${tmpdir};
	##	Dereplication
	echo -e "~~	Dereplicate	~~";
	find . -name "${suffix}*.fasta" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'vsearch --quiet --derep_fulllength '${tmpdir}/'${0} --sizeout --threads '${threads}' --relabel_sha1 --fasta_width 0 --minuniquesize 2 --output '${tmpdir}'/$(basename "${0/.fasta/}")_lin_der.fasta --log $(basename "${0/.fasta/}").log  &>> '${logfile}';'
	##	Chimerics removal
	echo -e "~~	ChimericSeqRemoval	~~";
	find . -name "${suffix}*_lin_der.fasta" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'vsearch --quiet --uchime_denovo '${tmpdir}/'${0} --sizein --threads '${threads}' --relabel $(basename "${0/_lin_der.fasta/}") --sizeout --xsize --nonchimeras '${tmpdir}/'$(basename "${0/_lin_der.fasta/}")_no_chim.fasta --log $(basename "${0/_lin_der.fasta/}")_chimeria.log &>> '${logfile}';';
	##	Clustering
	echo -e "~~	Clustering	~~";
	find . -name "${suffix}*_no_chim.fasta" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'vsearch --quiet --cluster_size '${tmpdir}/'${0} --id 0.98 --threads '${threads}' --sizein --clusterout_id --clusterout_sort --sizeout --xsize --relabel $(basename "${0/_no_chim.fasta/}") --centroids '${outputdir}/'$(basename "${0/_no_chim.fasta/}").fasta --log $(basename "${0/_no_chim.fasta/}")_clustering.log &>> '${logfile}';';
	
	cd ${working_dir};

	mv ${tmpdir} ${working_dir}/log;
	echo "Done";
else
	echo -e "Vsearch already done";
	echo -e "Vsearch already done" >> $logfile;
fi			

######################################
##	Sequence identification : BLAST	##
######################################
##	Nucleotide-Nucleotide BLAST

echo -e "##########################################\n##\tSequence identification : BLAST\t##\n##########################################" >> $logfile;
echo -e "##########################################\n##\tSequence identification : BLAST\t##\n##########################################";

blastdir=${working_dir}"/blast_result";

if [ ! -d ${working_dir}"/blast_result" ]
then

	mkdir ${blastdir} 2>> $logfile;

	cd ${outputdir};
	
	##	BLASTN
	find . -name "${suffix}*.fasta" | xargs --max-args=1 --max-procs=${threads} -- bash -c 'name=$(basename "${0/.fasta/}"); echo $name; blastn -task megablast -use_index false -db '${dbnt}' -query ${0} -out '${blastdir}'/${name}.blast -evalue 1e-05 -max_target_seqs 1 -num_threads '${threads}' -outfmt "6 qseqid sseqid evalue bitscore length pident frames staxids sskingdoms sscinames scomnames sblastnames stitle qseq qstart qend";'
	
	cd ${blastdir};
	
	for f in *.blast; do
		sed -i '1iQueryID\tSubjectID\tevalue\tbitscore\tlength query\tperc id\tframes\ttaxid\tkingdom\tscientifique name\tcommon name\tblast name\ttitle\tseq query\tstartq\tstopq' ${f};
	done
	
	echo "Done";
else
	echo -e "Sequence identification already done";
	echo -e "Sequence identification already done" >> $logfile;	
fi

##########################
##	Advanced Analysis	##
##########################
echo -e "##########################################\n##\tAdvanced analysis\t\t##\n##########################################" >> $logfile;
echo -e "##########################################\n##\tAdvanced analysis\t\t##\n##########################################";

cd ${working_dir};
cd ..;

if [ ! -d ${working_dir}"/analysis_new" ]
then
	chmod +x $PWD/analysis_summary_HPV_Vlast.pl;
	#~ echo ${blastdir};
	#~ echo ${working_dir};
	#~ echo ${suffix};
	#~ echo ${outputdir};
	#~ echo ${threads};
	#~ echo ${info};
	perl $PWD/analysis_summary_HPV_Vlast.pl -i ${blastdir} -o ${working_dir} -s ${suffix} -d ${outputdir} -t ${threads} -f ${info};
else
	echo -e "Advanced analysis already done";
	echo -e "Advanced analysis already done" >> $logfile;	
fi
