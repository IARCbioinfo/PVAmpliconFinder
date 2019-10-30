#/usr/bin/perl -s

##############################################################
##	analysis_summary_HPV_Vlast.pl							##
##	04/09/2017												##
##	Amplicon sequencing Illumina MiSeq	- 	Step 2			##
##	Alexis Robitaille : robitaillea@students.iarc.fr		##
##	IARC, LYON												##
##	Last modification = 08/10/2018							##
##	Version 1.1												##
##############################################################

##################
##	Libraries	##
##################
use strict;
use warnings;
use Getopt::Std;
use Math::Round ':all';
use File::Basename; # my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
use List::MoreUtils qw(uniq);
use List::MoreUtils qw(first_index);
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::Index::Fasta;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SearchIO;
use Bio::Seq;
use Text::CSV;
use Cwd;
use Bio::TreeIO;

##############
##	Options	##
##############
my $inputdirfasta="";
my $inputdirblast="";
my $outputdir="";
my $suffix="";
my $infos="";
my $threads="";
my %opts = ();
getopts( 'd:o:i:s:f:t:h', \%opts );
if ( defined( $opts{d} ) ) {
	$inputdirfasta=$opts{d};
}
else {
	print_usage();
	exit;
}
if ( defined( $opts{i} ) ) {
	$inputdirblast=$opts{i};
}
else {
	print_usage();
	exit;
}
if ( defined( $opts{o} ) ) {
	$outputdir=$opts{o};
}
else {
	print_usage();
	exit;
}
if ( defined( $opts{s} ) ) {
	$suffix=$opts{s};
}
else {
	print_usage();
	exit;
}
if ( defined( $opts{f} ) ) {
	$infos=$opts{f};
}
if ( defined( $opts{t} ) ) {
	$threads=$opts{t};
}
else {
	print_usage();
	exit;
}
if ( $opts{h} ) {
	print_usage();
	exit;
}

##############################
##	Command line example	##
##############################
#	time perl analysis_summary_HPV_Vlast.pl -i output/blast_result -o output/results -s pool -d output/vsearch -t 8 -f infofile.txt;


my $whereiam=getcwd;		##	where the analysis is done

my $dirname = dirname(__FILE__);	# where the script is located

if($inputdirfasta!~/^\//){
	$inputdirfasta=$whereiam."/".$inputdirfasta;
}

if($inputdirblast!~/^\//){
	$inputdirblast=$whereiam."/".$inputdirblast;
}

if($outputdir!~/^\//){
	$outputdir=$whereiam."/".$outputdir;
}
##############
##	Usage	##
##############
sub print_usage {
	printf "
analysis_summary.pl [-h] [-f infos_file] -i input_dir_blast -o output_dir -s suffix -d input_dir_fasta -t threads_number

	-i	input directory		{blastdir}
	-o	output directory	{working_dir}
	-s	suffix				{suffix}
	-d	fasta files directory	{outputdir}
	-f	infos file			{info}
	-t	number of threads	{threads}
	-h	show this help\n"
}


#Virus	Species Name	--> 	FILES (easy update but have to be given with the tool) OR HASHTABLE IN THE SCRIPT (less clean/hard to update but no extra file)??
#~ my $hpvinfos="/data/robitaillea/THESIS/HPV_type.txt";

##	CETTE INFO PEUT ETRE RECUPEREE DU BLASTN

#BLAST result for putative new
my $output_new=$outputdir."/analysis_new";
mkdir $output_new;
#BLAST result for putative known
my $output_known=$outputdir."/analysis_known";
mkdir $output_known;


##	Information of the type of primer used for each pool
my %primer=();	##	Keys=pools; Values = primer type
my %tissu=();	##	Keys=pools; Values = tissu type
my @pools=();	##	Table of pool (uniqu ID)
my @primers=();	##	Table uniq primer
my %cptprimer=();	##	Count number of time each primer is present
my %poolsid=();		##	Give an id to each pool ordered by "sort @pools"	Keys=pools; Values = id
my $boolean_infofile="true";

##########################################
##	Check the integrity of info file	##
##########################################
if($infos ne ""){
	my $cptlineinfo=`wc -l $infos | cut -d' ' -f 1`;chomp($cptlineinfo);
	my $cptlineblast=`ls -l $inputdirblast/"$suffix"*.blast | wc -l`;chomp($cptlineblast);
	my $cptlinefasta=`ls -l $inputdirfasta/"$suffix"*.fasta | wc -l`;chomp($cptlinefasta);
	#~ print $cptlineinfo."\t".$cptlineblast."\t".$cptlinefasta."\n";
	if($cptlineinfo=~ /^\d+$/){$cptlineinfo--;}else{print "Error : Info file badly formatted\n";exit;}
	if(($cptlinefasta eq $cptlineblast) && ($cptlinefasta eq $cptlineinfo)){
		my $res=readinfofile($infos);
		%primer=%{@$res[0]};		#  a ordoné par sort keys %primer;
		%tissu=%{@$res[1]};
		@pools=sort keys %primer;
		@primers=uniq sort values %primer;
	}
	elsif(($cptlineinfo ne $cptlinefasta) && ($cptlinefasta eq $cptlineblast)){
		print "Error : Info file badly formatted\n Did you forget the header (ID\tprimer\tpool)\n";exit;
	}
	elsif($cptlinefasta ne $cptlineblast){
		print "Error : Missing files of the blast result and/or fasta input\n";exit;
	}
}
else{
	print "No info file provided\n";
	
	my $names=`ls $inputdirblast | grep "^$suffix" | sed -e 's/.blast//g'`; chomp($names);
	my @tabnames=split("\n", $names);
	`touch $outputdir/infos_tmp.csv`;
	my $text="ID\tprimer\ttissue\n";
	foreach my $n (@tabnames){
		$text.=$n;
		$text.="\tNA\tNA\n";
		$primer{$n}="NA";
		$tissu{$n}="NA";
		push(@primers,"NA");
	}
	@pools=sort keys %primer;
	`echo "$text" > $outputdir/infos_tmp.csv`;
	$boolean_infofile="false";
}

##	Information of the type of primer used for each pool
#~ my %primer=();	##	Keys=pools; Values = primer type
#~ my %tissu=();	##	Keys=pools; Values = tissu type
#~ my @pools=();	##	Table of pool (uniqu ID)
#~ my @primers=();	##	Table uniq primer


##############################################
##	Check if empty file blast ou sequence	##						---> IF ISSUE WITH PRIMER IT'S HERE
##############################################						--> REMOVE PRIMER OF POOL IF EMPTY AND THE ONLY ONE IN %primer (no other pool has this primer)
foreach my $file (`ls -l $inputdirblast/*.blast`)
{
	chomp($file);
	my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
	#~ print $inputdirblast."/".$filename.$suffix."\n";
	if(-z $inputdirblast."/".$filename.$suffix){
		print "$filename is empty and thus will not be proceed by the program\n";
		my $tmp=delete $tissu{$filename};
		my $tmp2=delete $primer{$filename};
		
		my $index = 0;
		$index++ until $pools[$index] eq $filename;
		splice(@pools, $index, 1);
	}
}
#~ exit;

##	Initialize
foreach my $p (@primers){
	$cptprimer{$p}=0;
}

##	Compute the count
foreach my $p (values %primer){
	$cptprimer{$p}++;
}
my $cptpoolsid=1;
foreach my $p (sort @pools){
	$poolsid{$p}=$cptpoolsid;
	$cptpoolsid++;
}

sub readinfofile{
	my $file=shift;
	my %primer=();
	my %tissu=();
	
	chomp($file);
	my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
	open(F1, $file) or die "$!: $file\n";
	
	my $header = <F1>;
	my $tissucol=-1;
	my $primercol=-1;
	
	if((!($header =~ /tissue/i)) or (!($header =~ /primer/i))){
		print "Error : Infofile header badly formatted\nPlease follow the example : ID	primer	tissue\nSee README\n";exit;
	}
	else{
		my @head=split("\t",$header);
		foreach (my $col=0; $col<=$#head; $col++){
			if($head[$col]=~/tissu/i){
				$tissucol=$col;
			}
			elsif($head[$col]=~/primer/i){
				$primercol=$col;
			}
		}
	}
	while(<F1>)
	{
		my $line=$_;
		$_      =~ s/[\r\n]+$//;
		my @tab = split("\t", $_);
		foreach (my $col=0; $col<=$#tab; $col++){
			$primer{$tab[0]}=$tab[$primercol];
			$tissu{$tab[0]}=$tab[$tissucol];
		}
	}
	my @return=();
	push (@return, \%primer);
	push (@return, \%tissu);
	return (\@return);
}

##########################################
##	Define classification from PAVE		##
##########################################
##	NEED TO DOWNLOAD 2 FILES FROM PAVE (THEY ARN'T FORMATED THE SAME!!!!!)
	#	Human Ref Clone table : https://pave.niaid.nih.gov/#explore/reference_genomes/human_genomes --> export
	#	Animal Ref Clone Table : https://pave.niaid.nih.gov/#explore/reference_genomes/animal_genomes --> export
##	Both file has to be save as .csv or txt file (separator tabular)

#~ awk -v OFS="\t" -F"\t" '{print $1, $9, $3}' download_animal_RefClone_7896321e.csv
#~ awk -v OFS="\t" -F"\t" '{print $1, $5, $2}' download_human_RefClone_09156ea6.csv

#~	sudo apt-get install gnumeric

#~ cat <(sed -e '1d' download_animal_RefClone_*.csv | awk -v OFS="\t" -F"\t" '{print $1, $9, $3}') <(sed -e '1d' download_human_RefClone_*.csv | awk -v OFS="\t" -F"\t" '{print $1, $5, $2}') > pave_table.txt

my $table="$dirname/raxml/pave_table.txt";

my %hacctohpv=();
my %hacctotax=();
my %hhpvtotax=();
open(T, $table) or die "$! : $table\n";	#HPV1	V01116	Mupapillomavirus 1
while(<T>){
	chomp($_);
	my @tab=split(/\t/,$_);
	if((defined($tab[1])) && $tab[1]!~/^\s*$/){
		$hacctohpv{$tab[1]}=$tab[0];		#	BPV7	DQ217793	Dyoxipapillomavirus 1
	}
	if((defined($tab[2])) && $tab[2]!~/^\s*$/){
		$hacctotax{$tab[1]}=$tab[2];		#	BPV7	DQ217793	Dyoxipapillomavirus 1
	}
	if((defined($tab[0])) && $tab[0]!~/^\s*$/){
		$hhpvtotax{$tab[0]}=$tab[2];		#	BPV7	DQ217793	Dyoxipapillomavirus 1
	}
}
close(T);


my %htype;#number of reads by pool and organism category - key = pool & organism type
my %hhpv;#number of reads by pool and HPV category - key = tissue - VIRUS type
my %reads;#number of reads by pool
my %hprimer;#number of reads by primer and organism category - key = pool - organism type

my %hcountnew;#title new VIRUS by primer - key = primer | values = table of new virusname
my %hcountknow;#title known VIRUS by primer - key = primer | values = table of known virusname
my %hcountother;#title other by primer - key = primer | values = table of other seq

my %hnew;#title new VIRUS by pool - key = pool | values = number of uniqu new virusname
my %hknown;#title known VIRUS by pool - key = pool | values = number of uniqu known virusname
my %hother;#title other by pool - key = pool | values = number of uniqu other name

my %htarget;#{$pool}{known/new}{virusname}+=$reads
my %htargetcluster;#{$pool}{known/new}{virusname}+=$nbcluster

##For the NEW DiversityByTissu
my %hnewtable;	#{$tissu}{$virusfamily}{$virusgenus}{$virusspecies}{$virusname}=@(nbcluster;nbreads_pool1,nbreads_pool2,...,nbreads_pool8)

my %readstissue;	#{$tissue}=nb VIRUS reads
my %clusttissue;	#{$tissue}=nb VIRUS cluster

my %virusname2taxid;

##	Initialisation
my $j=1;
for my $pool (sort @pools){
	$htype{$pool}{'newVIRUS'}=0;
	$htype{$pool}{'VIRUS'}=0;
	$htype{$pool}{'human'}=0;
	$htype{$pool}{'bacteria'}=0;
	$htype{$pool}{'other'}=0;
	$reads{$pool}=0;
}
for my $p (values %primer){
	#~ print $p." here\n";
	$hprimer{$p}{'newVIRUS'}=0;
	$hprimer{$p}{'VIRUS'}=0;
	$hprimer{$p}{'human'}=0;
	$hprimer{$p}{'bacteria'}=0;
	$hprimer{$p}{'other'}=0;
}


for my $p (values %primer){
	@{$hcountnew{$p}}=();
	$hknown{$p}=0;
	$hnew{$p}=0;
	@{$hcountknow{$p}}=();
}



loadblastresult($inputdirblast);

sub loadblastresult{
	my $dir=shift;
	
	my %consider = map { $_ => 1 } @pools;
	
	######################################
	##	LOAD THE DATA FROM BLAST OUTPUT	##
	######################################
	foreach my $file (`ls $dir/*.blast`)
	{
		
		my $pool="";
		my $tissu="";
		chomp($file);
		my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
		
		##	Not consider if empty file
		if(!(exists($consider{$filename}))){
			next;
		}
		
		print $filename."\n";
		my $sampleName = $filename;
		
		open(F1, $file) or die "$!: $file\n";
		
		##	Write a file with only the HPV more than 10% dissimilarity sequence
		open(IN, ">".$output_new."/".$filename.".csv") or die "$!: $filename\n";
		##	Write a file with all the known HPV
		open(KNO, ">".$output_known."/".$filename.".csv") or die "$!: $filename\n";
		
		my $header = <F1>;
		print IN $header;
		print KNO $header;
		
		while(<F1>)
		{
			my $line=$_;
			$_      =~ s/[\r\n]+$//;
			my @tab = split("\t", $_);
			my $reads=0;
			
			
			
			my $qid=$tab[0];				##	Unique ID for each cluster 		| Pool1_S1_L001_R1_001.fastq1;size=228160;
			my $sid=$tab[1];				##	GI number of matching sequence	| gi|373158195|gb|JN231328.1|
			my $eval=$tab[2];				##	Evalue							| 3.15e-44
			my $bitscore=$tab[3];			##	Bit score						| 187
			my $qlen=$tab[4];				##	Length of the query				| 101
			my $percid=$tab[5];				##	Percentage of identity			| 100.000
			my $frame=$tab[6];				##	Frame query/subject				| 1/-1
			my $taxid=$tab[7];				##	Taxon Id						| 1034800
			my $kingdom=$tab[8];			##	Kingdom classification			| Viruses
			my $sname=$tab[9];				##	Scientifique name				| uncultured Papillomavirus
			my $cname=$tab[10];				##	Common name						| uncultured Papillomavirus
			my $bname=$tab[11];				##	blast name						| viruses
			my $title=$tab[12];				##	blast title						| Uncultured Papillomavirus contig01 putative L2 and putative L1 genes, complete cds
			my $qseq=$tab[13];				##	sequence
			my $qstart=$tab[14];			##	start align position query
			my $qend=$tab[15];				##	end align position query
			
			#print $title."\n";
			
			##	Parsing of the ID
			foreach my $p (sort @pools){
				if($qid=~/^($p)\d+;size=(\d+);*$/){
					$pool=$1;
					$reads=$2;
				}
				elsif($qid=~/^(Undertermined)\d+;size=(\d+);*$/){
					$pool=$1;
					$reads=$2;
				}
			}
			
			if($pool eq ""){
				print "Error in the format of the sequence ID : ".$qid."\n";
				exit;
			}
			
			$reads{$pool}+=$reads;
			
			##	Viruses
			if($kingdom=~/Viruses/i && $sname=~/papilloma/i){
				my $virusname="";
				$virusname=$title;
				
				if($taxid=~/;/){
					my @tabtax=split(';',$taxid);
					$taxid=$tabtax[0];
				}
				$virusname2taxid{$virusname}=$taxid;		#taxid LIPyV 1965344
				
				$readstissue{$tissu{$pool}}+=$reads;
				$clusttissue{$tissu{$pool}}+=1;
				
				#~ if($percid>=90.0){
				if($percid>=90.0){							#Test thershold 5%
					
					$htype{$pool}{'VIRUS'}+=$reads;
					$hprimer{$primer{$pool}}{'VIRUS'}+=$reads;
					
					
					push(@{$hcountknow{$pool}},$title);
					
					print KNO $_."\n";	
					
						
					if($virusname=~/(\S+)(\.|,|:|;)$/){
						chop($virusname);
					}
					
					$hhpv{$pool}{$virusname}+=$reads;
					
					
					$htarget{$pool}{'known'}{$virusname}+=$reads;
					$htargetcluster{$pool}{'known'}{$virusname}+=1;
				}
				##For the new one
				else{
					$htype{$pool}{'newVIRUS'}+=$reads;	#OK
					$hprimer{$primer{$pool}}{'newVIRUS'}+=$reads;
					
					push(@{$hcountnew{$pool}},$title);
					
					print IN $_."\n";

					if($virusname=~/(\S+)(\.|,|:|;)$/){
						chop($virusname);
					}
					
					$hhpv{$pool}{$virusname}+=$reads;
					
					$htarget{$pool}{'new'}{$virusname}+=$reads;
					$htargetcluster{$pool}{'new'}{$virusname}+=1;

				}		
				
			}
			##	Human
			elsif($kingdom=~/Eukaryota/i && $sname=~/Homo sapiens/i){
				$htype{$pool}{'human'}+=$reads;
				$hprimer{$primer{$pool}}{'human'}+=$reads;
				push(@{$hcountother{$pool}},$title);
			}
			##	Bacteria
			elsif($kingdom=~/Bacteria/i){
				$htype{$pool}{'bacteria'}+=$reads;
				$hprimer{$primer{$pool}}{'bacteria'}+=$reads;
				push(@{$hcountother{$pool}},$title);
			}
			##	Other
			else{
				$htype{$pool}{'other'}+=$reads;
				$hprimer{$primer{$pool}}{'other'}+=$reads;
				push(@{$hcountother{$pool}},$title);
			}
			
		}
		close(KNO);
		close(IN);
		close(F1);
	}
	return;
}


##	compute the number of uniqu VIRUS type (new/known) in the different primer type
##############################################
##	Count of new / uniqu known HPV by pool	##
##############################################
my @tmp=();
my %hVIRUStot=();
foreach my $pool (sort @pools){
	#KNOWN
	@tmp=();
	@tmp=uniq @{$hcountknow{$pool}};
	push(@{$hknown{$pool}},@tmp);
	#~ $hHPVtot{"known"}+=scalar @tmp;
	push(@{$hVIRUStot{"known"}},@tmp);
	
	#NEW
	@tmp=();
	@tmp=uniq @{$hcountnew{$pool}};
	push(@{$hnew{$pool}},@tmp);
	#~ $hHPVtot{"new"}+=scalar @tmp;
	push(@{$hVIRUStot{"new"}},@tmp);

	#OTHER
	@tmp=();
	@tmp=uniq @{$hcountother{$pool}};
	push(@{$hother{$pool}},@tmp);
	#~ $hHPVtot{"other"}+=scalar @tmp;
	push(@{$hVIRUStot{"other"}},@tmp);
}


##Reduce at uniqu knwown and other in all the experiment --> not for new because 2 new with same blast results in 2 differents pools are not necessary the same new HPV
my %hVIRUStot_uniq=();
foreach my $k (keys %hVIRUStot){
	if($k ne "new"){
		$hVIRUStot_uniq{$k}=uniq @{$hVIRUStot{$k}};
	}
}


##################################################
##	Write table_summary_MegaBlast_results.csv	##
##################################################
open(OUT,">".$outputdir."/table_summary_MegaBlast_results.csv") or die "$!";
print OUT "MegaBlast results\n";

##################################################
##	Table 1 - Sequencing statistics per samples	##
##################################################
print OUT "Table 1 - Sequencing statistics per samples\n";
if($boolean_infofile eq "true"){
	print OUT "Pool\tPrimer\tTissue\tAllreads\tVIRUS_reads\tpVIRUS_reads\tOther_reads\tpOther_reads\tVIRUSnew_reads\tpVIRUSnew_onTot\tpVIRUSnew_onVIRUSknown\tNumber of new VIRUS\tNumber of known VIRUS\n";
}else{
	print OUT "Pool\tAllreads\tVIRUS_reads\tpVIRUS_reads\tOther_reads\tpOther_reads\tVIRUSnew_reads\tpVIRUSnew_onTot\tpVIRUSnew_onVIRUSknown\tNumber of new VIRUS\tNumber of known VIRUS\n";	
}
for my $p (sort keys %htype){
	
	my $othersreads=$htype{$p}{'human'}+$htype{$p}{'bacteria'}+$htype{$p}{'other'};
	#~ print $othersreads." other\n";
	my $allreads=$othersreads+$htype{$p}{'VIRUS'}+$htype{$p}{'newVIRUS'};
	#~ print $allreads." all\n";
	my $virusreads=$htype{$p}{'VIRUS'}+$htype{$p}{'newVIRUS'};
	my $fracnew;
	if($virusreads == 0){
		$fracnew=0;
	}
	else{
		$fracnew=nearest(.001,($htype{$p}{'newVIRUS'}/$virusreads)*100)
	}
	#~ print $virusreads." virus\n";
	
	if($allreads!=$reads{$p}){
		print "Error\n";
		exit;
	}
	if($boolean_infofile eq "true"){
		print OUT $p."\t".$primer{$p}."\t".$tissu{$p}."\t".$reads{$p}."\t".$virusreads."\t".nearest(.001,($virusreads/$allreads)*100)."\t".$othersreads."\t".nearest(.001,($othersreads/$allreads)*100)."\t".$htype{$p}{'newVIRUS'}."\t".nearest(.001,($htype{$p}{'newVIRUS'}/$allreads)*100)."\t".$fracnew."\t".scalar(@{$hnew{$p}})."\t".scalar(@{$hknown{$p}})."\n";
	}else{
		print OUT $p."\t".$reads{$p}."\t".$virusreads."\t".nearest(.001,($virusreads/$allreads)*100)."\t".$othersreads."\t".nearest(.001,($othersreads/$allreads)*100)."\t".$htype{$p}{'newVIRUS'}."\t".nearest(.001,($htype{$p}{'newVIRUS'}/$allreads)*100)."\t".$fracnew."\t".scalar(@{$hnew{$p}})."\t".scalar(@{$hknown{$p}})."\n";
	}	
}

##################################################
##	Table 2 - Sequencing statistics per primers	##
##################################################
if($boolean_infofile eq "true"){
	print OUT "\n\n";
	print OUT "Table 2 - Sequencing statistics per primers\n";
	print OUT "Primer\tAllreads\tVIRUS_reads\tpVIRUS_reads\tOther_reads\tpOther_reads\tVIRUSnew_reads\tpVIRUSnew_onTot\tpVIRUSnew_onVIRUSknown\n";
}


foreach my $primer (sort keys %hprimer){
	if($cptprimer{$primer}>1){
		$hprimer{$primer}{'other'}=nearest(1,($hprimer{$primer}{'other'}/$cptprimer{$primer}));
		$hprimer{$primer}{'bacteria'}=nearest(1,($hprimer{$primer}{'bacteria'}/$cptprimer{$primer}));
		$hprimer{$primer}{'human'}=nearest(1,($hprimer{$primer}{'human'}/$cptprimer{$primer}));
		$hprimer{$primer}{'newVIRUS'}=nearest(1,($hprimer{$primer}{'newVIRUS'}/$cptprimer{$primer}));
		$hprimer{$primer}{'VIRUS'}=nearest(1,($hprimer{$primer}{'VIRUS'}/$cptprimer{$primer}));
		$hknown{$primer}=nearest(1,($hknown{$primer}/$cptprimer{$primer}));
		$hnew{$primer}=nearest(1,($hnew{$primer}/$cptprimer{$primer}));
	}
	my $othersreads=$hprimer{$primer}{'human'}+$hprimer{$primer}{'bacteria'}+$hprimer{$primer}{'other'};
	my $allreads=$othersreads+$hprimer{$primer}{'VIRUS'}+$hprimer{$primer}{'newVIRUS'};
	my $virusreads=$hprimer{$primer}{'VIRUS'}+$hprimer{$primer}{'newVIRUS'};
	my $fracnew;
	
	if($virusreads == 0){
		$fracnew=0;
	}
	else{
		$fracnew=nearest(.001,($hprimer{$primer}{'newVIRUS'}/$virusreads)*100)
	}
	
	if($allreads==0){
		print "Unexpected error in the count of reads for ".$primer." primer\n";
		exit;
	}
	if($boolean_infofile eq "true"){
		print OUT $primer."\t".$allreads."\t".$virusreads."\t".nearest(.001,($virusreads/$allreads)*100)."\t".$othersreads."\t".nearest(.001,($othersreads/$allreads)*100)."\t".$hprimer{$primer}{'newVIRUS'}."\t".nearest(.001,($hprimer{$primer}{'newVIRUS'}/$allreads)*100)."\t".$fracnew."\n";
	}
	else{
		print OUT "Mean\t".$allreads."\t".$virusreads."\t".nearest(.001,($virusreads/$allreads)*100)."\t".$othersreads."\t".nearest(.001,($othersreads/$allreads)*100)."\t".$hprimer{$primer}{'newVIRUS'}."\t".nearest(.001,($hprimer{$primer}{'newVIRUS'}/$allreads)*100)."\t".$fracnew."\t".nearest(1,(scalar(@{$hVIRUStot{"new"}})/scalar(@pools)))."\t".nearest(1,(scalar(@{$hVIRUStot{"known"}})/scalar(@pools)))."\n";
	}
}


print OUT "\n\n";

my %hprimnew=();
my %hprimknown=();
my %hprimother=();

foreach my $pool (@pools){
	$hprimnew{$primer{$pool}}+=scalar(@{$hnew{$pool}});
	$hprimknown{$primer{$pool}}+=scalar(@{$hknown{$pool}});
	$hprimother{$primer{$pool}}+=scalar(@{$hother{$pool}});
}
##############################################################
##	Table 3 - Sequencing statistics per primers & tissue	##
##############################################################
print OUT "Table 3 - Sequencing statistics per primers & tissue\n";
if($boolean_infofile eq "true"){
	print OUT "Primer\tTissue\tOther\tVIRUSknown\tVIRUSnew\n";
}
else{print OUT "Pool\tOther\tVIRUSknown\tVIRUSnew\n";
	
}
foreach my $pool (sort @pools){
	if($boolean_infofile eq "true"){
		print OUT $primer{$pool}."\t".$tissu{$pool}."\t".scalar(@{$hother{$pool}})."\t".scalar(@{$hknown{$pool}})."\t".scalar(@{$hnew{$pool}})."\n";
	}
	else{
		print OUT $pool."\t".scalar(@{$hother{$pool}})."\t".scalar(@{$hknown{$pool}})."\t".scalar(@{$hnew{$pool}})."\n";	
	}
}

######################################################
##	Table 4 - Species identification per primers	##
######################################################
if($boolean_infofile eq "true"){
	print OUT "\n\n";
	print OUT "Table 4 - Species identification per primers\n";
	print OUT "Primer\tOther\tVIRUSknown\tVIRUSnew\n";
	foreach my $primer (sort @primers){
		print OUT $primer."\t".$hprimother{$primer}."\t".$hprimknown{$primer}."\t".$hprimnew{$primer}."\n";
	}
}
print OUT "Total species\t".scalar(@{$hVIRUStot{"other"}})."\t".scalar(@{$hVIRUStot{"known"}})."\t".scalar(@{$hVIRUStot{"new"}})."\n";
print OUT "Total unique species\t".scalar(uniq @{$hVIRUStot{"other"}})."\t".scalar(uniq @{$hVIRUStot{"known"}})."\t".scalar(@{$hVIRUStot{"new"}})."\n";
print OUT "\n\n";

my %known=();	#{pool}{family}{genus}=nb reads
my %new=();	#{pool}{family}{genus}=nb reads

my @famnew=();	#all family of new HPVs
my @famknown=();#all family of known HPVs

my %hgennew=();		#all genus of new HPVs per family : keys= family/values=tab of genus
my %hgenknown=();		#all genus of kown HPVs per family : keys= family/values=tab of genus



##############################
##	GET LINEAGE INFORMATION	##
##############################
##	A partir du taxid recupérer la classification (dans la boucle initial), faire par taxid
##	https://github.com/zyxue/ncbitax2lin

#cd ~/app/ncbitax2lin
#source activate venv/
#make

##	grep -i "virus" lineages* > lineages-virus.csv

##	gunzip lineages*.csv.gz

##	Utilisation de la structure précédemment construite :	$polyoname2taxid{$virusname}=$taxid;

##########################
##	!!	ATTENTION	!!	##
##########################
##	Utilisation d'un chemin absolu!!!!

my %taxid2gen=();
my %taxid2spe=();
my %taxid2fam=();

my $lineage="$dirname/databases/lineages/lineagesVirus.csv";		

my $csv = Text::CSV->new ( { binary => 1 } ) or die "Cannot use CSV: ".Text::CSV->error_diag ();

open my $io, "<", $lineage or die "$!: $lineage";

while(my $row = $csv->getline ($io)){
	chomp($row);
	my @tab = @$row;
	chomp($tab[0]);	#taxid
	chomp($tab[5]);	#family
	chomp($tab[6]);	#genus
	chomp($tab[7]);	#species
	chomp($tab[14]);	#complement info
	
	if($tab[5]!~/^\s*$/){
		$taxid2fam{$tab[0]}=$tab[5];
	}
	else{
		$taxid2fam{$tab[0]}="Unclassified";
	}
	
	if($tab[6]!~/^\s*$/){
		$taxid2gen{$tab[0]}=$tab[6];
	}
	else{
		$taxid2gen{$tab[0]}="Unclassified";
	}
	if($tab[7]!~/^\s*$/ && $tab[7]!~/^Human papillomavirus$/i){
		$taxid2spe{$tab[0]}=$tab[7];
	}
	else{
		$taxid2spe{$tab[0]}="Unclassified";
	}
}

close($io);

my %hoverall=();	# {Family}{Genus}{Species}{virusname}{tissu}=nb reads

foreach my $pool (sort keys %htarget){		#HERE sorted par sample name 
							
	foreach my $virusname (keys %{$htarget{$pool}{'known'}}){
		if((!(defined($virusname2taxid{$virusname}))) or $virusname2taxid{$virusname} eq ""){
			print $pool." error line 652\n";
			#~ exit;
		}
		
		if(!(defined($taxid2gen{$virusname2taxid{$virusname}}))){
			print $virusname."\n";
			print $virusname2taxid{$virusname}."\n";
			exit;
		}
		
		my $genus=$taxid2gen{$virusname2taxid{$virusname}};
		my $species=$taxid2spe{$virusname2taxid{$virusname}};
		my $family=$taxid2fam{$virusname2taxid{$virusname}};
		
		$known{$pool}{$family}{$genus}+=$htarget{$pool}{'known'}{$virusname};
		
		##Table diversity by site (tissue)
		if(defined($hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname})){
			##	Indice 0 is kept for number of cluster
			$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[0]+=$htargetcluster{$pool}{'known'}{$virusname};	#nbcluster

			##	Use of the pool id attribuated at the begining, from 1 to nb of pools --> indice 0 is save for nb of cluster (see line above)
			$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[$poolsid{$pool}]+=$htarget{$pool}{'known'}{$virusname};					#ICI sorted poolsid{pool} donc id
		}
		else{
			##	Initialization
			@${$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}=();
			$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[0]=$htargetcluster{$pool}{'known'}{$virusname};	#nbcluster
			
			##	Nombre de fois retrouvé dans le pool X (see %poolsid hashtable for poolsname to corresponding id (begin at 1 not at 0!! --> see above))
			$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[$poolsid{$pool}]=$htarget{$pool}{'known'}{$virusname};
		}
		
		$hoverall{$family}{$genus}{$species}{$virusname}{$tissu{$pool}}+=$htarget{$pool}{'known'}{$virusname};
		
	}
		
	foreach my $virusname (keys %{$htarget{$pool}{'new'}}){
		if((!(defined($virusname2taxid{$virusname}))) or $virusname2taxid{$virusname} eq ""){
			print $pool." error line 633\n";
			#~ exit;
		}
		my $genus=$taxid2gen{$virusname2taxid{$virusname}};
		my $species=$taxid2spe{$virusname2taxid{$virusname}};
		my $family=$taxid2fam{$virusname2taxid{$virusname}};
		
		$new{$pool}{$family}{$genus}+=$htarget{$pool}{'new'}{$virusname};
		
		##Table diversity by site (tissue)
		if(defined($hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname})){
			##	Indice 0 is kept for number of cluster
			$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[0]+=$htargetcluster{$pool}{'new'}{$virusname};	#nbcluster

			##	Use of the pool id attribuated at the begining, from 1 to nb of pools --> indice 0 is save for nb of cluster (see line above)
			$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[$poolsid{$pool}]+=$htarget{$pool}{'new'}{$virusname};
		}
		else{
			##	Initialization
			@${$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}=();
			$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[0]=$htargetcluster{$pool}{'new'}{$virusname};	#nbcluster
			
			##	Nombre de fois retrouvé dans le pool X (see %poolsid hashtable for poolsname to corresponding id (begin at 1 not at 0!! --> see above))
			$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[$poolsid{$pool}]=$htarget{$pool}{'new'}{$virusname};
		}
		
		$hoverall{$family}{$genus}{$species}{$virusname}{$tissu{$pool}}+=$htarget{$pool}{'new'}{$virusname};
		
	}
	
	foreach my $fam (keys %{$known{$pool}}){
		if (! grep {$_ eq $fam} @famknown){
			push(@famknown,$fam);
		}
		foreach my $gen (keys %{$known{$pool}{$fam}}){
			if (! grep {$_ eq $gen} @{$hgenknown{$fam}}){
				push(@{$hgenknown{$fam}},$gen);
			}
		}
	}
	foreach my $fam (keys %{$new{$pool}}){
		if (! grep {$_ eq $fam} @famnew){
			push(@famnew,$fam);
		}
		foreach my $gen (keys %{$new{$pool}{$fam}}){
			if (! grep {$_ eq $gen} @{$hgennew{$fam}}){
				push(@{$hgennew{$fam}},$gen);
			}
		}
	}
}

##########################################
##	Table 5 - KNOWN virus family level	##
##########################################
#~ print $#catnew."\t".$#catknown."\n";
##	WRITE TABLE KNOWN VIRUS CATEGORY
print OUT "Table 5 - KNOWN virus family level\n";
print OUT "Pool";
foreach my $fam (sort @famknown){     
	print OUT "\t".$fam;
}
print OUT "\n";
foreach my $pool (sort keys %known){
	print OUT $pool;
	foreach my $fam (sort @famknown){
		if(defined($known{$pool}{$fam})){
			my $sum=0;
			foreach my $gen (sort keys %{$known{$pool}{$fam}}){
				$sum+=$known{$pool}{$fam}{$gen};
			}
			print OUT "\t".$sum;
		}
		else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
##########################################
##	Table 6 - KNOWN virus genus level	##
##########################################
print OUT "\nTable 6 - KNOWN virus genus level\n";
print OUT "Family\t";
foreach my $fam (sort @famknown){
	my $nbr = $#{$hgenknown{$fam}};
	print OUT $fam;
	for(my $i=0; $i<=$nbr; $i++){
		print OUT "\t";
	}
	#~ print OUT "\t".$fam;
}
print OUT "\nPool\\Genus\t";
foreach my $fam (sort @famknown){
	foreach my $gen (sort @{$hgenknown{$fam}}){
		print OUT $gen."\t";
	}
}
print OUT "\n";
foreach my $pool (sort keys %known){
	print OUT $pool;
	foreach my $fam (sort @famknown){
		foreach my $gen (sort @{$hgenknown{$fam}}){
			if(defined($known{$pool}{$fam}{$gen})){
				print OUT "\t".$known{$pool}{$fam}{$gen};
			}
			else{
				print OUT "\t0";
			}
		}
	}
	print OUT "\n";
}

##########################################
##	Table 7 - NEW virus family level	##
##########################################
print OUT "\nTable 7 - NEW virus family level\n";
print OUT "Pool";
foreach my $fam (sort @famnew){     
	print OUT "\t".$fam;
}
print OUT "\n";
foreach my $pool (sort keys %new){
	print OUT $pool;
	foreach my $fam (sort @famnew){
		if(defined($new{$pool}{$fam})){
			my $sum=0;
			foreach my $gen (sort keys %{$new{$pool}{$fam}}){
				$sum+=$new{$pool}{$fam}{$gen};
			}
			print OUT "\t".$sum;
		}
		else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}

######################################
##	Table 8 - NEW virus genus level	##
######################################
print OUT "\nTable 8 - NEW virus genus level\n";
print OUT "Family\t";
foreach my $fam (sort @famnew){
	my $nbr = $#{$hgennew{$fam}};
	print OUT $fam;
	for(my $i=0; $i<=$nbr; $i++){
		print OUT "\t";
	}
	#~ print OUT "\t".$fam;
}
print OUT "\nPool\\Genus\t";
foreach my $fam (sort @famnew){
	foreach my $gen (sort @{$hgennew{$fam}}){
		print OUT $gen."\t";
	}
}
print OUT "\n";
foreach my $pool (sort keys %new){
	print OUT $pool;
	foreach my $fam (sort @famnew){
		foreach my $gen (sort @{$hgennew{$fam}}){
			if(defined($new{$pool}{$fam}{$gen})){
				print OUT "\t".$new{$pool}{$fam}{$gen};
			}
			else{
				print OUT "\t0";
			}
		}
	}
	print OUT "\n";
}

close(OUT);
######################################################
##	END Write table_summary_MegaBlast_results.csv	##
######################################################

my @listprimer=uniq sort values %primer;


##################################################
##	WRITE TABLE HPV DIVERSITY BY TISSUE TYPE	##
##################################################
##				!!	NEW VERSION	!!				##
##################################################
##	Independant tables per tissue type
##	ADD KRONA files

my $krona=$outputdir."/KRONA_MegaBlast";
mkdir $krona;

my $ff=0;
my $fg=0;
my $fs=0;
my $fn=0;
foreach my $t (sort keys %hnewtable){		#For each tissue
	
	my $kronaname="";
	
	my @primertmp_no_uniq=();
	my @poolsidtmp=();
	foreach my $p (sort @pools){			# Pour chaque nom ordered
		if($tissu{$p} eq $t){				# Si le tissu du pollname correspond a celui procedé
			push(@primertmp_no_uniq,$primer{$p});		#le primer de ce pool										$p = NGSDEEP_DL2A1_S11_L001 --> $primer{$p} = 11
			push(@poolsidtmp,$poolsid{$p});		###	Both table are in the same order	$p = NGSDEEP_DL2A1_S11_L001 --> $poolsid{$p} = 	
		}
	}
	my @primertmp=();
	@primertmp=uniq @primertmp_no_uniq;
	
	if($boolean_infofile eq "true"){
		open(OUT,">".$outputdir."/diversityByTissu_".$t."_MegaBlast.csv") or die "$!";
		print OUT "TISSUE\tFamily\tGenus\tSpecies\tRelated\tNb Cluster (%)\t".join("\t",sort @primertmp)."\tPool\n";		#vérifier ordre primertmp
		print OUT $t."\t";
		$kronaname="krona_".$t."_MegaBlast";
		open(KRO,">".$krona."/krona_".$t."_MegaBlast.txt") or die "$!";
	}
	else{
		open(OUT,">".$outputdir."/OverallDiversity_MegaBlast.csv") or die "$!";
		print OUT "Family\tFamily\tGenus\tSpecies\tRelated\tNb Cluster (%)\t".join("\t",sort @primertmp)."\tPool\n";		#vérifier ordre primertmp
		print OUT "Virus\t";
		
		$kronaname="krona_ByTissue_MegaBlast";
		open(KRO,">".$krona."/krona_ByTissue_MegaBlast") or die "$!";
	}
	foreach my $f (sort keys %{$hnewtable{$t}}){
		if($ff==0){
			print OUT $f."\t";
			$ff=1;
		}
		else{
			print OUT "\t".$f."\t";
		}
		foreach my $g (sort keys %{$hnewtable{$t}{$f}}){		#For each genus	$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}
			if($fg==0){
				print OUT $g."\t";
				$fg=1;
			}
			else{
				print OUT "\t\t".$g."\t";
			}
			foreach my $s (sort keys %{$hnewtable{$t}{$f}{$g}}){		#For each species
				if($fs==0){
					print OUT $s."\t";
					$fs=1;
				}
				else{
					print OUT "\t\t\t".$s."\t";
				}
				foreach my $n (sort keys %{$hnewtable{$t}{$f}{$g}{$s}}){		#For each virusname
					my @tmp_pool=();
					my %tmp_primer=();
					if($fn==0){
						print OUT $n."\t";
						$fn=1;
					}
					else{
						print OUT "\t\t\t\t".$n."\t";
					}
					my $tabtmp=$@{$hnewtable{$t}{$f}{$g}{$s}{$n}};				##Get the corresponding table :	0 --> cluster number ; then initialize only for pool number having this HPV
																			##	0	3	4	5
					
					my $cluster=@$tabtmp[0]." (".nearest(.01,((@$tabtmp[0]/$clusttissue{$t})*100)).")";		#Get cluster 

					for my $i (@poolsidtmp){	#	3	4	5	6	7	8	9	10	11	
						if((defined(@$tabtmp[$i])) && @$tabtmp>0){		#	3	4	5
							my %hpoolsidrev= reverse %poolsid;			#	poolid vers poolname	
							my $p=$hpoolsidrev{$i};						#	$p=NGSDEEP_DL2A1_S11_L001	--> $i=3
							$tmp_primer{$primer{$p}}+=@$tabtmp[$i];		#	$primer{$p} = nb reads of pool number X having this primer		
							push(@tmp_pool,$p);			#	get pool name - $p = NGSDEEP_DL2A1_S11_L001
						}
					}
					my $tmp_towrite="";
					my $total_krona=0;
					foreach my $prim (sort @primertmp){#	all primer concerned
						if(exists($tmp_primer{$prim})){		# tmp_primer de NGSDEEP_DL2A1_S11_L001 --> nb reads of pool number (3)	
							$tmp_towrite.=$tmp_primer{$prim}."\t";		##on ecrit
							$total_krona+=$tmp_primer{$prim};
						}
						else{
							$tmp_towrite.="-\t";
						}
					}
					chomp($tmp_towrite);
					print OUT $cluster."\t".$tmp_towrite.join(",",@tmp_pool)."\n";
					
					print KRO $total_krona."\t".$t."\t".$g."\t".$s."\n";
				}
				$fn=0;
			}
			$fs=0;
		}
		$fg=0;
	}
	$ff=0;
	
	close(OUT);
	close(KRO);
	
	`ktImportText $krona/$kronaname.txt -o $krona/$kronaname.html`;
	
}



##########################################
##	WRITE TABLE HPV DIVERSITY OVERALL	##
##########################################
$ff=0;
$fg=0;
$fs=0;
$fn=0;
my @sorttissu=();
foreach my $t (sort keys %hnewtable){		#For each tissu (sorted)
	push(@sorttissu,$t);
}

open(OUTALL,">".$outputdir."/DiversityByTissu_MegaBlast.csv") or die "$!";
print OUTALL "TISSUE\tFamily\tGenus\tSpecies\tRelated\t".join("\t",@sorttissu)."\tPool\n";		#vérifier ordre primertmp
print OUTALL "ALL\t";

open(KRO,">".$krona."/krona_OverAll_MegaBlast.txt") or die "$!";

#	time perl analysis_summary_polyo.pl -i /home/robitaillea/ICB/DATA_DEEP/NGS_POLYO/blast_result -o /home/robitaillea/ICB/DATA_DEEP/NGS_POLYO/newtable -s NGSDEEP_ -d /home/robitaillea/ICB/DATA_DEEP/NGS_POLYO/vsearch -t 8 -f /home/robitaillea/ICB/DATA_DEEP/NGS_POLYO/infofile.csv

foreach my $f (sort keys %hoverall){		#For each family
	if($ff==0){
		print OUTALL $f."\t";
		$ff=1;
	}
	else{
		print OUTALL "\t".$f."\t";
	}
	foreach my $g (sort keys %{$hoverall{$f}}){		#For each genus (sorted)
		#~ print OUTALL $g."\t";
		if($fg==0){
			print OUTALL $g."\t";
			$fg=1;
		}
		else{
			print OUTALL "\t\t".$g."\t";
		}
		foreach my $s (sort keys %{$hoverall{$f}{$g}}){	#For each species
			if($fs==0){
				print OUTALL $s."\t";
				$fs=1;
			}
			else{
				print OUTALL "\t\t\t".$s."\t";
			}
			foreach my $n (sort keys %{$hoverall{$f}{$g}{$s}}){		#For each polyoname
				if($fn==0){
					print OUTALL $n."\t";
					$fn=1;
				}
				else{
					print OUTALL "\t\t\t\t".$n."\t";
				}
				
				my @pool_list=();
				
				my $total_krona=0;
				foreach my $t (@sorttissu){
					my $nbread="-";
					if(defined($hoverall{$f}{$g}{$s}{$n}{$t})){
						$nbread=$hoverall{$f}{$g}{$s}{$n}{$t};
						$total_krona+=$hoverall{$f}{$g}{$s}{$n}{$t};
						
							## To display the pool concerned	$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[$poolsid{$pool}]=$htarget{$pool}{'new'}{$virusname};
						my $tabtmp;
						if(defined($hnewtable{$t}{$f}{$g}{$s}{$n})){
							$tabtmp=$@{$hnewtable{$t}{$f}{$g}{$s}{$n}};
						}
						#~ else{
							#~ print $hoverall{$f}{$g}{$s}{$n}{$t}."\n";
							#~ print $t."\n";
							#~ print $f."\n";
							#~ print $g."\n";
							#~ print $s."\n";
							#~ print $n."\n";
							#~ print "Error\n";
							#~ exit;
						#~ }
						for(my $i=1; $i<=$#{$tabtmp}; $i++){		# 0 is kept for the number of cluster
							if((defined(@$tabtmp[$i])) && @$tabtmp>0){			#	3	4	5
								my %hpoolsidrev= reverse %poolsid;				#	poolid vers poolname	
								my $p=$hpoolsidrev{$i};							
								push(@pool_list,$p);
							}
						}
					
					}
					print OUTALL $nbread."\t";	
					
				}
				print OUTALL join(',',@pool_list);
				print OUTALL "\n";
				
				print KRO $total_krona."\tALL\t".$g."\t".$s."\n";
			}
			$fn=0;
		}
		$fs=0;
	}
	$fg=0;
}
$ff=0;
close(KRO);
close(OUTALL);

`ktImportText $krona/krona_OverAll_MegaBlast.txt -o $krona/krona_OverAll_MegaBlast.html`;

##########
##	NEW	##
##########
print "Concat NEW VIRUS sequence\n";
my @res=concat_sequence($output_new,$inputdirfasta);
my %hdata=%{$res[0]};
my %hseq=%{$res[1]};

##############
##	KNOWN	##
##############
print "Concat KNOWN VIRUS sequence\n";
my @res_known=concat_sequence($output_known,$inputdirfasta);
my %hdata_known=%{$res_known[0]};
my %hseq_known=%{$res_known[1]};


##################################
##	Concat_sequence function	##
##################################
##	This function take as input a directory path of blast result file and identify the common cluster in each file, try to concatene the common information for this clusters.
##	INPUT :	- Directory Path of blast results (csv format)
##			- Directory Path of fasta sequence (fasta format)	!!	The name of the file as to be the same as the name of the output blast file	!!
##	OUTPUT:	- \$hdata{title}{pool}=@infos([0]=uniqu_id,[1]=percentage_id,[2]=frame,[3]=reads,[4]=ginumber,[5]=pool,[6]=tissue,[7]=primer,[8]=seq_len,[9]=query_seq,[10]=query_id)
##			- \$hseq{title}{pool}=@seq(id1|seq1,id2|seq2,id3|seq3,...)

sub concat_sequence{
	my $dir=shift;
	my $inputdirfasta=shift;
	
	my $pool="";
		## TREAT THE NEW FILE WRITE WITH ONLY INFORMATIVE READS
	#	hash table {title}{pool}=table of infos
	my %hdatafun=();

	#	hash table {title}{pool}=table of seq if several seq for this title
	my %hseqfun=();

	foreach my $file (`ls $dir/*.csv`)
	{
		chomp($file);
		my ($filename, $directories, $suffix) = fileparse($file, qr/\.[^.]*/);
		print $filename."\n";
		my $sampleName = $filename;
		
		my $db       = Bio::DB::Fasta->new($inputdirfasta."/".$filename.".fasta", -maxopen=>1000, -debug);
		
		my @ids = $db->get_all_primary_ids;
			
		#~ print "Le nombre de sequence dans $file est $#ids\n";
		
		open(F1, $file) or die "$!: $file\n";
		
		my $header = <F1>;
		
		while(<F1>)
		{
			my $line=$_;
			$_      =~ s/[\r\n]+$//;
			my @tab = split("\t", $_);
			
			my @infos=();
			
			my $reads=0;
			my $qid=$tab[0];
			my $sid=$tab[1];
			my $eval=$tab[2];
			my $bitscore=$tab[3];
			my $qlen=$tab[4];
			my $percid=$tab[5];
			my $frame=$tab[6];
			my $taxid=$tab[7];
			my $kingdom=$tab[8];
			my $sname=$tab[9];
			my $cname=$tab[10];
			my $bname=$tab[11];
			my $title=$tab[12];
			my $qseq=$tab[13];
			my $qstart=$tab[14];			##	start align position query
			my $qend=$tab[15];				##	end align position query
			
			chomp($qid);
			

			my $seq     = $db->get_Seq_by_id("$qid");

			if(!(defined($seq))){
				print "Pour chaque sequence dans $file \n La sequence "; 
				print $qid." est pas trouver dans $inputdirfasta\t $filename\n";
				exit;
			}
			my $seqstr  = $seq->seq;
			chomp($seqstr);
			chomp($qid);
			
			##	Parsing of the ID
			if($qid=~/^($sampleName)\S+;size=(\d+);*$/){
				$pool=$1;
				$reads=$2;
			}
			else{
				print "Error in the format of the sequence ID : ".$qid."\n";
				exit;
			}
			my $tissue=$tissu{$pool};
			
			
			if(defined($hdatafun{$title}{$pool})){	#Déja vu dans le pool
				if($frame ne ${$hdatafun{$title}{$pool}}[2]){	#Si les frame sont opposés
					my $revcomp=reverse_complement_IUPAC($seqstr);
					if(reverse_complement_IUPAC($seqstr) eq ${$hdatafun{$title}{$pool}}[9]){	#Et si la nouvelle sequence est l'exact complémentaire de celle déjà connue
						##	Concatenation du nombre de reads dans la table de hash
						${$hdatafun{$title}{$pool}}[3]+=$reads;
					}
					elsif($revcomp=~/${$hdatafun{$title}{$pool}}[9]/){		#Alors nouvelle sequence + grande que l'ancienne --> MAJ des infos par rapport à la plus grande sequence
						${$hdatafun{$title}{$pool}}[1]=$percid;
						${$hdatafun{$title}{$pool}}[2]=$frame;
						${$hdatafun{$title}{$pool}}[3]+=$reads;
						${$hdatafun{$title}{$pool}}[8]=$qlen;
						${$hdatafun{$title}{$pool}}[9]=$seqstr;
						
						${$hdatafun{$title}{$pool}}[11]=$qstart;
						${$hdatafun{$title}{$pool}}[12]=$qend;
						
					}
					elsif(${$hdatafun{$title}{$pool}}[9]=~/$revcomp/){		#Alors ancienne sequence + grande que la nouvelle --> garde ancienne infos (ajout du nombre de read)
						${$hdatafun{$title}{$pool}}[3]+=$reads;
					}
					else{	##Try to do a consensus with cap3
						#~ print "ici\t".$#{$hseqfun{$title}{$pool}}."\n";
						
						${$hdatafun{$title}{$pool}}[3]+=$reads;
						
						if(!($hseqfun{$title}{$pool})){					# Une seule sequence deja presente  pour ce cluster
							push(@{$hseqfun{$title}{$pool}},${$hdatafun{$title}{$pool}}[10]."|".${$hdatafun{$title}{$pool}}[9]); #On lui donne la séquence d'avant que si la table est vide
							push(@{$hseqfun{$title}{$pool}},$qid."|".$seqstr);				# Et on lui ajoute la nouvelle
						}
						else{												# D'autres sequence deja connu pour ce clsuter
							my $similar="false";
							foreach my $prevseq (@{$hseqfun{$title}{$pool}}){			## Sinon pour chaque sequence de ce cluster déja presente
								if(reverse_complement_IUPAC($seqstr) eq $prevseq){		## Exact complementaire d'une sequence de ce clsuter
									#~ ${$hdatafun{$title}{$pool}}[3]+=$reads;
									$similar="true";
									last;
								}
								elsif($revcomp=~/$prevseq/){							#Alors nouvelle sequence identique et + grande que une previous sequence de ce cluster
									$similar="true";
									last;
								}
								elsif($prevseq=~/$revcomp/){							#Alors ancienne sequence identique et + grande que la nouvelle
									#~ ${$hdatafun{$title}{$pool}}[3]+=$reads;
									$similar="true";
									last;
								}
							}
							if($similar eq "false"){					## Auncune similarité trouver avec es autres sequences du cluster
								push(@{$hseqfun{$title}{$pool}},$qid."|".$seqstr);			# Ajout comme new sequence du cluster
							}
						}
											#Et la nouvelle
						#~ print reverse_complement_IUPAC($qseq)."\n".${$hdatafun{$title}{$pool}}[9]."\n\n";
						#~ `cd $inputdir`;
						#~ `mkdir tmp`;
						#~ open(">".$inputdir."/tmp/
					}
				}
				else{	#Les frames sont les mêmes
					if($seqstr eq ${$hdatafun{$title}{$pool}}[9]){	#Et si la nouvelle sequence est l'exact complémentaire de celle déjà connue
						##	Concatenation du nombre de reads dans la table de hash
						${$hdatafun{$title}{$pool}}[3]+=$reads;
					}
					elsif($seqstr=~/${$hdatafun{$title}{$pool}}[9]/){		#Alors nouvelle sequence + grande que l'ancienne --> MAJ des infos
						${$hdatafun{$title}{$pool}}[1]=$percid;
						${$hdatafun{$title}{$pool}}[2]=$frame;
						${$hdatafun{$title}{$pool}}[3]+=$reads;
						${$hdatafun{$title}{$pool}}[8]=$qlen;
						${$hdatafun{$title}{$pool}}[9]=$seqstr;
						
						${$hdatafun{$title}{$pool}}[11]=$qstart;
						${$hdatafun{$title}{$pool}}[12]=$qend;
					}
					elsif(${$hdatafun{$title}{$pool}}[9]=~/$seqstr/){		#Alors ancienne sequence + grande que la nouvelle --> garde ancienne infos (ajout du nombre de read)
						${$hdatafun{$title}{$pool}}[3]+=$reads;
					}
					else{	##Try to do a consensus with cap3
						#~ print "ici\t".$#{$hseqfun{$title}{$pool}}."\n";
						${$hdatafun{$title}{$pool}}[3]+=$reads;
						
						if(!($hseqfun{$title}{$pool})){						# Une seule sequence deja presente  pour ce cluster
							push(@{$hseqfun{$title}{$pool}},${$hdatafun{$title}{$pool}}[10]."|".${$hdatafun{$title}{$pool}}[9]); #On lui donne la séquence d'avant que si la table est vide
							push(@{$hseqfun{$title}{$pool}},$qid."|".$seqstr);				# Et on lui ajoute la nouvelle
						}
						else{												# D'autres sequence deja connu pour ce clsuter
							my $similar="false";
							foreach my $prevseq (@{$hseqfun{$title}{$pool}}){			## Sinon pour chaque sequence de ce cluster déja presente
								if($seqstr eq $prevseq){		## Exact complementaire d'une sequence de ce clsuter
									#~ ${$hdatafun{$title}{$pool}}[3]+=$reads;
									$similar="true";
									last;
								}
								elsif($seqstr=~/$prevseq/){							#Alors nouvelle sequence identique et + grande que une previous sequence de ce cluster
									$similar="true";
									last;
								}
								elsif($prevseq=~/$seqstr/){							#Alors ancienne sequence identique et + grande que la nouvelle
									#~ ${$hdatafun{$title}{$pool}}[3]+=$reads;
									$similar="true";
									last;
								}
							}
							if($similar eq "false"){					## Auncune similarité trouver avec es autres sequences du cluster
								push(@{$hseqfun{$title}{$pool}},$qid."|".$seqstr);			# Ajout comme new sequence du cluster
							}
						}
					}
				}
			}
			else{
				#infos[0]=uniqu id, infos[1]=percentage id,infos[2]=frame, infos[3]=reads, infos[4]=ginumber, infos[5]=pool, infos[6]=tissue, infos[7]=primer,infos[8]=seq len, infos[9]=query seq, infos[10]=query id, infos[11]=qstart, infos[12]=qend
				
				#	List info of the line to keep in the final output
				my @infos=();
				push(@infos,($j."VIRUSput",$percid,$frame,$reads,$sid,$pool,$tissue,$primer{$pool},$qlen,$seqstr,$qid,$qstart,$qend));
				$hdatafun{$title}{$pool}=\@infos;
				$j++;
			}
			
		}
	}
	my @return=();
	push(@return,\%hdatafun);
	push(@return,\%hseqfun);
	return(@return);
}

######################################################################
##	HERE WE HAVE TO DOWNLOAD, OR GIVE WITH THE TOOL, THE PAVE DB	##
## fasta ---> Be carfeul of taxonomy ID that are all 0, blast hate it
#	solution :
#				sed -i 's/>gi|/>IDgi|/g' PaVE.fasta


#	makeblastdb -in PaVE.fasta -input_type fasta -hash_index -dbtype nucl (-parse_seqids)

##################################
##	BLASTN DATABASE DEFINITION	##
##################################
my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
   -db_data => $dirname.'/databases/PaVE/PaVE.fasta',
   -create => 1,
   -alphabet=>'dna'
		);
		
##########
##	NEW	##
##########
my $tmp=$output_new."/tmp";
mkdir $tmp;

chdir $output_new;
print "Create contigs NEW VIRUS\n";
my @aecrire=contig_build(\%hdata,\%hseq,$inputdirfasta,$tmp);

##############
##	KNOWN	##
##############
my $tmp3=$output_known."/tmp";
mkdir $tmp3;

chdir $output_known;
print "Create contigs KNOWN VIRUS\n";
my @aecrire3=contig_build(\%hdata_known,\%hseq_known,$inputdirfasta,$tmp3);


##############################
##	contig_build function	##
##############################
##	This function take the output of Concat_sequence function (%hdata and %hseq) and try to construct contig (cap3) based on the different sequence that get the same megablast match.
##	Then it launch a blastn of this contig against a PaVE database to get the closest HPV complete genome, and then gat the classification of this best HPV match.
##	INPUT :	- %hdata{title}{pool}=@infos([0]=uniqu_id,[1]=percentage_id,[2]=frame,[3]=reads,[4]=ginumber,[5]=pool,[6]=tissue,[7]=primer,[8]=seq_len,[9]=query_seq,[10]=query_id)
##			- %hseq{title}{pool}=@seq(id1|seq1,id2|seq2,id3|seq3,...)
##			- Directory Path of fasta sequence (fasta format)	!!	The name of the file as to be the same as the name of the output blast file	!!
##			- Directory path for tmp file of blast and cap3
##	OUTPUT:	- Table of information to write in a file 
sub contig_build{
	my $hdatafun=shift;
	my $hseqfun=shift;
	my $inputdirfasta=shift;
	my $tmpfun=shift;

	my @aecrirefun=();

	my $k=0;
	foreach my $pool (sort @pools){
		##########################
		##	Pour chaque POOL	##
		##########################
		
		##Préparation du blastn
		my $p=$pool;																##Transformer Pool en pool car nom de fichier uniquement en minuscule
				 
		foreach my $title (keys %$hdatafun){
			##########################
			##	Pour chaque cluster	##
			##########################		
			
			if(defined($$hseqfun{$title}{$pool})){	
				
				##############################################
				##	Si plusieurs sequences dans ce cluster	##
				##############################################
				
				##	Put the sequence in a tmpfun fasta file
				my $s=1;
				open(TMP, ">".$tmpfun."/$k$k$k.fasta") or die "$!: $k$k$k.fasta\n";
				foreach my $seq (@{$$hseqfun{$title}{$pool}}){
					my @tabtmp=split(/\|/,$seq);
					my $idtmp=$tabtmp[0];
					my $l=$tabtmp[1];
					print TMP ">".$idtmp."\n".$l."\n";
					#~ print ">".$idtmp."\n".$l."\n";
					$s++;
				}
				close(TMP);
				
				##############################
				##	Try to create a contig	##
				##############################
				
				`touch $tmpfun/$k$k$k.fasta.cap.contigs`;
				my $nb=1;

				my $first=-1;
				#~ my $toto="";

				## 
				
				
				while($nb!=0){
					if($first>0){
						`cat $tmpfun/$k$k$k.fasta.cap.contigs $tmpfun/$k$k$k.fasta.cap.singlets > $tmpfun/$k$k$k.fasta`;
					}
					#~ print "cap3 $tmpfun/$k$k$k.fasta > $tmpfun/$k$k$k.res\n";
					my $toto = `cap3 $tmpfun/$k$k$k.fasta > $tmpfun/$k$k$k.res`;
					$nb = `wc -l $tmpfun/$k$k$k.fasta.cap.contigs | cut -d ' ' -f 1`;
					chomp($nb);
					$first=1;
				}

				##	On calcule l'abondance et le percentage d'identité
				
				my $perc_dis=nearest(.01,(100.00-(${$$hdatafun{$title}{$pool}}[1])));
				my $abundance=nearest(.0001,((${$$hdatafun{$title}{$pool}}[3]/$reads{$pool})*100));
				
				##	On prépare les tables ou seront stocké les résutlats a écrire (longueur des séquences et les séquences)
				my $seqf="";
				my @tseqtmpfun=();
				my @lentemp=();
				
				##	Initialisation des varaible récupérant l'HPV le plus proche ainsi que sa classification
				my @hpvclosest=();
				my @classification=();
				my $hpvclosest="";
				my $classification="";
				my @start_stop_len_hit=();
				
				##################################
				##	Plusieurs contigs construit	##
				##################################
				##	Open contigs file			
				my $dbtmpfun	=	Bio::SeqIO->new( -format => 'Fasta', -file => $tmpfun."/$k$k$k.fasta");

				while((my $seqtmpfun = $dbtmpfun->next_seq())) {
					
					my $seqstrtmpfun  = $seqtmpfun->seq;
					
					#~ print "BEFORE ALL\n";
					#~ if(${$$hdatafun{$title}{$pool}}[0]=~/^51VIRUSput/){
						#~ print $seqstrtmpfun."\n";
					#~ }
					
					my $len		   = $seqtmpfun->length;
					
					#~ print $len." ici\t";
					##	Get length and sequence
					
					push(@lentemp,$len);
					
					my $id  = $seqtmpfun->display_id();
						
					my $title_tmp=$title;
					
					$title_tmp=~ s/\//_/g;
					$title_tmp=~ s/,/_/g;
					$title_tmp=~ s/\s/_/g;
					
					my $out="$tmpfun/$pool.$title_tmp.$id.blast"; #HERE TO CORRECT
					
					my $in2=HPVblastn($seqtmpfun,$out);


					my $result=$in2->next_result;									
					my $var=$result->query_name();
					my $hit=$result->next_hit();
					
					my $frame_query="";
					my $frame_hit="";
						
					if(defined $hit){
						my $desc = $hit->description();
						my $hit_name = $hit->name();
						#~ print $desc." la\n";		##	By modifying structure of naming for fasta file, allow to do the blast, but not to get gi, so have to take it from name or description
						#~ print $hit_name." ici\n";
						
						my $classification_tmp="";
						
						$classification_tmp=$hhpvtotax{$hit_name};
						
						if(!(defined($classification_tmp)) or ($classification_tmp eq "")){
							$classification_tmp="Unclassified";
						}
						#~ print $hpv."\n";
						#~ print $classification_tmp."\n";
						
						#~ my @tab_name=split(/\|/,$hit_name);
						#~ my $gi = $tab_name[1];
						#~ chomp($gi);

						#~ my $taxid="";

						#~ $taxid=$gi_tax{$gi};
						#~ chomp($taxid);
						#~ my $classification_tmp="";
						#~ if(defined($taxid2spe{$taxid})){
							#~ $classification_tmp=$taxid2spe{$taxid};
						#~ }
						#~ elsif(defined($taxid2gen{$taxid})){
							#~ $classification_tmp=$taxid2gen{$taxid};
						#~ }
						#~ else{
							#~ $classification_tmp="Unclassified";
						#~ }
						
						
						my $percent_id_moy="NA";
						my $hsp=$hit->next_hsp();
						my @percent_id=();
						my $compte=0;
						my @start_stop_len_hsp=();

						if(defined $hsp){
							do{
								$compte++;
								push @percent_id, $hsp->percent_identity;
								
								my $queryaln=$hsp->start."-".$hsp->end."(".$hsp->length('query').")";
								push @start_stop_len_hsp, $queryaln;
								
								
								$frame_query=$hsp->query->strand;
								$frame_hit=$hsp->hit->strand;
								
								#~ print "HSP here : $compte\tFrame query : $frame_query\tFrame hit : $frame_hit\n";

										
							}while(my $hsp=$hit->next_hsp());
						}
						#~ if($compte>1){
							#~ print $title." LA\n";
						#~ }
						
						@start_stop_len_hsp=uniq(@start_stop_len_hsp);
						
						my $sum=0;
						map { $sum += $_ } @percent_id;
						$percent_id_moy=nearest(.01,(($sum/$compte)));
						push(@start_stop_len_hit,join(";",@start_stop_len_hsp));
						
						push(@hpvclosest,$hit_name."($percent_id_moy%)");
						push(@classification,$classification_tmp);
						
						if($frame_query!=$frame_hit){
							push(@tseqtmpfun,reverse_complement_IUPAC($seqstrtmpfun));	## ICI frame a changer en fct de blastn results
						}
						else{
							push(@tseqtmpfun,$seqstrtmpfun);	## ICI frame a changer en fct de blastn results
						}
						
					}				##	Si pas de hit
					else{
						push(@hpvclosest,"NA");
						push(@classification,"NA");
						push(@start_stop_len_hit,"NA");
						
						push(@tseqtmpfun,$seqstrtmpfun);
						#~ print "Pas de hit\n";
						#~ print ${$$hdatafun{$title}{$pool}}[0]."\n";
						#~ print $seqstrtmpfun."\n\n";
					}
					
					#~ if(${$$hdatafun{$title}{$pool}}[0]=~/^122VIRUSput/){
						#~ print $k.$k.$k."\n";
						#~ print $seqstrtmpfun."\n";
					#~ }
				}
				##	Fin traitement toutes les sequences du clusters
				
				
				##	Concatene information
				#~ foreach my $lt (@lentemp){
					#~ print $lt." here \n";
				#~ }
				######################################	HERE	######################################	
				##	Here rank having longuest sequence first
				my $len="";
				my $infoalign="";
				
				#~ if($#hpvclosest != $#classification){
					#~ print ${$$hdatafun{$title}{$pool}}[0]."\n";
					#~ print $k.$k.$k."fasta";
					#~ exit;
				#~ }
				
				#~ my @lenfun=@lentemp;
				
				#~ if(${$$hdatafun{$title}{$pool}}[0]=~/^51VIRUSput/){
					#~ my $tmp=0;
					#~ foreach my $lt (@lentemp){
						#~ print "Avant function : ".$lt."\n";
						#~ print $tseqtmpfun[$tmp]."\n";
						#~ $tmp++;
					#~ }
				#~ }
				
				
				#~ if(\@lenfun eq \@lentemp){
					#~ print "SAME REFERENCE\n";exit;
				#~ }
				
				if($#lentemp > 0){		## Plus d'un contigs/singlets au moins un match			 && (@start_stop_len_hsp != 0)
					my $res=longuestFirst(\@lentemp,\@hpvclosest,\@classification,\@start_stop_len_hit,\@tseqtmpfun);		# Give a copy to avoid replacement @len
					$len=$$res[0];
					$hpvclosest=$$res[1];
					$classification=$$res[2];
					$infoalign=$$res[3];
					$seqf=$$res[4];
				}
				else{				##	Un seul contigs/singlets au moins un match
					$len=$lentemp[0];
					$hpvclosest=$hpvclosest[0];
					$classification=$classification[0];
					$infoalign=$start_stop_len_hit[0];
					$seqf=$tseqtmpfun[0];
				}
				
				#~ print "\n".$len."\n";
				#~ else{												##	pas de match
					#~ print "error\n";exit;
					#~ $len=join("/",@lentemp);
					#~ $hpvclosest=join("/",@hpvclosest);
					#~ $classification=join("/",@classification);
				#~ }
				
				#~ if(${$$hdatafun{$title}{$pool}}[0]=~/^51VIRUSput/){
					#~ print "\nApès function : ".$len."\n";
					#~ print $seqf."\n";
					#~ print $k.$k.$k."\n";
					#~ exit;
				#~ }
				
				
				#~ my $len=join("/",@lentemp);
				#~ $hpvclosest=join("/",@hpvclosest);
				#~ $classification=join("/",@classification);
				
				#~ my $infoalign=join("/",@start_stop_len_hit);
				
				#~ #print $hpvclosest."par ici \n";
				
				#~ $seqf=join("\t",@tseqtmpfun);
				
				if($boolean_infofile eq "true"){
					#~ push(@aecrirefun,${$$hdatafun{$title}{$pool}}[0]."\t".$perc_dis."\t".$abundance."\t".${$$hdatafun{$title}{$pool}}[3]."\t".${$$hdatafun{$title}{$pool}}[4]."\t".$pool."\t".$title."\t".${$$hdatafun{$title}{$pool}}[6]."\t".${$$hdatafun{$title}{$pool}}[7]."\t".join("/",@lentemp)."\t".$infoalign."\t".$hpvclosest."\t".$classification."\t".$seqf);
					push(@aecrirefun,${$$hdatafun{$title}{$pool}}[0]."\t".$perc_dis."\t".$abundance."\t".${$$hdatafun{$title}{$pool}}[3]."\t".${$$hdatafun{$title}{$pool}}[4]."\t".${$$hdatafun{$title}{$pool}}[11]."-".${$$hdatafun{$title}{$pool}}[12]."\t".$title."\t".$pool."\t".${$$hdatafun{$title}{$pool}}[6]."\t".${$$hdatafun{$title}{$pool}}[7]."\t".$len."\t".$infoalign."\t".$hpvclosest."\t".$classification."\t".$seqf);
				
				}
				else{
					push(@aecrirefun,${$$hdatafun{$title}{$pool}}[0]."\t".$perc_dis."\t".$abundance."\t".${$$hdatafun{$title}{$pool}}[3]."\t".${$$hdatafun{$title}{$pool}}[4]."\t".${$$hdatafun{$title}{$pool}}[11]."-".${$$hdatafun{$title}{$pool}}[12]."\t".$title."\t".$pool."\t".$len."\t".$infoalign."\t".$hpvclosest."\t".$classification."\t".$seqf);
				}
				
			
				##On increment le nombre du prochain fichier temporaire créé
				$k++;
			}
			elsif(defined($$hdatafun{$title}{$pool})){
				##############################################################################
				##	Et si une seule sequence pour ce cluster et qu'il existe pour ce pool	##
				##############################################################################
				
				##	Need to be intialized for each sequence otherwise it's finish when it's raised the end of the file
				my $in=Bio::SeqIO->new( -file => $inputdirfasta."/".$p.".fasta",		##Initialisation du fichier contenant les contigs
				 -format => 'Fasta');
				
				my $abundance=nearest(.0001,((${$$hdatafun{$title}{$pool}}[3]/$reads{$pool})*100));
				my $perc_dis=nearest(.01,(100.00-(${$$hdatafun{$title}{$pool}}[1])));
				my $len = length(${$$hdatafun{$title}{$pool}}[9]);
				
				my $hpvclosest="";
				my $classification="";
				my @start_stop_len_hit=();
				
				##blastn HPV
				while (my $seq = $in->next_seq){
				
					my $id=$seq->display_id();
					
					if($id eq ${$$hdatafun{$title}{$pool}}[10]){
						#print $id."\n";
						
						my $out="$tmpfun/$id.blast";
						
						my $in2=HPVblastn($seq,$out);
				
						my $result=$in2->next_result;									
						my $var=$result->query_name();
						my $hit=$result->next_hit();
						
						#~ if($hit->num_hsps>1){
							#~ print $title."\n";
							#~ foreach my $hsp ($hit->next_hsp()){
								#~ print $hsp->percent_identity."\n";
								#~ print $hsp->score."\n";
							#~ }
						#~ }
						
						my $frame_query="";
						my $frame_hit="";
						
						if(defined $hit){
							my $hit_name = $hit->name();
							my $desc = $hit->description();
			
							#~ print $desc." la\n";		##	By modifying structure of naming for fasta file, allow to do the blast, but not to get gi, so have to take it from name or description
							#~ print $hit_name." ici\n";

							my $classification_tmp="";
							
							$classification_tmp=$hhpvtotax{$hit_name};
						
							if(!(defined($classification_tmp)) or ($classification_tmp eq "")){
								$classification_tmp="Unclassified";
							}
							#~ my @tab_name=split(/\|/,$hit_name);
							#~ my $gi = $tab_name[1];
							#~ chomp($gi);

							#~ my $taxid="";

							#~ $taxid=$gi_tax{$gi};
							#~ chomp($taxid);
							#~ my $classification_tmp="";
							#~ if(defined($taxid2spe{$taxid})){
								#~ $classification_tmp=$taxid2spe{$taxid};
							#~ }
							#~ elsif(defined($taxid2gen{$taxid})){
								#~ $classification_tmp=$taxid2gen{$taxid};
							#~ }
							#~ else{
								#~ $classification_tmp="Unclassified";
							#~ }
							
							my $percent_id_moy="NA";
							my $hsp=$hit->next_hsp();
							my @percent_id=();
							my $compte=0;
							my @start_stop_len_hsp=();
							
							#~ print "new hit\n";
							
							if(defined $hsp){
								do{
									$compte++;
									push @percent_id, $hsp->percent_identity;
									
									my $queryaln=$hsp->start."-".$hsp->end."(".$hsp->length('query').")";
									#~ print $hit_name."\n";
									#~ print $queryaln."\n";
									#~ print $queryaln."\n";
									#~ print $hsp->rank." rank\n";
									push @start_stop_len_hsp, $queryaln;

									$frame_query=$hsp->query->strand;
									$frame_hit=$hsp->hit->strand;
									
									#~ print "HSP : $compte\tFrame query : $frame_query\tFrame hit : $frame_hit\tFrame subject : $frame_subject\t$frame_subject_hit\n";
									
									if($frame_query!=$frame_hit){			##Reverse the sequence in order to have the good orientation before RaxML
										${$$hdatafun{$title}{$pool}}[9]=reverse_complement_IUPAC($seq->seq);
									}
									else{
										${$$hdatafun{$title}{$pool}}[9]=$seq->seq;
									}
									
								}while(my $hsp=$hit->next_hsp());
							}
							
							#~ if($compte>1){
								#~ print $title." dernier\n";
							#~ }
							
							@start_stop_len_hsp=uniq(@start_stop_len_hsp);
							
							my $sum=0;
							map { $sum += $_ } @percent_id;
							$percent_id_moy=nearest(.01,(($sum/$compte)));
							push(@start_stop_len_hit,join(";",@start_stop_len_hsp));

							$hpvclosest=$hit_name."($percent_id_moy%)";
							$classification=$classification_tmp;
							
						}
						else{
							$hpvclosest="NA";
							$classification="NA";
						}
						#~ push(@start_stop_len_hit,join(";",@start_stop_len_hsp));
					}
				}
				my $infoalign=join("/",@start_stop_len_hit);
				
				#print $hpvclosest."par la \n";
				if($boolean_infofile eq "true"){
					push(@aecrirefun,${$$hdatafun{$title}{$pool}}[0]."\t".$perc_dis."\t".$abundance."\t".${$$hdatafun{$title}{$pool}}[3]."\t".${$$hdatafun{$title}{$pool}}[4]."\t".${$$hdatafun{$title}{$pool}}[11]."-".${$$hdatafun{$title}{$pool}}[12]."\t".$title."\t".$pool."\t".${$$hdatafun{$title}{$pool}}[6]."\t".${$$hdatafun{$title}{$pool}}[7]."\t".$len."\t".$infoalign."\t".$hpvclosest."\t".$classification."\t".${$$hdatafun{$title}{$pool}}[9]);
				
				}else{
					#~ push(@aecrirefun,${$$hdatafun{$title}{$pool}}[0]."\t".$perc_dis."\t".$abundance."\t".${$$hdatafun{$title}{$pool}}[3]."\t".${$$hdatafun{$title}{$pool}}[4]."\t".$pool."\t".$title."\t".$len."\t".$infoalign."\t".$hpvclosest."\t".$classification."\t".${$$hdatafun{$title}{$pool}}[9]);
					push(@aecrirefun,${$$hdatafun{$title}{$pool}}[0]."\t".$perc_dis."\t".$abundance."\t".${$$hdatafun{$title}{$pool}}[3]."\t".${$$hdatafun{$title}{$pool}}[4]."\t".${$$hdatafun{$title}{$pool}}[11]."-".${$$hdatafun{$title}{$pool}}[12]."\t".$title."\t".$pool."\t".$len."\t".$infoalign."\t".$hpvclosest."\t".$classification."\t".${$$hdatafun{$title}{$pool}}[9]);
				}
				#~ print "$title\t$perc_dis\n";
			}
		}
	}
	$fac->cleanup;
	return(@aecrirefun);
}

my $bol_new_empty="F";
if(@aecrire != 0){		##Si des nouvelles sequences ont été trouvées
	##########
	##	NEW	##
	##########
	open(OUT,">".$outputdir."/table_putative_new_VIRUS.txt") or die "$! : $outputdir/table_putative_new_VIRUS.txt\n";
	open(FA,">".$outputdir."/putative_new_VIRUS.fa") or die "$! : $outputdir/putative_new_VIRUS.fa\n";
	if($boolean_infofile eq "true"){
		print OUT "VIRUSname\t%dissimilarity\tAbundance\tN°reads\tGInum\tAlignmentPosition_MegaBlast\tVIRUS_closest_MegaBlast\tPool\tTissu\tPrimer\tLength\tAlignmentPositionBlastN_start:stop(length)\tVIRUS_closest_Blast\tClassification\tSequence(s)\n";
	}
	else{
		print OUT "VIRUSname\t%dissimilarity\tAbundance\tN°reads\tGInum\tAlignmentPosition_MegaBlast\tVIRUS_closest_MegaBlast\tPool\tLength\tAlignmentPosition_start:stop(length)\tVIRUS_closest_Blast\tClassification\tSequence(s)\n";
	}
	foreach my $line (@aecrire){
		print OUT $line."\n";
		my @tab=split /\t/, $line;
		if($#tab>14){
			my $cpt=1;
			for (my $i=14; $i<=$#tab; $i++){
				print FA ">".$tab[0].".".$cpt."\n".$tab[$i]."\n";
				$cpt++;
			}
		}
		else{
			print FA ">".$tab[0]."\n".$tab[$#tab]."\n";
		}
	}
	close(FA);
	close(OUT);
}
else{
	$bol_new_empty="T";
}

my $bol_known_empty="F";
if(@aecrire3 != 0){		##Si des nouvelles sequences ont été trouvées
	##############
	##	KNOWN	##
	##############
	open(OUT,">".$outputdir."/table_putative_known_VIRUS.txt") or die "$! : $outputdir/table_putative_known_VIRUS.txt\n";
	open(FA,">".$outputdir."/putative_known_VIRUS.fa") or die "$! : $outputdir/putative_known_VIRUS.fa\n";
	if($boolean_infofile eq "true"){
		print OUT "VIRUSname\t%dissimilarity\tAbundance\tN°reads\tGInum\tAlignmentPosition_MegaBlast\tVIRUS_closest_MegaBlast\tPool\tTissu\tPrimer\tLength\tAlignmentPositionBlastN_start:stop(length)\tVIRUS_closest_Blast\tClassification\tSequence(s)\n";
	}
	else{
		print OUT "VIRUSname\t%dissimilarity\tAbundance\tN°reads\tGInum\tAlignmentPosition_MegaBlast\tVIRUS_closest_MegaBlast\tPool\tLength\tAlignmentPosition_start:stop(length)\tVIRUS_closest_Blast\tClassification\tSequence(s)\n";
	}
	foreach my $line (@aecrire3){
		print OUT $line."\n";
		my @tab=split /\t/, $line;
		if($#tab>14){
			my $cpt=1;
			for (my $i=14; $i<=$#tab; $i++){
				print FA ">".$tab[0].".".$cpt."\n".$tab[$i]."\n";
				$cpt++;
			}
		}
		else{
			print FA ">".$tab[0]."\n".$tab[$#tab]."\n";
		}
	}
	close(FA);
	close(OUT);
}
else{
	$bol_known_empty="T";
}


##############
##	RaxML	##
##############
##

print "\n Raxml identification\n";
my $newfasta=$outputdir."/putative_new_VIRUS.fa";
my $knownfasta=$outputdir."/putative_known_VIRUS.fa";

#~ `raxml.sh $whereiam $newfasta $knownfasta`;

mkdir $whereiam."/raxml" or die ($!);

mkdir $whereiam."/raxml/new" or die ($!);
mkdir $whereiam."/raxml/known" or die ($!);


if($bol_new_empty eq "F"){
	##	New
	chdir "$whereiam/raxml/new";

	`papara -t $dirname/raxml/L1_All_genome_NUC_GTRGAMMA_tree_newick.nwk -s $dirname/raxml/L1_All_genome_NUC_alignment.phylip -q $newfasta -j $threads`;

	`raxmlHPC-PTHREADS-AVX2 -f v -s papara_alignment.default $dirname/raxml/L1_All_genome_NUC_GTRGAMMA_tree_newick.nwk -m GTRGAMMA -n new -T $threads --epa-keep-placements=1`;
}

if($bol_known_empty eq "F"){
	##	Known
	chdir "$whereiam/raxml/known";

	`papara -t $dirname/raxml/L1_All_genome_NUC_GTRGAMMA_tree_newick.nwk -s $dirname/raxml/L1_All_genome_NUC_alignment.phylip -q $knownfasta -j $threads`;

	`raxmlHPC-PTHREADS-AVX2 -f v -s papara_alignment.default -t $dirname/raxml/L1_All_genome_NUC_GTRGAMMA_tree_newick.nwk -m GTRGAMMA -n known -T $threads --epa-keep-placements=1`;
}

chdir "$whereiam/raxml";


my $class_new="$whereiam/raxml/new/RAxML_classification.new";
my $class_known="$whereiam/raxml/known/RAxML_classification.known";

my $new_tree="$whereiam/raxml/new/RAxML_labelledTree.new";
my $known_tree="$whereiam/raxml/known/RAxML_labelledTree.known";

my $info_new="$whereiam/raxml/new/RAxML_info.new";
my $info_known="$whereiam/raxml/known/RAxML_info.known";

my $out_new="$whereiam/raxml/new/Raxml_new.csv";
my $out_known="$whereiam/raxml/known/Raxml_known.csv";

#~ my $knwon_json_tree="$whereiam/raxml/known/RAxML_portableTree.known.jplace";
#~ my $new_json_tree="$whereiam/raxml/new/RAxML_portableTree.new.jplace";


#~ parseRaxML($class_new,$new_tree,$info_new,$out_new);
#~ parseRaxML($class_known,$known_tree,$info_known,$out_known);
if($bol_new_empty eq "F"){
	parseRaxMLv2($new_tree,$info_new,$out_new);
}
if($bol_known_empty eq "F"){
	parseRaxMLv2($known_tree,$info_known,$out_known);
}

my %querytaxall=();
my %querycloseall=();

if($bol_new_empty eq "F"){
	readraxml($out_new);
}
if($bol_known_empty eq "F"){
	readraxml($out_known);
}

my $putative_new=$outputdir."/table_putative_new_VIRUS.txt";
my $putative_known=$outputdir."/table_putative_known_VIRUS.txt";

##	Re-intitalisation of variable
%hnewtable=();		#For RaxML
my %hnewtable_blastn=();
%hseq=();
my @catknown=();
my @catnew=();
%known=();
%new=();

my @catknown_blastn=();
my @catnew_blastn=();
my %known_blastn=();
my %new_blastn=();

if($bol_new_empty eq "F"){
	my $outnew=$outputdir."/table_putative_new_VIRUS_RaxML.txt";
	open(T,">".$outnew) or die "$! : $outnew\n";
	treatputative($putative_new,$out_new);
	close(T);
}

if($bol_known_empty eq "F"){
	my $outknown=$outputdir."/table_putative_known_VIRUS_RaxML.txt";
	open(T,">".$outknown) or die "$! : $outknown\n";
	treatputative($putative_known,$out_known);
	close(T);
}


##	This rewrite the final output, and include the RaxML classification
##	IF SEVERAL SEQUENCE FOR A CLUSTER, TAKE THE CLASSIFICATION OF THE LONGUEST ONE
sub treatputative{
	my $file=shift;
	my $fileraxml=shift;
	
	my $type="";
	if($file=~/new/){
		$type="new";
	}
	else{
		$type="known";
	}
	
	if($boolean_infofile eq "true"){
		#print T "HPVname\t%dissimilarity\tAbundance\tN°reads\tGInum\tPool\tHPV_closest_MegaBlast\tTissu\tPrimer\tLength\tAlignmentPosition_start:stop(length)\tHPV_closest_PaVE_Blast\tClassification\tClassificationRaxML\tSequence(s)\n";
		print T "VIRUSname\t%dissimilarity\tAbundance\tN°reads\tGInum\tAlignmentPosition_MegaBlast\tVIRUS_closest_MegaBlast\tPool\tTissu\tPrimer\tLength\tAlignmentPositionBlastN_start:stop(length)\tVIRUS_closest_Blast\tBlastN_Classification\tRaxML_closest_PV\tRaxML_Classification\tSequence(s)\n";
	}
	else{
		print T "VIRUSname\t%dissimilarity\tAbundance\tN°reads\tGInum\tAlignmentPosition_MegaBlast\tVIRUS_closest_MegaBlast\tPool\tLength\tAlignmentPositionBlastN_start:stop(length)\tVIRUS_closest_Blast\tBlastN_Classification\tRaxML_closest_PV\tRaxML_Classification\tSequence(s)\n";
	}

	
	my %querytax=();
	my %queryclose=();
	open(RAX,$fileraxml) or die "$! : $fileraxml\n";
	while(<RAX>){
		chomp($_);
		my $l=$_;
		my @tab=split /\t/,$l;
		if(defined($hacctohpv{$tab[1]})){
			$queryclose{$tab[0]}=$hacctohpv{$tab[1]};	#1005VIRUSput.1	MF588684	Betapapillomavirus
		}
		else{
			$queryclose{$tab[0]}="NA";
		}

		$querytax{$tab[0]}=$tab[2];
	}
	close(RAX);
	
	
	open(IN,$file) or die "$! : $file\n";
	my $header=<IN>;
	
	
	while(<IN>){
		chomp($_);
		my $l=$_;
		my @tab=split /\t/,$l;
		
		my $idselected="";			##	The id selected in case of several sequence in the cluster
		
		my $id="";
		my $perc_dis="";
		my $abun="";
		my $nbread="";
		my $gi="";
		my $alig_pos_mega="";
		my $hpv_mega="";
		my $pool="";
		my $tissu="";
		my $primer="";
		my $len="";
		my $align_pos_blastn="";
		my $hpv_blast="";
		my $classification_blastn="";
		my $seq="";
		
		my $raxml_class="";
		my $raxml_closest="";
		
		my $blastn_species="";
		my $blastn_closest="";
		
		if($boolean_infofile eq "true"){
			
			$id=$tab[0];
			$perc_dis=$tab[1];
			$abun=$tab[2];
			$nbread=$tab[3];
			$gi=$tab[4];
			$alig_pos_mega=$tab[5];
			$hpv_mega=$tab[6];
			$pool=$tab[7];
			$tissu=$tab[8];
			$primer=$tab[9];
			$len=$tab[10];
			$align_pos_blastn=$tab[11];
			$hpv_blast=$tab[12];
			$classification_blastn=$tab[13];
			$seq="";
			$raxml_class="";
			
			if($#tab>14){				##	Plus d'une sequence pour ce cluster
				my $cpt=1;
				my @tlen=split /\//,$len;
				my $imax=getindexmax(\@tlen);
				
				my $index=$imax+1;
				
				$idselected=$id.".".$index;
				
				my @taxtmp=();
				my @closetmp=();
				
				for(my $i=14; $i<=$#tab;$i++){
					$seq.=$tab[$i]."\t";
					$hseq{$id.$cpt}=$tab[$i];
					
					my $idtmp=$id.".".$cpt;
					
					#~ print $idtmp."\n";
					#~ print $querytax{$idtmp}."\n";
					
					push(@taxtmp,$querytax{$idtmp});
					push(@closetmp,$queryclose{$idtmp});	#Here is missing
					
					#~ if(!(defined($queryclose{$idtmp}))){
						#~ print $idtmp."\n";
						#~ print $l."\n";
						#~ print $#tab."\n";
						#~ exit;
					#~ }
					
					$cpt++;
				}
				$raxml_class=join('/',@taxtmp);
				$raxml_closest=join('/',@closetmp);
				
				if($hpv_blast=~/\//){	#more than one seq
					my @tabtmp=split(/\//,$hpv_blast);
					$blastn_closest=$tabtmp[0];
				}
				if($classification_blastn=~/\//){	#more than one seq
					my @tabtmp=split(/\//,$classification_blastn);
					$blastn_species=$tabtmp[0];
				}
				
			}
			else{						##	Une seule sequence pour ce cluser
				
				$seq=$tab[14];
				$hseq{$id}=$tab[14];
				$idselected=$id;
				$raxml_class=$querytax{$idselected};
				$raxml_closest=$queryclose{$idselected};
				$blastn_closest=$hpv_blast;
				$blastn_species=$classification_blastn;
			}
		}
		else{	##NO INFOFILE, LESS COLUMN		
			
			$id=$tab[0];
			$perc_dis=$tab[1];
			$abun=$tab[2];
			$nbread=$tab[3];
			$gi=$tab[4];
			$alig_pos_mega=$tab[5];
			$hpv_mega=$tab[6];
			$pool=$tab[7];
			$len=$tab[8];
			$align_pos_blastn=$tab[9];
			$hpv_blast=$tab[10];
			$classification_blastn=$tab[11];
			$seq="";
			$raxml_class="";
			
			if($#tab>12){				##	Plus d'une sequence pour ce cluster
				my $cpt=1;
				my @tlen=split /\//,$len;
				my $imax=getindexmax(\@tlen);
				
				my $index=$imax+1;
				
				$idselected=$id.".".$index;
				
				my @taxtmp=();
				my @closetmp=();
				
				for(my $i=12; $i<=$#tab;$i++){
					$seq.=$tab[$i]."\t";
					$hseq{$id.$cpt}=$tab[$i];
					
					my $idtmp=$id.".".$cpt;
					
					#~ print $idtmp."\n";
					#~ print $querytax{$idtmp}."\n";
					
					push(@taxtmp,$querytax{$idtmp});
					push(@closetmp,$queryclose{$idtmp});
					
					$cpt++;
				}
				$raxml_class=join('/',@taxtmp);
				$raxml_closest=join('/',@closetmp);
				
				if($hpv_blast=~/\//){	#more than one seq
					my @tabtmp=split(/\//,$hpv_blast);
					$blastn_closest=$tabtmp[0];
				}
				if($classification_blastn=~/\//){	#more than one seq
					my @tabtmp=split(/\//,$classification_blastn);
					$blastn_species=$tabtmp[0];
				}
				
			}
			else{						##	Une seule sequence pour ce cluser
				
				$seq=$tab[12];
				$hseq{$id}=$tab[12];
				$idselected=$id;
				$raxml_class=$querytax{$idselected};
				$raxml_closest=$queryclose{$idselected};
				$blastn_closest=$hpv_blast;
				$blastn_species=$classification_blastn;
			}
			
		}
		
		##to modifiy
		
		#~ $@{$hnewtable{$tissu{$pool}}{$cattmp}{$hpvcat{$hpvname}}{$hpvname}}[$poolsid{$pool}]+=$htarget{$pool}{'known'}{$hpvname};
		
		my $genus="";
		#~ print $querytax{$idselected}." rax\n";
		if($querytax{$idselected}=~/^(\S+)papillomavirus (\d+)/){
			$genus=$1;
		}
		elsif($querytax{$idselected}=~/(\S+)papillomavirus/){
			if($1=~/alpha/i or $1=~/beta/i or $1=~/gamma/i or $1=~/mu/i or $1=~/nu/i){
				$genus="Unreferenced";
			}
			else{
				$genus=$1;
			}
		}
		else{
			$genus=$querytax{$idselected};
		}
		#~ print $genus."\n";
		
		my $genus_blastn="";
		#~ print $blastn_species." blast\n";
		if($blastn_species=~/(\S+)papillomavirus (\d+)/){
			$genus_blastn=$1;
		}
		elsif($blastn_species=~/(\S+)papillomavirus/){
			if($1=~/alpha/i or $1=~/beta/i or $1=~/gamma/i or $1=~/mu/i or $1=~/nu/i){
				$genus_blastn="Unreferenced";
			}
			else{
				$genus_blastn=$1;
			}
		}
		else{
			$genus_blastn="Unclassified";
		}
		
		#{$tissu}{$hpvgenus}{$hpvspecies}{$hpvname}=@(nbreads_pool1,nbreads_pool2,...,nbreads_pool8)
		#~ print $tissu."\n";
		#~ print $genus."\n";
		#~ if(!(defined($genus))){
			#~ print $idselected."\n";
		#~ }
		#~ print $querytax{$idselected}."\n";
		#~ print $queryclose{$idselected}."\n";
		#~ print $poolsid{$pool}."\n";
		#~ print $nbread."\n";
		
		if(!(defined($queryclose{$idselected}))){
			print $idselected."\n";
			print $l."\n";
			exit;
		}
		
		
		if(!(defined($hnewtable{$tissu}{$genus}{$querytax{$idselected}}{$queryclose{$idselected}}))){
			@{$hnewtable{$tissu}{$genus}{$querytax{$idselected}}{$queryclose{$idselected}}}=();
		}
		
		my $blastn_tmp=$blastn_closest;
		if(($blastn_tmp=~/(.*)\(\d+\.\d+%\)/) or ($blastn_tmp=~/(.*)\(\d+%\)/)){
			$blastn_closest=$1;
		}
		
		##	Pour Blastn
		if(!(defined($hnewtable_blastn{$tissu}{$genus_blastn}{$blastn_species}{$blastn_closest}))){
			@{$hnewtable_blastn{$tissu}{$genus_blastn}{$blastn_species}{$blastn_closest}}=();
		}
		
		$@{$hnewtable{$tissu}{$genus}{$querytax{$idselected}}{$queryclose{$idselected}}}[$poolsid{$pool}]+=$nbread;
		
		$@{$hnewtable_blastn{$tissu}{$genus_blastn}{$blastn_species}{$blastn_closest}}[$poolsid{$pool}]+=$nbread;
		
		my $simplegenus=$genus;
		#~ if($querytax{$idselected}=~/^(\S+)papillomavirus/){
			#~ $simplegenus=$1;
			#~ if($simplegenus=~/human/i){
				#~ $simplegenus=$querytax{$idselected};
			#~ }
		#~ }
		#~ else{
			#~ $simplegenus=$querytax{$idselected};
		#~ }
		
		
		#
		#	BlastN
		#
		if($type eq "new"){
			$new_blastn{$pool}{$genus_blastn}+=$nbread;
		}
		else{
			$known_blastn{$pool}{$genus_blastn}+=$nbread;
		}
		
		foreach my $cat (keys %{$known_blastn{$pool}}){
			if (! grep {$_ eq $cat} @catknown_blastn){
				push(@catknown_blastn,$cat);
			}
			
		}
		foreach my $cat (keys %{$new_blastn{$pool}}){
			if (! grep {$_ eq $cat} @catnew_blastn){
				#~ print "Category new :".$cat."\n";
				push(@catnew_blastn,$cat);
			}
		}
		
		#
		#	RaxML
		#
		if(defined($querytaxall{$type}{$idselected})){
			if($type eq "known"){
				$known{$pool}{$simplegenus}+=$nbread;
				#~ push(@catknown,$querytax{$idselected});
			}
			else{
				$new{$pool}{$simplegenus}+=$nbread;
				#~ push(@catnew,$querytax{$idselected});
			}
		}
		
		foreach my $cat (keys %{$known{$pool}}){
			if (! grep {$_ eq $cat} @catknown){
				push(@catknown,$cat);
			}
			
		}
		foreach my $cat (keys %{$new{$pool}}){
			if (! grep {$_ eq $cat} @catnew){
				#~ print "Category new :".$cat."\n";
				push(@catnew,$cat);
			}
		}
		
		if($boolean_infofile eq "true"){
			print T "$id\t$perc_dis\t$abun\t$nbread\t$gi\t$alig_pos_mega\t$hpv_mega\t$pool\t$tissu\t$primer\t$len\t$align_pos_blastn\t$hpv_blast\t$classification_blastn\t$raxml_closest\t$raxml_class\t$seq\n";
		}
		else{
			print T "$id\t$perc_dis\t$abun\t$nbread\t$gi\t$alig_pos_mega\t$hpv_mega\t$primer\t$len\t$align_pos_blastn\t$hpv_blast\t$classification_blastn\t$raxml_closest\t$raxml_class\t$seq\n";
		}
		
	}
	close(IN);
}


#################################################			###		RAXML		###
##	WRITE TABLE HPV DIVERSITY BY TISSUE TYPE	##
##################################################
##				!!	NEW VERSION	!!				##
##################################################
##	2 independant table : 1 for SKIN and 1 for ORAL

$krona=$outputdir."/KRONA_RaxML";
mkdir $krona;

%hoverall=();	# {Family}{Genus}{Species}{virusname}{tissu}=nb reads

$fg=0;
$fs=0;
$fn=0;
foreach my $t (sort keys %hnewtable){		##TISSU
	
	my $kronaname="";
	
	my @primertmp_no_uniq=();
	my @poolsidtmp=();
	foreach my $p (sort @pools){			# Pour chaque nom ordered
		if($tissu{$p} eq $t){				# Si le tissu du pollname correspond a celui procedé
			push(@primertmp_no_uniq,$primer{$p});		#le primer de ce pool										$p = NGSDEEP_DL2A1_S11_L001 --> $primer{$p} = 11
			push(@poolsidtmp,$poolsid{$p});		###	Both table are in the same order	$p = NGSDEEP_DL2A1_S11_L001 --> $poolsid{$p} = 	
		}
	}
	my @primertmp=();
	@primertmp=uniq @primertmp_no_uniq;
	
	if($boolean_infofile eq "true"){
		open(OUT,">".$outputdir."/diversityByTissu_".$t."_RaxML.csv") or die "$!";
		print OUT "TISSUE\tFamily\tGenus\tSpecies\tRelated\t".join("\t",sort @primertmp)."\tPool\n";		#vérifier ordre primertmp
		print OUT $t."\tPapillomaviridae\t";
		$kronaname="krona_".$t."_RaxML";
		open(KRO,">".$krona."/krona_".$t."_RaxML.txt") or die "$!";
	}
	else{
		open(OUT,">".$outputdir."/OverallDiversity_RaxML.csv") or die "$!";
		print OUT "Family\tFamily\tGenus\tSpecies\tRelated\t".join("\t",sort @primertmp)."\tPool\n";		#vérifier ordre primertmp
		print OUT "Virus\t";
		
		$kronaname="krona_ByTissue_RaxML";
		open(KRO,">".$krona."/krona_ByTissue_RaxML.txt") or die "$!";
	}
	
	
	
	foreach my $g (sort keys %{$hnewtable{$t}}){			##GENUS
		if($fg==0){
			print OUT $g."\t";
			$fg=1;
		}
		else{
			print OUT "\t\t".$g."\t";
		}
		foreach my $s (sort keys %{$hnewtable{$t}{$g}}){		##SPECIES
			if($fs==0){
				print OUT $s."\t";
				$fs=1;
			}
			else{
				print OUT "\t\t\t".$s."\t";
			}
			foreach my $n (sort keys %{$hnewtable{$t}{$g}{$s}}){		##NAME
				my @tmp_pool=();
				my %tmp_primer=();
				if($fn==0){
					print OUT $n."\t";
					$fn=1;
				}
				else{
					print OUT "\t\t\t\t".$n."\t";
				}
				my $tabtmp=$@{$hnewtable{$t}{$g}{$s}{$n}};
				#~ my $cluster=@$tabtmp[0]." (".nearest(.01,((@$tabtmp[0]/$clusttissue{$t})*100)).")";

				#~ for my $i (@poolsidtmp){
					#~ if((defined(@$tabtmp[$i])) && @$tabtmp>0){
						#~ my %hpoolsidrev= reverse %poolsid;
						#~ my $p=$hpoolsidrev{$i};
						#~ $tmp_primer{$p}+=@$tabtmp[$i];
						#~ push(@tmp_pool,$p);
					#~ }
				#~ }
				#~ my $tmp_towrite="";
				#~ my $total_krona=0;
				#~ foreach my $i (sort @poolsidtmp){
					#~ my %hpoolsidrev= reverse %poolsid;
					#~ my $p=$hpoolsidrev{$i};
					#~ if(exists($tmp_primer{$p})){
						#~ $tmp_towrite.=$tmp_primer{$p}."\t";
						#~ $total_krona+=$tmp_primer{$p};
					#~ }
					#~ else{
						#~ $tmp_towrite.="-\t";
					#~ }
				#~ }
				
				for my $i (@poolsidtmp){	#	3	4	5	6	7	8	9	10	11	
					if((defined(@$tabtmp[$i])) && @$tabtmp>0){		#	3	4	5
						my %hpoolsidrev= reverse %poolsid;			#	poolid vers poolname	
						my $p=$hpoolsidrev{$i};						#	$p=NGSDEEP_DL2A1_S11_L001	--> $i=3
						$tmp_primer{$primer{$p}}+=@$tabtmp[$i];		#	$primer{$p} = nb reads of pool number X having this primer		
						push(@tmp_pool,$p);			#	get pool name - $p = NGSDEEP_DL2A1_S11_L001
					}
				}
				my $tmp_towrite="";
				my $total_krona=0;
				foreach my $prim (sort @primertmp){#	all primer concerned
					if(exists($tmp_primer{$prim})){		# tmp_primer de NGSDEEP_DL2A1_S11_L001 --> nb reads of pool number (3)	
						$tmp_towrite.=$tmp_primer{$prim}."\t";		##on ecrit
						$total_krona+=$tmp_primer{$prim};
					}
					else{
						$tmp_towrite.="-\t";
					}
				}
				
				chomp($tmp_towrite);
				print OUT $tmp_towrite.join(",",@tmp_pool)."\n";
				
				print KRO $total_krona."\t".$t."\t".$g."\t".$s."\n";
				
				$hoverall{'Papillomaviridae'}{$g}{$s}{$n}{$t}+=$total_krona;
				
				#~ my %hoverall=();	# {Family}{Genus}{Species}{virusname}{tissu}=nb reads


			}
			$fn=0;
		}
		$fs=0;
	}
	$fg=0;
	close(OUT);
	
	`ktImportText $krona/$kronaname.txt -o $krona/$kronaname.html`;
}


##########################################
##	WRITE TABLE HPV DIVERSITY OVERALL	##
##########################################
$ff=0;
$fg=0;
$fs=0;
$fn=0;
@sorttissu=();
foreach my $t (sort keys %hnewtable){		#For each tissu (sorted)
	push(@sorttissu,$t);
}

open(OUTALL,">".$outputdir."/DiversityByTissu_RaxML.csv") or die "$!";
print OUTALL "TISSUE\tFamily\tGenus\tSpecies\tRelated\t".join("\t",@sorttissu)."\tPool\n";		#vérifier ordre primertmp
print OUTALL "ALL\t";

open(KRO,">".$krona."/krona_OverAll_RaxML.txt") or die "$!";

foreach my $f (sort keys %hoverall){		#For each family
	if($ff==0){
		print OUTALL $f."\t";
		$ff=1;
	}
	else{
		print OUTALL "\t".$f."\t";
	}
	foreach my $g (sort keys %{$hoverall{$f}}){		#For each genus (sorted)
		#~ print OUTALL $g."\t";
		if($fg==0){
			print OUTALL $g."\t";
			$fg=1;
		}
		else{
			print OUTALL "\t\t".$g."\t";
		}
		foreach my $s (sort keys %{$hoverall{$f}{$g}}){	#For each species
			if($fs==0){
				print OUTALL $s."\t";
				$fs=1;
			}
			else{
				print OUTALL "\t\t\t".$s."\t";
			}
			foreach my $n (sort keys %{$hoverall{$f}{$g}{$s}}){		#For each polyoname
				if($fn==0){
					print OUTALL $n."\t";
					$fn=1;
				}
				else{
					print OUTALL "\t\t\t\t".$n."\t";
				}
				
				my @pool_list=();
				
				my $total_krona=0;
				foreach my $t (@sorttissu){
					my $nbread="-";
					if(defined($hoverall{$f}{$g}{$s}{$n}{$t})){
						$nbread=$hoverall{$f}{$g}{$s}{$n}{$t};
						$total_krona+=$hoverall{$f}{$g}{$s}{$n}{$t};
						
							## To display the pool concerned	$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[$poolsid{$pool}]=$htarget{$pool}{'new'}{$virusname};
						my $tabtmp;
						#~ $@{$hnewtable{$tissu}{$genus}{$querytax{$idselected}}{$queryclose{$idselected}}}[$poolsid{$pool}]+=$nbread;
						#{$tissu}{$hpvgenus}{$hpvspecies}{$hpvname}=@(nbreads_pool1,nbreads_pool2,...,nbreads_pool8)
						
						#~ my $tabtmp=$@{$hnewtable{$t}{$g}{$s}{$n}};
						
						if(defined($hnewtable{$t}{$g}{$s}{$n})){
							$tabtmp=$@{$hnewtable{$t}{$g}{$s}{$n}};
						}
						#~ else{
							#~ print $hoverall{$f}{$g}{$s}{$n}{$t}."\n";
							#~ print $t."\n";
							#~ print $f."\n";
							#~ print $g."\n";
							#~ print $s."\n";
							#~ print $n."\n";
							#~ print "Error\n";
							#~ exit;
						#~ }
						for(my $i=1; $i<=$#{$tabtmp}; $i++){		# 0 is kept for the number of cluster
							if((defined(@$tabtmp[$i])) && @$tabtmp>0){			#	3	4	5
								my %hpoolsidrev= reverse %poolsid;				#	poolid vers poolname	
								my $p=$hpoolsidrev{$i};							
								push(@pool_list,$p);
							}
						}
					
					}
					print OUTALL $nbread."\t";	
					
				}
				print OUTALL join(',',@pool_list);
				print OUTALL "\n";
				
				print KRO $total_krona."\tALL\t".$g."\t".$s."\n";
			}
			$fn=0;
		}
		$fs=0;
	}
	$fg=0;
}
$ff=0;
close(KRO);
close(OUTALL);

`ktImportText $krona/krona_OverAll_RaxML.txt -o $krona/krona_OverAll_RaxML.html`;




#################################################			###		BLASTN		###
##	WRITE TABLE HPV DIVERSITY BY TISSUE TYPE	##
##################################################
##				!!	NEW VERSION	!!				##
##################################################
##	2 independant table : 1 for SKIN and 1 for ORAL

$krona=$outputdir."/KRONA_BlastN";
mkdir $krona;

my %hoverall_blastn=();	# {Family}{Genus}{Species}{virusname}{tissu}=nb reads

$fg=0;
$fs=0;
$fn=0;
foreach my $t (sort keys %hnewtable_blastn){		##TISSU
	
	my $kronaname="";
	
	my @primertmp_no_uniq=();
	my @poolsidtmp=();
	foreach my $p (sort @pools){			# Pour chaque nom ordered
		if($tissu{$p} eq $t){				# Si le tissu du pollname correspond a celui procedé
			push(@primertmp_no_uniq,$primer{$p});		#le primer de ce pool										$p = NGSDEEP_DL2A1_S11_L001 --> $primer{$p} = 11
			push(@poolsidtmp,$poolsid{$p});		###	Both table are in the same order	$p = NGSDEEP_DL2A1_S11_L001 --> $poolsid{$p} = 	
		}
	}
	my @primertmp=();
	@primertmp=uniq @primertmp_no_uniq;
	
	if($boolean_infofile eq "true"){
		open(OUT,">".$outputdir."/diversityByTissu_".$t."_BlastN.csv") or die "$!";
		print OUT "TISSUE\tFamily\tGenus\tSpecies\tRelated\t".join("\t",sort @primertmp)."\tPool\n";		#vérifier ordre primertmp
		print OUT $t."\tPapillomaviridae\t";
		$kronaname="krona_".$t."_BlastN";
		open(KRO,">".$krona."/krona_".$t."_BlastN.txt") or die "$!";
	}
	else{
		open(OUT,">".$outputdir."/OverallDiversity_BlastN.csv") or die "$!";
		print OUT "Family\tFamily\tGenus\tSpecies\tRelated\t".join("\t",sort @primertmp)."\tPool\n";		#vérifier ordre primertmp
		print OUT "Virus\t";
		
		$kronaname="krona_ByTissue_BlastN";
		open(KRO,">".$krona."/krona_ByTissue_BlastN.txt") or die "$!";
	}
	
	
	
	foreach my $g (sort keys %{$hnewtable_blastn{$t}}){			##GENUS
		if($fg==0){
			print OUT $g."\t";
			$fg=1;
		}
		else{
			print OUT "\t\t".$g."\t";
		}
		foreach my $s (sort keys %{$hnewtable_blastn{$t}{$g}}){		##SPECIES
			if($fs==0){
				print OUT $s."\t";
				$fs=1;
			}
			else{
				print OUT "\t\t\t".$s."\t";
			}
			foreach my $n (sort keys %{$hnewtable_blastn{$t}{$g}{$s}}){		##NAME
				my @tmp_pool=();
				my %tmp_primer=();
				if($fn==0){
					print OUT $n."\t";
					$fn=1;
				}
				else{
					print OUT "\t\t\t\t".$n."\t";
				}
				my $tabtmp=$@{$hnewtable_blastn{$t}{$g}{$s}{$n}};
				#~ my $cluster=@$tabtmp[0]." (".nearest(.01,((@$tabtmp[0]/$clusttissue{$t})*100)).")";

				#~ for my $i (@poolsidtmp){
					#~ if((defined(@$tabtmp[$i])) && @$tabtmp>0){
						#~ my %hpoolsidrev= reverse %poolsid;
						#~ my $p=$hpoolsidrev{$i};
						#~ $tmp_primer{$p}+=@$tabtmp[$i];
						#~ push(@tmp_pool,$p);
					#~ }
				#~ }
				#~ my $tmp_towrite="";
				#~ my $total_krona=0;
				#~ foreach my $i (sort @poolsidtmp){
					#~ my %hpoolsidrev= reverse %poolsid;
					#~ my $p=$hpoolsidrev{$i};
					#~ if(exists($tmp_primer{$p})){
						#~ $tmp_towrite.=$tmp_primer{$p}."\t";
						#~ $total_krona+=$tmp_primer{$p};
					#~ }
					#~ else{
						#~ $tmp_towrite.="-\t";
					#~ }
				#~ }
				
				for my $i (@poolsidtmp){	#	3	4	5	6	7	8	9	10	11	
					if((defined(@$tabtmp[$i])) && @$tabtmp>0){		#	3	4	5
						my %hpoolsidrev= reverse %poolsid;			#	poolid vers poolname	
						my $p=$hpoolsidrev{$i};						#	$p=NGSDEEP_DL2A1_S11_L001	--> $i=3
						$tmp_primer{$primer{$p}}+=@$tabtmp[$i];		#	$primer{$p} = nb reads of pool number X having this primer		
						push(@tmp_pool,$p);			#	get pool name - $p = NGSDEEP_DL2A1_S11_L001
					}
				}
				my $tmp_towrite="";
				my $total_krona=0;
				foreach my $prim (sort @primertmp){#	all primer concerned
					if(exists($tmp_primer{$prim})){		# tmp_primer de NGSDEEP_DL2A1_S11_L001 --> nb reads of pool number (3)	
						$tmp_towrite.=$tmp_primer{$prim}."\t";		##on ecrit
						$total_krona+=$tmp_primer{$prim};
					}
					else{
						$tmp_towrite.="-\t";
					}
				}
				
				chomp($tmp_towrite);
				print OUT $tmp_towrite.join(",",@tmp_pool)."\n";
				
				print KRO $total_krona."\t".$t."\t".$g."\t".$s."\n";
				
				$hoverall_blastn{'Papillomaviridae'}{$g}{$s}{$n}{$t}+=$total_krona;
				
				#~ my %hoverall=();	# {Family}{Genus}{Species}{virusname}{tissu}=nb reads


			}
			$fn=0;
		}
		$fs=0;
	}
	$fg=0;
	close(OUT);
	
	`ktImportText $krona/$kronaname.txt -o $krona/$kronaname.html`;
}



##########################################			###		BLASTN		###
##	WRITE TABLE HPV DIVERSITY OVERALL	##
##########################################
$ff=0;
$fg=0;
$fs=0;
$fn=0;
@sorttissu=();
foreach my $t (sort keys %hnewtable_blastn){		#For each tissu (sorted)
	push(@sorttissu,$t);
}

open(OUTALL,">".$outputdir."/DiversityByTissu_BlastN.csv") or die "$!";
print OUTALL "TISSUE\tFamily\tGenus\tSpecies\tRelated\t".join("\t",@sorttissu)."\tPool\n";		#vérifier ordre primertmp
print OUTALL "ALL\t";

open(KRO,">".$krona."/krona_OverAll_BlastN.txt") or die "$!";

foreach my $f (sort keys %hoverall_blastn){		#For each family
	if($ff==0){
		print OUTALL $f."\t";
		$ff=1;
	}
	else{
		print OUTALL "\t".$f."\t";
	}
	foreach my $g (sort keys %{$hoverall_blastn{$f}}){		#For each genus (sorted)
		#~ print OUTALL $g."\t";
		if($fg==0){
			print OUTALL $g."\t";
			$fg=1;
		}
		else{
			print OUTALL "\t\t".$g."\t";
		}
		foreach my $s (sort keys %{$hoverall_blastn{$f}{$g}}){	#For each species
			if($fs==0){
				print OUTALL $s."\t";
				$fs=1;
			}
			else{
				print OUTALL "\t\t\t".$s."\t";
			}
			foreach my $n (sort keys %{$hoverall_blastn{$f}{$g}{$s}}){		#For each polyoname
				if($fn==0){
					print OUTALL $n."\t";
					$fn=1;
				}
				else{
					print OUTALL "\t\t\t\t".$n."\t";
				}
				
				my @pool_list=();
				
				my $total_krona=0;
				foreach my $t (@sorttissu){
					my $nbread="-";
					if(defined($hoverall_blastn{$f}{$g}{$s}{$n}{$t})){
						$nbread=$hoverall_blastn{$f}{$g}{$s}{$n}{$t};
						$total_krona+=$hoverall_blastn{$f}{$g}{$s}{$n}{$t};
						
							## To display the pool concerned	$@{$hnewtable{$tissu{$pool}}{$family}{$genus}{$species}{$virusname}}[$poolsid{$pool}]=$htarget{$pool}{'new'}{$virusname};
						my $tabtmp;
						#~ $@{$hnewtable{$tissu}{$genus}{$querytax{$idselected}}{$queryclose{$idselected}}}[$poolsid{$pool}]+=$nbread;
						#{$tissu}{$hpvgenus}{$hpvspecies}{$hpvname}=@(nbreads_pool1,nbreads_pool2,...,nbreads_pool8)
						
						#~ my $tabtmp=$@{$hnewtable{$t}{$g}{$s}{$n}};
						
						if(defined($hnewtable_blastn{$t}{$g}{$s}{$n})){
							$tabtmp=$@{$hnewtable_blastn{$t}{$g}{$s}{$n}};
						}
						#~ else{
							#~ print $hoverall{$f}{$g}{$s}{$n}{$t}."\n";
							#~ print $t."\n";
							#~ print $f."\n";
							#~ print $g."\n";
							#~ print $s."\n";
							#~ print $n."\n";
							#~ print "Error\n";
							#~ exit;
						#~ }
						for(my $i=1; $i<=$#{$tabtmp}; $i++){		# 0 is kept for the number of cluster
							if((defined(@$tabtmp[$i])) && @$tabtmp>0){			#	3	4	5
								my %hpoolsidrev= reverse %poolsid;				#	poolid vers poolname	
								my $p=$hpoolsidrev{$i};							
								push(@pool_list,$p);
							}
						}
					
					}
					print OUTALL $nbread."\t";	
					
				}
				print OUTALL join(',',@pool_list);
				print OUTALL "\n";
				
				print KRO $total_krona."\tALL\t".$g."\t".$s."\n";
			}
			$fn=0;
		}
		$fs=0;
	}
	$fg=0;
}
$ff=0;
close(KRO);
close(OUTALL);

`ktImportText $krona/krona_OverAll_BlastN.txt -o $krona/krona_OverAll_BlastN.html`;








##################################
##	WRITE TABLE SUMMARY RAXML	##
##################################
my $outsummary="$outputdir/table_summary_RaxML.csv";

##########################################
##	Table 9 - KNOWN virus genus level	##
##########################################
open(OUT,">".$outsummary) or die "$! : $outsummary\n";
print OUT "\nTable 9 - KNOWN virus genus level - RaxML\n";
print OUT "Pool";
foreach my $cat (sort @catknown){
	print OUT "\t".$cat;
}
print OUT "\n";
foreach my $pool (sort keys %known){
	print OUT $pool;
	foreach my $cat (sort @catknown){
		if(defined($known{$pool}{$cat})){
			print OUT "\t".$known{$pool}{$cat};
		}
		else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
######################################
##	Table 10 - NEW virus genus level	##
######################################
print OUT "\nTable 10 - NEW virus genus level - RaxML\n";
print OUT "Pool";
##	WRITE TABLE NEW HPV CATEGORY
foreach my $cat (sort @catnew){
	print OUT "\t".$cat;
}
print OUT "\n";
foreach my $pool (sort keys %new){
	print OUT $pool;
	foreach my $cat (sort @catnew){
		if(defined($new{$pool}{$cat})){
			print OUT "\t".$new{$pool}{$cat};
		}
		else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
print OUT "\n";


##########################################
##	Table 11 - KNOWN virus genus level	##
##########################################
print OUT "\nTable 11 - KNOWN virus genus level - BlastN\n";
print OUT "Pool";
foreach my $cat (sort @catknown_blastn){
	print OUT "\t".$cat;
}
print OUT "\n";
foreach my $pool (sort keys %known_blastn){
	print OUT $pool;
	foreach my $cat (sort @catknown_blastn){
		if(defined($known_blastn{$pool}{$cat})){
			print OUT "\t".$known_blastn{$pool}{$cat};
		}
		else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
######################################
##	Table 12 - NEW virus genus level	##
######################################
print OUT "\nTable 12 - NEW virus genus level - BlastN\n";
print OUT "Pool";
##	WRITE TABLE NEW HPV CATEGORY
foreach my $cat (sort @catnew_blastn){
	print OUT "\t".$cat;
}
print OUT "\n";
foreach my $pool (sort keys %new_blastn){
	print OUT $pool;
	foreach my $cat (sort @catnew_blastn){
		if(defined($new_blastn{$pool}{$cat})){
			print OUT "\t".$new_blastn{$pool}{$cat};
		}
		else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close(OUT);


##	Get index of max value in a table
sub getindexmax{
	my $myarray=shift;
	my $index = 0;
	my $maxval = @$myarray[$index];
	
	#~ print $maxval."\n";
	#~ print $maxval;
	
	for(my $i=0; $i<=$#{$myarray}; $i++){
		if($maxval < @$myarray[$i]){
			$maxval = @$myarray[$i];
			$index = $i;
		}
	}
	return($index);
}


## This function read the file produce et the rpeviosu step (see parse RaxML result function)
sub readraxml{
	my $file=shift;
	my $type="";
	if($file=~/new/){
		$type="new";
	}
	else{
		$type="known";
	}
	open(IN,$file) or die "$! : $file\n";
	while(<IN>){
		chomp($_);
		my $l=$_;
		my @tab=split /\t/,$l;
		$querycloseall{$type}{$tab[0]}=$tab[1];
		$querytaxall{$type}{$tab[0]}=$tab[2];
	}
	close(IN);
}


##################################################
##	Second Version or parsing RaxML EPA result	##
##################################################
sub parseRaxMLv2{
	my ($raxtree,$info,$out)=@_;
	
	my %hidentical=();
	open(INFO, $info) or die "$! : $info\n";	#IMPORTANT WARNING: Sequences 97HPVput and 44HPVput are exactly identical
	while(<INFO>){
		chomp($_);
		my $line=$_;
		if($line=~/IMPORTANT WARNING: Sequences (.*) and (.*) are exactly identical/){
			$hidentical{"QUERY___".$1}="QUERY___".$2;
			#IMPORTANT WARNING: Sequences 757VIRUSput and 1069VIRUSput.2 are exactly identical
			#~ $hidentical{757VIRUSput}=1069VIRUSput.2;
			
			
		}
	}
	close(INFO);

	my $treeio = Bio::TreeIO->new(-format => 'newick',
					  -file => $raxtree);

	my $tree = $treeio->next_tree;


	my @taxa = $tree->get_leaf_nodes;

	my $ancestor;
	undef($ancestor);

	my %hquerytohpv=();

	foreach my $node (@taxa){
		if($node->id=~/^QUERY/){
			#~ print "\t".$node->id."\n";
			my $ancestor=$node->ancestor;	#ancetre de QUERY___X (1 seulement)
			while(defined($ancestor)){
				my @childs=();
				foreach my $child ($ancestor->each_Descendent){			#	Pour chaque enfant de l'ancetre de QUERY___X
					my $node_name=$child->id;
					if(($child->is_Leaf) && ($node_name!~/^QUERY/)){	# si un de ces frere est connu et n'est pas une query
						#~ print $child->id." ici\n";
						push(@childs,$child);
					}
				}
				if($#childs<0){					## Pas de childs known				
					$ancestor=$ancestor->ancestor;
					
					if(!(defined($ancestor))){
						#~ print "No ancestor for ".$node->id."\n";
						$hquerytohpv{$node->id}="NA";		# Alors pas dans l'arbre...
					}
					
				}
				else{
					$hquerytohpv{$node->id}=$childs[0]->id;
					undef($ancestor);
				}
			}
		}
	}
	
	open(OUT,">".$out) or die "$! : $out\n";
	foreach my $q (sort keys %hquerytohpv){
		my $acc="";
		my $tax="";
		if((defined($hquerytohpv{$q}))){
			$acc=$hquerytohpv{$q};
			if((defined($hacctotax{$hquerytohpv{$q}}))){
				$tax=$hacctotax{$hquerytohpv{$q}};				#Modifier ce n'est plus acc
			}
			else{
				$tax="Unclassified";
			}
		}
		else{
			$acc="NA";
		}
		my $id="";
		if($q=~/^QUERY___(.+)/){
			$id=$1;
		}
		print OUT $id."\t".$acc."\t".$tax."\n";
	}
	foreach my $iden (sort keys %hidentical){
		my $q=$hidentical{$iden};
		my $acc="";
		my $tax="";
		if((defined($hquerytohpv{$iden}))){
			$acc=$hquerytohpv{$iden};
			if((defined($hacctotax{$hquerytohpv{$iden}}))){
				$tax=$hacctotax{$hquerytohpv{$iden}};
			}
			else{
				$tax="Unclassified";
			}
		}
		else{
			$acc="NA";
		}
		my $id="";
		if($q=~/^QUERY___(.+)/){
			$id=$1;
		}
		print OUT $id."\t".$acc."\t".$tax."\n";
	}
	close(OUT);
}


##
#SequIO object to blast
#Blast output file name
#return Bio::SearchIO object result
sub HPVblastn{
	my $seqfun=shift;
	my $outfun=shift;
	my $resultfun = $fac->blastn( -query => $seqfun,										#Faire le blast avec la sequence $seq -db_data => '/data/robitaillea/blastHPVdb/HPV.fasta',	
			 -outfile => $outfun, 
			 #~ -db_data => '/data/robitaillea/blastHPVdb/HPV.fasta',
			 -method_args => [	-num_alignments => '1',			#	recupération de 50 alignements
							-show_gis => 'yes',					#	le nom du resultat doit contenir le numéro GI
							-word_size => '11',					#	taille du mot pour effectuer la recherche
							-evalue => '1e-1',					#	e-value attendu est de 1e-5							-->1e-4
							-num_threads => $threads,				#	Nombre de CPU a utiliser				#	Pourcentage d'identité minimum attendu
							-penalty => '-3',					#	Penalitée du score du blast pour un mismatch		-->	-1
							-reward => '2',						#	Gain du score du blast pour un match				-->	1
							-gapopen => '5',					#														-->	0
							-gapextend => '2',					#														-->	2
							] );
							
	my $in2fun=Bio::SearchIO->new(	-format => 'blast',								#Ouverture en lecture du fichier de resultat blast
								-file => $outfun);
	return $in2fun;
}


sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

sub reverse_complement_IUPAC {
        my $dna = shift;

        # reverse the DNA sequence
        my $revcomp = reverse($dna);

        # complement the reversed DNA sequence
        #~ $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        $revcomp =~ tr/ATGCatgc/TACGtacg/;
        return $revcomp;
}

#~ my $res=longuestFirst(\@lentemp,\@hpvclosest,\classification,\@start_stop_len_hit,\@tseqtmpfun);

sub longuestFirst {
	my $tlen=shift;
	my $thpv=shift;
	my $tclass=shift;
	my $tpos=shift;
	my $tseq=shift;
	
	my @tindex=();
	
	my @sorted_len = sort { $a <=> $b } @$tlen;	# Plus petit en premier
	
	my @tlentemp = @$tlen;
	
	while(@sorted_len){
		my $longuest = pop @sorted_len;
		#~ print "longuest ".$longuest."\n";
		##	Traiter cas ou plusieurs sequence même taille --> pour le moment first index prend tjrs la meme sequence
		push @tindex, first_index { $_ eq $longuest } @tlentemp;
		
		if($tindex[$#tindex] == -1){
			$tindex[$#tindex] = 0;
		}
		
		#~ print "last index ".$#tindex." and associated value ".$tindex[$#tindex]."\n";
		$tlentemp[$tindex[$#tindex]]="undef";
	}
	
	#~ foreach my $i (@tindex){
		#~ print $i."\n";
	#~ }
	
	
	#~ while(@$tlen){
		#~ my $indextop=getindexmax($tlen);
		#~ print $indextop." ici \n";
		#~ push @tindex, $indextop;
		#splice @$tlen, $indextop, 1;	## Pas bon car cahnge index des restant// use undef?
		#undef $$tlen[$indextop];	## Pas bon car cahnge index des restant// use undef?
	#~ }
	
	my @tlenfun=();
	my @thpvfun=();
	my @tclassfun=();
	my @posfun=();
	my @tseqfun=();
	
	foreach my $i (@tindex){
		#~ print $i."\n";
		#~ print $$tlen[$i]." la \n";
		push(@tlenfun,$$tlen[$i]);			## 
		push(@thpvfun,$$thpv[$i]);
		push(@tclassfun,$$tclass[$i]);
		push(@posfun,$$tpos[$i]);
		push(@tseqfun,$$tseq[$i]);
	}
	
	my $lenfun=join("/",@tlenfun); 
	my $hpvfun=join("/",@thpvfun); 
	my $classfun=join("/",@tclassfun); 
	my $posfun=join("/",@posfun); 
	my $seqfun=join("\t",@tseqfun);
	
	my @res=();
	push(@res,$lenfun);
	push(@res,$hpvfun);
	push(@res,$classfun);
	push(@res,$posfun);
	push(@res,$seqfun);
	 
	return(\@res);
}
