#/usr/bin/perl -s

##############################################################
##	format_PaVE_db.pl										##
##	06/06/2019												##
##	Amplicon sequencing Illumina MiSeq						##
##	Alexis Robitaille : robitaillea@students.iarc.fr		##
##	IARC, LYON												##
##	Last modification = 06/06/2019							##
##	Version 1.0												##
##############################################################

#	This script is used to format the fasta file from PaVE db containing all L1 gene from PVs
#	This format is needed to correctly create the tree from RaxML
#	Indeed, the fasta header like ">gi|9627074.L1|lcl|AaPV1_L1.1| Alces alces Papillomavirus 1 (AaPV1), L1 gene" is not convenient for phylogeny analysis
#	It should be transfomed in "AaPV1", shorter and explicite

##################
##	Libraries	##
##################
use strict;
use warnings;
use Cwd;


my $whereiam=getcwd;

my $pave_fasta=$whereiam."/../databases/PaVE/PaVE.fasta";
my $formated_pave_fasta=$whereiam."/../databases/PaVE/PaVE_formated.fasta";


open(OUT, ">".$formated_pave_fasta) or die "$!: $formated_pave_fasta\n";

open(F1, $pave_fasta) or die "$!: $pave_fasta\n";
		
while(<F1>)
{
	my $line=$_;
	chomp($line);
	if($line=~/^>/){
		my @tab1=split('\|',$line);
		my $val=$tab1[3];
		my @tab2=split('_',$val);
		my $hpv=$tab2[0];
		print OUT ">".$hpv."\n";
	}
	else{
		print OUT $line."\n";
	}
}
		
close(F1);
close(OUT);
