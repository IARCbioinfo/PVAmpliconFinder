#/usr/bin/perl -s

##############################################################
##	Name2ACC.pl												##
##	18/06/2019												##
##	Amplicon sequencing Illumina MiSeq						##
##	Alexis Robitaille : robitaillea@students.iarc.fr		##
##	IARC, LYON												##
##	Last modification = 18/06/2019							##
##	Version 1.0												##
##############################################################

#	This script is used to format the MSA and RT to contains accession number instead of HPV number using Pave.table

##################
##	Libraries	##
##################
use strict;
use warnings;
use Cwd;
use Bio::TreeIO;

my $whereiam=getcwd;

my $pave_table=$whereiam."/../raxml/pave_table.txt";

my $old_tree=$whereiam."/../raxml/PaVE-4059-40949_consensus.nwk";
my $old_align=$whereiam."/../raxml/PaVE_align.fasta";

my $new_tree=$whereiam."/../raxml/L1_All_genome_NUC_GTRGAMMA_tree_newick.nwk";
my $new_align=$whereiam."/../raxml/L1_All_genome_NUC_alignment.fasta";

my %hconversion=();

open(F1, $pave_table) or die "$!: $pave_table\n";

while(<F1>)
{
	my $line=$_;
	chomp($line);
	
	my @tab=split("\t",$line);	#	BPV19	KU519394	Unclassified
	$hconversion{$tab[0]}=$tab[1];
	#~ print $tab[0]."\t".$tab[1]."\n";
}
		
close(F1);


my @tree=();
my @subtree=();
my $treetmp="";
my $subtreetmp="";
my $subsubtreetmp="";

open(OLDTREE, $old_tree) or die "$!: $old_tree\n";
		
while(<OLDTREE>)
{
	my $line=$_;
	$treetmp=$line;
	#~ chomp($line);
	#~ print $line."\n";
	#~ @tree=split(',',$line);
	#~ foreach my $k (keys %hconversion){
		
	#~ }
}
		
close(OLDTREE);


my @sorted_by_length = sort { length($b) <=> length($a) } keys %hconversion;

foreach my $name (@sorted_by_length){
	$treetmp =~ s/$name/$hconversion{$name}/i;
}
#~ print $treetmp;


open(NEWTREE, ">".$new_tree) or die "$!: $new_tree\n";
print NEWTREE $treetmp;
close(NEWTREE);


my $align="";

#~ my $value="";

#~ print $hconversion{'AmPV3'};
#~ exit;

open(OLDALIGN, $old_align) or die "$!: $old_align\n";
while(<OLDALIGN>)
{
	my $line=$_;
	chomp($line);
	#~ print $line."\n";
	#~ print $line;
	
	my $towrite="";
	
	if($line =~ /^>(\S+)/){
		my $piece=$1;
		#~ print $piece;
		$towrite=">".$hconversion{$piece};
		#~ print $towrite;
		#~ foreach my $name (keys %hconversion){
			#~ if($line=~/$name/){
				#~ $line =~ s/$name/$hconversion{$name}/i;
				#~ print $line;
				#~ last;
			#~ }
		#~ }
	}
	else{
		$towrite=$line;
	}
	#~ if($line !~ /^\s+/){
		#~ print $line."\n";
		#~ my @tab=split('-',$line);
		#~ chomp($tab[0]);
		#~ print $tab[0];
		#~ my $value=$hconversion{$tab[0]};
		#~ print $value."\n";
		#~ foreach my $name (keys %hconversion){
			#~ if($line=~/$name/){
				#~ print $line;
				#~ $line =~ s/$name/$hconversion{$name}/i;
				#~ print $line;
				#~ push
				#~ last;
			#~ }
		#~ }
	#~ }
	
	
	$align.=$towrite."\n";
}
		
close(OLDALIGN);


open(NEWALIGN, ">".$new_align) or die "$!: $new_align\n";
print NEWALIGN $align;
close(NEWALIGN);


