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
#~ use Statistics::R;
use List::MoreUtils qw(uniq);
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::Index::Fasta;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SearchIO;
use Bio::Seq;
use Text::CSV;
use Cwd;
use JSON::Parse 'parse_json';
use JSON::Parse 'json_file_to_perl';
use JSON;
use Data::Dumper;
use Bio::TreeIO;


my $table="/home/robitaillea/ICB/NGS_luisa/raxml/pave_table.txt";

my %hacctohpv=();
my %hacctotax=();
open(T, $table) or die "$! : $table\n";	#HPV1	V01116	Mupapillomavirus 1
while(<T>){
	chomp($_);
	my @tab=split(/\t/,$_);
	if((defined($tab[1])) && $tab[1]!~/^\s*$/){
		$hacctohpv{$tab[1]}=$tab[2];		#	BPV7	DQ217793	Dyoxipapillomavirus 1
	}
	if((defined($tab[2])) && $tab[2]!~/^\s*$/){
		$hacctotax{$tab[1]}=$tab[2];		#	BPV7	DQ217793	Dyoxipapillomavirus 1
	}
}
close(T);

#~ print $hacctotax{"EU493092"}."\n";
#~ exit;

my $class_known="/home/robitaillea/ICB/NGS_luisa/raxml/known/RAxML_classification.known";

my $known_tree="/home/robitaillea/ICB/NGS_luisa/raxml/known/RAxML_labelledTree.known";

my $info_known="/home/robitaillea/ICB/NGS_luisa/raxml/known/RAxML_info.known";

my $out_known="/home/robitaillea/ICB/NGS_luisa/raxml/known/Raxml_known.csv";

my $known_labeltree="/home/robitaillea/ICB/NGS_luisa/raxml/known/RAxML_originalLabelledTree.known";

#~ parseRaxML($class_known,$known_tree,$info_known,$out_known);
my $known_tree_json="/home/robitaillea/ICB/NGS_luisa/raxml/known/RAxML_portableTree.known.jplace";

parseRaxMLv2($known_tree,$info_known,$out_known);

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
						$hquerytohpv{$node->id}="NA";
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
				$tax=$hacctotax{$hquerytohpv{$q}};
			}
			else{
				$tax="Unclassified";
			}
		}
		else{
			$acc="NA";
		}
		print OUT $q."\t".$acc."\t".$tax."\n";
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
		print OUT $q."\t".$acc."\t".$tax."\n";
	}
	close(OUT);
}


exit;

my %hhpvtoacc;
my %hhpvtotax;
my $treeio = Bio::TreeIO->new(-format => 'newick',
			      -file => $known_tree);

my $tree = $treeio->next_tree;


my @taxa = $tree->get_leaf_nodes;

my $node = $tree->find_node(-id => 'QUERY___70VIRUSput.7');


my $ancestor=$node->ancestor;	#ancetre de QUERY___X (1 seulement)
while($ancestor){
	my @childs=();
	foreach my $child ($ancestor->each_Descendent){			#	Pour chaque enfant de l'ancetre de QUERY___70VIRUSput.7
		my $node_name=$child->id;
		if(($child->is_Leaf) && ($node_name!~/^QUERY/)){	# si un de ces frere est connu
			print $child->id." ici\n";
			push(@childs,$child);
		}
	}
	if($#childs<0){
		$ancestor=$ancestor->ancestor;
	}
	else{
		undef($ancestor);
	}
}



undef($ancestor);

foreach my $node (@taxa){
	if($node->id=~/^QUERY/){
		print "\t".$node->id."\n";
		my $ancestor=$node->ancestor;	#ancetre de QUERY___X (1 seulement)
		while(defined($ancestor)){
			my @childs=();
			foreach my $child ($ancestor->each_Descendent){			#	Pour chaque enfant de l'ancetre de QUERY___70VIRUSput.7
				my $node_name=$child->id;
				if(($child->is_Leaf) && ($node_name!~/^QUERY/)){	# si un de ces frere est connu et n'est pas une query
					print $child->id." ici\n";
					push(@childs,$child);
				}
			}
			if($#childs<0){
				$ancestor=$ancestor->ancestor;
				
				if(!(defined($ancestor))){
					print "No ancestor for ".$node->id."\n";
				}
				
			}
			else{
				undef($ancestor);
			}
		}
	}
}
#~ for my $child ($ancestor->each_Descendent){
	#~ print $child->id." lala\n";
#~ }

#~ my $ancestor2=$ancestor->ancestor;
exit;

#~ if($ancestor2->is_Leaf){
	#~ print $ancestor2->id."\n";
#~ }
#~ else{
	

#~ for my $child2 ($ancestor2->each_Descendent){
	#~ if($child2->is_Leaf){
		#~ print $child2->id." lala\n";
	#~ }
	#~ else{
		
#~ }
	
	
#~ my $lca = $tree->get_lca(-nodes => @nodes);

#~ print $lca->id."\n";


exit;

      
#~ while( my $tree = $treeio->next_tree ) {
	#~ for my $node ( $tree->get_nodes ) {
	#~ printf "id: %s bootstrap: %s\n", 
		 #~ $node->id || '', 
		 #~ $node->bootstrap || '', 
		 #~ "\n";
	#~ }
#~ }		      
			      
#~ exit;		      
			      
#~ my $t="";
#~ open(JSON, $known_tree_json) or die "$! : $known_tree_json\n";	#IMPORTANT WARNING: Sequences 97HPVput and 44HPVput are exactly identical
#~ while(<JSON>){
	#~ chomp($_);
	#~ $t.=$_;
#~ }
#~ close(JSON);




#~ my $p = json_file_to_perl ($known_tree_json);
#~ my $text = decode_json($t);
#~ print Dumper($text);

#~ foreach my $key (sort keys $p){
	#~ print $key."\t".$p->{$key}."\n";
#~ }

#~ my $fields=$p->{'fields'};	#ARRAY
#~ my $metadata=$p->{'metadata'};	#HASH
#~ my $placements=$p->{'placements'};	#ARRAY
#~ my $tree =$p->{'tree'};	#STRING


#~ foreach my $elem (@$fields){
	#~ print $elem."\n";
#~ }

#~ foreach my $elem (@$placements){
	#~ print $elem."(elem)\n";		#HASH
	#~ foreach my $k (sort keys $elem){
		#~ print $k."(k)\t";		#n ou p
		#~ my $tab=$elem->{$k};	# ARRAY
		#~ print $tab."(tab)\t";		
		#~ foreach my $tmp (@$tab){
			#~ print $tmp."(tmp)\n";	#query seq ou Array
			#~ my @tab2=$tmp;
		#~ }
	#~ }
#~ }

#~ foreach my $elem (@$placements){
	#~ print Dumper $elem."\n";
	#~ my $toto=Dumper($elem);
	#~ print $toto."\n";
#~ }

#~ exit;
#~ foreach my $elem (@$placements){
	#~ my $toto=Dumper($elem);
	
	#~ my $data = $elem;
	
	#~ my @fields=();
	#~ foreach my $key (sort keys $elem){
		#~ print "toto";
		#~ push(@fields,$key);
	#~ }

	#~ my $csv = Text::CSV->new({auto_diag=>2,binary=>1, eol=>"\n", always_quote=>1 });

	#~ $csv->print(select, \@fields);
	#~ for my $datapoint ( @{ $data->{Datapoints} } ) {
		#~ $csv->print(select, [ map {$datapoint->{$_}} @fields ]);
	#~ }
#~ }


#~ foreach my $elem (@$placements){
	#~ foreach my $key (sort keys $elem){
		#~ print $elem->{'p'}[0][0]."\t";	#placement
		#~ print $elem->{'n'}[0]."\n";		#names
	#~ }
#~ }


sub parseRaxML{
	my ($class,$raxtree,$info,$out)=@_;
	
	my %hidentical=();
	open(INFO, $info) or die "$! : $info\n";	#IMPORTANT WARNING: Sequences 97HPVput and 44HPVput are exactly identical
	while(<INFO>){
		chomp($_);
		my $line=$_;
		if($line=~/IMPORTANT WARNING: Sequences (.*) and (.*) are exactly identical/){
			$hidentical{$1}=$2;
			#IMPORTANT WARNING: Sequences 757VIRUSput and 1069VIRUSput.2 are exactly identical
			#~ $hidentical{757VIRUSput}=1069VIRUSput.2;
			
			
		}
	}
	close(INFO);
	
	
	my %hhpvpos=();
	open(CLA, $class) or die "$! : $class\n";
	while(<CLA>){
		chomp($_);
		my @tab=split(/\s/,$_);
		$hhpvpos{$tab[0]}=$tab[1];
		#	EX : 757VIRUSput	I315
		# hhpvpos{757VIRUSput}  =	I315
	}
	close(CLA);

	my $tree="";
	open(TREE, $raxtree) or die "$! : $raxtree\n";
	while(<TREE>){
		chomp($_);
		my $line=$_;
		if($line!~/^\s*$/){
			$tree=$line;
		}
	}
	close(TREE);
	#~ print $tree;

	my %hidtoclose=();

	my @tab=split(/\[I\d+\]/, $tree);		##THIS produce a table with the different insertion position

	my %hquerytohpv=();

	foreach my $elem (@tab){
		if($elem=~/QUERY___/){	##Si il y a une query a cet endroit	
			foreach my $hpv (keys %hhpvtoacc){
				my $acc=$hhpvtoacc{$hpv};
				if($elem=~/$acc/){			#On identify the closest hpv acc number pour ces QUERYS
					my @q=split(/QUERY___/,$elem);
					foreach my $piece (@q){					
						if($piece=~/(\d+VIRUSput\.*\d*):/){
							my $toto=$1;
							$hquerytohpv{$toto}=$hpv;							
							#~ if($piece=~/757VIRUSput/){
								#~ print $toto."\t".$hpv." ICICICICI\n";
							#~ }
						}
					}
				}
			}
			#~ push(@tabquery,$elem);
		}
		else{
			print $elem."\n";
		}
	}
	
	open(OUT,">".$out) or die "$! : $out\n";
	foreach my $q (sort keys %hquerytohpv){
		print OUT $q."\t".$hquerytohpv{$q}."\t".$hhpvtotax{$hquerytohpv{$q}}."\n";
	}
	foreach my $iden (sort keys %hidentical){
		my $q=$hidentical{$iden};
		#~ print $iden." is like ".$hidentical{$iden}."\n";
		#~ print $iden."\t".$hquerytohpv{$iden}."\t".$hhpvtotax{$hquerytohpv{$iden}}."\n";
		print OUT $q."\t".$hquerytohpv{$iden}."\t".$hhpvtotax{$hquerytohpv{$iden}}."\n";
	}
	close(OUT);
}



##############################################
##	Old Version or parsing RaxML EPA result	##
##############################################
## Produce a table that include the id of the sequence, the closest PV based on the tree, and the classification associated with this PV
sub parseRaxML{
	my ($class,$raxtree,$info,$out)=@_;
	
	my %hidentical=();
	open(INFO, $info) or die "$! : $info\n";	#IMPORTANT WARNING: Sequences 97HPVput and 44HPVput are exactly identical
	while(<INFO>){
		chomp($_);
		my $line=$_;
		if($line=~/IMPORTANT WARNING: Sequences (.*) and (.*) are exactly identical/){
			$hidentical{$1}=$2;
			#IMPORTANT WARNING: Sequences 757VIRUSput and 1069VIRUSput.2 are exactly identical
			#~ $hidentical{757VIRUSput}=1069VIRUSput.2;
			
			
		}
	}
	close(INFO);
	
	
	my %hhpvpos=();
	open(CLA, $class) or die "$! : $class\n";
	while(<CLA>){
		chomp($_);
		my @tab=split(/\s/,$_);
		$hhpvpos{$tab[0]}=$tab[1];
		#	EX : 757VIRUSput	I315
		# hhpvpos{757VIRUSput}  =	I315
	}
	close(CLA);

	my $tree="";
	open(TREE, $raxtree) or die "$! : $raxtree\n";
	while(<TREE>){
		chomp($_);
		my $line=$_;
		if($line!~/^\s*$/){
			$tree=$line;
		}
	}
	close(TREE);
	#~ print $tree;

	my %hidtoclose=();

	my @tab=split(/\[I\d+\]/, $tree);		##THIS produce a table with the different insertion position

	my %hquerytohpv=();

	foreach my $elem (@tab){
		if($elem=~/QUERY___/){	##Si il y a une query a cet endroit	
			foreach my $hpv (keys %hhpvtoacc){
				my $acc=$hhpvtoacc{$hpv};
				if($elem=~/$acc/){			#On identify the closest hpv acc number pour ces QUERYS
					my @q=split(/QUERY___/,$elem);
					foreach my $piece (@q){					
						if($piece=~/(\d+VIRUSput\.*\d*):/){
							my $toto=$1;
							$hquerytohpv{$toto}=$hpv;							
							if($toto=~/757VIRUSput/){
								print $toto."\t".$hpv." ICICICICI\n";
							}
						}
					}
				}
			}
			#~ push(@tabquery,$elem);
		}
	}
	
	open(OUT,">".$out) or die "$! : $out\n";
	foreach my $q (sort keys %hquerytohpv){
		print OUT $q."\t".$hquerytohpv{$q}."\t".$hhpvtotax{$hquerytohpv{$q}}."\n";
	}
	foreach my $iden (sort keys %hidentical){
		my $q=$hidentical{$iden};
		print $iden." is like ".$hidentical{$iden}."\n";
		print $iden."\t".$hquerytohpv{$iden}."\t".$hhpvtotax{$hquerytohpv{$iden}}."\n";
		print OUT $q."\t".$hquerytohpv{$iden}."\t".$hhpvtotax{$hquerytohpv{$iden}}."\n";
	}
	close(OUT);
}
