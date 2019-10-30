#/usr/bin/perl -s

##############################################################
##	create_table_HPV_presence.pl							##
##	22/01/2019												##
##	Amplicon sequencing Illumina MiSeq	- 	Step 3			##
##	Alexis Robitaille : robitaillea@students.iarc.fr		##
##	IARC, LYON												##
##	Last modification = 30/10/2019							##
##	Version 1.0												##
##############################################################

#	The goal is to create a table with presence/absence of HPV type stratified by tissu type and known/new PV
#	HPV num	|	Tissue1	|	tissue2
#ex:HPV1	|	yes	|	no

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $in_folder="/home/robitaillea/ICB/NGS_luisa/output_correction_Luisa_2_final/HPV_presence_known_new";	

my $out_folder="/home/robitaillea/ICB/NGS_luisa/output_correction_Luisa_2_final/HPV_presence_known_new";	


chdir "$in_folder";

my $files=`ls table_putative_*_VIRUS_RaxML.txt`;

#~ print $files."\n";

my %hreadsmega=();		#megablast
my %hreadsrax=();		#raxml

my %hinfomega=();
my %hinforax=();


my %hreadsmega_genus=();		#megablast
my %hreadsrax_genus=();		#raxml

my %hseqcountmega_genus=();
my %hseqcountrax_genus=();

my %htypemega_genus=();		#megablast
my %htyperax_genus=();	

	
foreach my $file (`ls table_putative_*_VIRUS_RaxML.txt`){
	chomp($file);
	print $file."\n";
	my $group="";
	if($file=~/table_putative_(\S+)_VIRUS_RaxML.txt/){
		$group=$1;
		
		open(IN, $in_folder."/".$file) or die "$! : $in_folder/$file\n";	#
		my $header=<IN>;
		while(<IN>){
			chomp($_);
			my $l=$_;
			my @tab=split /\t/,$l;
			
			my $id=$tab[0];
			my $perc_dis=$tab[1];
			my $abun=$tab[2];
			my $nbread=$tab[3];
			my $gi=$tab[4];
			my $alig_pos_mega=$tab[5];
			my $hpv_mega=$tab[6];
			my $pool=$tab[7];
			my $tissu=$tab[8];
			my $primer=$tab[9];
			#~ my $primer=$tab[8];
			#~ my $tissu=$tab[9];
			my $len=$tab[10];	#OK
			my $align_pos_blastn=$tab[11];
			my $hpv_blast=$tab[12];
			my $classification_blastn=$tab[13];
			my $raxml_closest=$tab[14];	#OK
			my $raxml_class=$tab[15];
			my $seq=$tab[16];
			
			my @thpv_blast=split(/\//,$hpv_blast);
			my @tclassification_blastn=split(/\//,$classification_blastn);
			
			my @traxml_closest=split(/\//,$raxml_closest);
			my @traxml_class=split(/\//,$raxml_class);
			
			my @tlen=split(/\//,$len);
			
			my $max=getindexmax(\@tlen);
			
			#~ if($id=~/202VIRUSput/){
					#~ print $tlen[$max]."\n";
					#~ exit;
			#~ }
			
			$hpv_blast=$thpv_blast[$max];
			$classification_blastn=$tclassification_blastn[$max];
			$raxml_closest=$traxml_closest[$max];
			$raxml_class=$traxml_class[$max];
			

			my $t=$tissu;
			
			my $hpv="";
			if($hpv_blast=~/(.+)\(.+\)$/i){
				$hpv=$1;
			}
			else{
				$hpv=$hpv_blast;
			}
			#~ print $hpv."\n";
			#~ exit;
			
			my $genusrax=getgenus($raxml_class);
			my $genusblast=getgenus($classification_blastn);
			
			
			##	BLASTN
						#NEW	HPV1	AK	Primer
			$hreadsmega{$group}{$hpv}{$t}{$primer}+=$nbread;
			
			$hreadsmega_genus{$group}{$genusblast}{$t}{$primer}+=$nbread;
			
			if(!(defined($hseqcountmega_genus{$group}{$genusblast}{$t}{$primer}))){
				$hseqcountmega_genus{$group}{$genusblast}{$t}{$primer}=1;
			}
			else{
				$hseqcountmega_genus{$group}{$genusblast}{$t}{$primer}++;
			}
			
			#~ print $group."\t".$hpv."\t".$t."\t".$primer."\n";
			
			$hinfomega{$hpv}=$classification_blastn;
			
			push(@{$htypemega_genus{$group}{$genusblast}},$hpv);
			
			
			
			##	RAXML-EPA
			
			$hreadsrax{$group}{$raxml_closest}{$t}{$primer}+=$nbread;
			
			$hreadsrax_genus{$group}{$genusrax}{$t}{$primer}+=$nbread;
			
			if(!(defined($hseqcountrax_genus{$group}{$genusrax}{$t}{$primer}))){
				$hseqcountrax_genus{$group}{$genusrax}{$t}{$primer}=1;
			}
			else{
				$hseqcountrax_genus{$group}{$genusrax}{$t}{$primer}++;
			}
			
			$hinforax{$raxml_closest}=$raxml_class;
			
			#~ print $group."\t".$hpv."\t".$t."\t".$nbread."\t".$raxml_closest."\n";
			
			push(@{$htyperax_genus{$group}{$genusrax}},$raxml_closest);
			
		}
		close(IN);
	}
}

sub getgenus{
	my $hpvname=shift;
	my $genus="";
	if($hpvname=~/^alpha/i){
		$genus="Alphapapillomavirus";
	}
	elsif($hpvname=~/^beta/i){
		$genus="Betapapillomavirus";
	}
	elsif($hpvname=~/^gamma/i){
		$genus="Gammapapillomavirus";
	}
	elsif($hpvname=~/^mu/i){
		$genus="Mupapillomavirus";
	}
	elsif($hpvname=~/^nu/i){
		$genus="Nupapillomavirus";
	}
	elsif($hpvname=~/^unclass/i){
		$genus="Unclassified";
	}
	else{
		$genus="Non-human papillomavirus";
	}
	return($genus);
}


################################
################################	HSAK primer sep	$hreadsmega_genus{$group}{$genusblast}{$t}{$primer}+=$nbread;
################################
for my $group (sort keys %hreadsmega_genus){
	open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_tissu_BlastN_newVersion.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_tissu_BlastN_newVersion.txt\n";
	
	print OUT "HPV genus\tUnique virus species\tHS\t\tAK\t\n";
	print OUT "\t\tFAP\tFAPM1\tFAP\tFAPM1\n";
	
	for my $genus (sort keys %{$hreadsmega_genus{$group}}){

		my $fapreadshs="-";
		my $fapm1readshs="-";
		my $fapreadsak="-";
		my $fapm1readsak="-";
		
		my $fapseqshs="-";
		my $fapm1seqshs="-";
		my $fapseqsak="-";
		my $fapm1seqsak="-";
		
		##Reads
		if((defined($hreadsmega_genus{$group}{$genus}{'HS'}{'FAP'})) && ($hreadsmega_genus{$group}{$genus}{'HS'}{'FAP'}!~/^\s*$/)){
			$fapreadshs=$hreadsmega_genus{$group}{$genus}{'HS'}{'FAP'};
		}
		if((defined($hreadsmega_genus{$group}{$genus}{'HS'}{'FAPM1'})) && ($hreadsmega_genus{$group}{$genus}{'HS'}{'FAPM1'}!~/^\s*$/)){
			$fapm1readshs=$hreadsmega_genus{$group}{$genus}{'HS'}{'FAPM1'};
		}
		if((defined($hreadsmega_genus{$group}{$genus}{'AK'}{'FAP'})) && ($hreadsmega_genus{$group}{$genus}{'AK'}{'FAP'}!~/^\s*$/)){
			$fapreadsak=$hreadsmega_genus{$group}{$genus}{'AK'}{'FAP'};
		}
		if((defined($hreadsmega_genus{$group}{$genus}{'AK'}{'FAPM1'})) && ($hreadsmega_genus{$group}{$genus}{'AK'}{'FAPM1'}!~/^\s*$/)){
			$fapm1readsak=$hreadsmega_genus{$group}{$genus}{'AK'}{'FAPM1'};
		}
		
		##Seqs
		if((defined($hseqcountmega_genus{$group}{$genus}{'HS'}{'FAP'})) && ($hseqcountmega_genus{$group}{$genus}{'HS'}{'FAP'}!~/^\s*$/)){
			$fapseqshs=$hseqcountmega_genus{$group}{$genus}{'HS'}{'FAP'};
		}
		if((defined($hseqcountmega_genus{$group}{$genus}{'HS'}{'FAPM1'})) && ($hseqcountmega_genus{$group}{$genus}{'HS'}{'FAPM1'}!~/^\s*$/)){
			$fapm1seqshs=$hseqcountmega_genus{$group}{$genus}{'HS'}{'FAPM1'};
		}
		if((defined($hseqcountmega_genus{$group}{$genus}{'AK'}{'FAP'})) && ($hseqcountmega_genus{$group}{$genus}{'AK'}{'FAP'}!~/^\s*$/)){
			$fapseqsak=$hseqcountmega_genus{$group}{$genus}{'AK'}{'FAP'};
		}
		if((defined($hseqcountmega_genus{$group}{$genus}{'AK'}{'FAPM1'})) && ($hseqcountmega_genus{$group}{$genus}{'AK'}{'FAPM1'}!~/^\s*$/)){
			$fapm1seqsak=$hseqcountmega_genus{$group}{$genus}{'AK'}{'FAPM1'};
		}
		
		my @tabcount=uniq(@{$htypemega_genus{$group}{$genus}});
		my $countuniq=scalar @tabcount;
		
		print OUT $genus."\t".$countuniq."\t".$fapseqshs." (".$fapreadshs.")\t".$fapm1seqshs." (".$fapm1readshs.")\t".$fapseqsak." (".$fapreadsak.")\t".$fapm1seqsak." (".$fapm1readsak.")\n";
	}
	close(OUT);
}

for my $group (sort keys %hreadsrax_genus){
	open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_tissu_RaxML_newVersion.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_tissu_RaxML_newVersion.txt\n";
	print OUT "HPV genus\tUnique virus species\tHS\t\tAK\t\n";
	print OUT "\t\tFAP\tFAPM1\tFAP\tFAPM1\n";
	
	for my $genus (sort keys %{$hreadsrax_genus{$group}}){

		my $fapreadshs="-";
		my $fapm1readshs="-";
		my $fapreadsak="-";
		my $fapm1readsak="-";
		
		my $fapseqshs="-";
		my $fapm1seqshs="-";
		my $fapseqsak="-";
		my $fapm1seqsak="-";
		
		##Reads
		if((defined($hreadsrax_genus{$group}{$genus}{'HS'}{'FAP'})) && ($hreadsrax_genus{$group}{$genus}{'HS'}{'FAP'}!~/^\s*$/)){
			$fapreadshs=$hreadsrax_genus{$group}{$genus}{'HS'}{'FAP'};
		}
		if((defined($hreadsrax_genus{$group}{$genus}{'HS'}{'FAPM1'})) && ($hreadsrax_genus{$group}{$genus}{'HS'}{'FAPM1'}!~/^\s*$/)){
			$fapm1readshs=$hreadsrax_genus{$group}{$genus}{'HS'}{'FAPM1'};
		}
		if((defined($hreadsrax_genus{$group}{$genus}{'AK'}{'FAP'})) && ($hreadsrax_genus{$group}{$genus}{'AK'}{'FAP'}!~/^\s*$/)){
			$fapreadsak=$hreadsrax_genus{$group}{$genus}{'AK'}{'FAP'};
		}
		if((defined($hreadsrax_genus{$group}{$genus}{'AK'}{'FAPM1'})) && ($hreadsrax_genus{$group}{$genus}{'AK'}{'FAPM1'}!~/^\s*$/)){
			$fapm1readsak=$hreadsrax_genus{$group}{$genus}{'AK'}{'FAPM1'};
		}
		
		##Seqs
		if((defined($hseqcountrax_genus{$group}{$genus}{'HS'}{'FAP'})) && ($hseqcountrax_genus{$group}{$genus}{'HS'}{'FAP'}!~/^\s*$/)){
			$fapseqshs=$hseqcountrax_genus{$group}{$genus}{'HS'}{'FAP'};
		}
		if((defined($hseqcountrax_genus{$group}{$genus}{'HS'}{'FAPM1'})) && ($hseqcountrax_genus{$group}{$genus}{'HS'}{'FAPM1'}!~/^\s*$/)){
			$fapm1seqshs=$hseqcountrax_genus{$group}{$genus}{'HS'}{'FAPM1'};
		}
		if((defined($hseqcountrax_genus{$group}{$genus}{'AK'}{'FAP'})) && ($hseqcountrax_genus{$group}{$genus}{'AK'}{'FAP'}!~/^\s*$/)){
			$fapseqsak=$hseqcountrax_genus{$group}{$genus}{'AK'}{'FAP'};
		}
		if((defined($hseqcountrax_genus{$group}{$genus}{'AK'}{'FAPM1'})) && ($hseqcountrax_genus{$group}{$genus}{'AK'}{'FAPM1'}!~/^\s*$/)){
			$fapm1seqsak=$hseqcountrax_genus{$group}{$genus}{'AK'}{'FAPM1'};
		}
		
		my @tabcount=uniq(@{$htyperax_genus{$group}{$genus}});
		my $countuniq=scalar @tabcount;
		
		print OUT $genus."\t".$countuniq."\t".$fapseqshs." (".$fapreadshs.")\t".$fapm1seqshs." (".$fapm1readshs.")\t".$fapseqsak." (".$fapreadsak.")\t".$fapm1seqsak." (".$fapm1readsak.")\n";
	}
	close(OUT);
}
#~ exit;



	###################################################################################################################################################
	##############################################					Previous version					###############################################
	###################################################################################################################################################
	
################################
################################	HSAK primer sep	$hreadsrax_genus{$group}{$genusrax}{$t}{$primer}+=$nbread;
################################
for my $group (sort keys %hreadsmega){
	open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_tissu_BlastN.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_tissu_BlastN.txt\n";
	print OUT "Name\tSpecies\tHS\t\tAK\t\n";
	print OUT "\t\tFAP\tFAPM1\tFAP\tFAPM1\n";
	for my $hpv (sort keys %{$hreadsmega{$group}}){

		my $fapreadshs="-";
		my $fapm1readshs="-";
		my $fapreadsak="-";
		my $fapm1readsak="-";
		
		if((defined($hreadsmega{$group}{$hpv}{'HS'}{'FAP'})) && ($hreadsmega{$group}{$hpv}{'HS'}{'FAP'}!~/^\s*$/)){
			$fapreadshs=$hreadsmega{$group}{$hpv}{'HS'}{'FAP'};
		}
		if((defined($hreadsmega{$group}{$hpv}{'HS'}{'FAPM1'})) && ($hreadsmega{$group}{$hpv}{'HS'}{'FAPM1'}!~/^\s*$/)){
			$fapm1readshs=$hreadsmega{$group}{$hpv}{'HS'}{'FAPM1'};
		}
		
		if((defined($hreadsmega{$group}{$hpv}{'AK'}{'FAP'})) && ($hreadsmega{$group}{$hpv}{'AK'}{'FAP'}!~/^\s*$/)){
			$fapreadsak=$hreadsmega{$group}{$hpv}{'AK'}{'FAP'};
		}
		if((defined($hreadsmega{$group}{$hpv}{'AK'}{'FAPM1'})) && ($hreadsmega{$group}{$hpv}{'AK'}{'FAPM1'}!~/^\s*$/)){
			$fapm1readsak=$hreadsmega{$group}{$hpv}{'AK'}{'FAPM1'};
		}
		
		print OUT $hpv."\t".$hinfomega{$hpv}."\t".$fapreadshs."\t".$fapm1readshs."\t".$fapreadsak."\t".$fapm1readsak."\n";
		print $hpv."\t".$hinfomega{$hpv}."\t".$fapreadshs."\t".$fapm1readshs."\t".$fapreadsak."\t".$fapm1readsak."\n";
		
	}
	close(OUT);
}

for my $group (sort keys %hreadsrax){
	open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_tissu_RaxML.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_tissu_RaxML.txt\n";
	print OUT "Name\tSpecies\tHS\t\tAK\t\n";
	print OUT "\t\tFAP\tFAPM1\tFAP\tFAPM1\n";
	for my $hpv (sort keys %{$hreadsrax{$group}}){
		
		my $fapreadshs="-";
		my $fapm1readshs="-";
		my $fapreadsak="-";
		my $fapm1readsak="-";
		
		if((defined($hreadsrax{$group}{$hpv}{'HS'}{'FAP'})) && ($hreadsrax{$group}{$hpv}{'HS'}{'FAP'}!~/^\s*$/)){
			$fapreadshs=$hreadsrax{$group}{$hpv}{'HS'}{'FAP'};
		}
		if((defined($hreadsrax{$group}{$hpv}{'HS'}{'FAPM1'})) && ($hreadsrax{$group}{$hpv}{'HS'}{'FAPM1'}!~/^\s*$/)){
			$fapm1readshs=$hreadsrax{$group}{$hpv}{'HS'}{'FAPM1'};
		}
		if((defined($hreadsrax{$group}{$hpv}{'AK'}{'FAP'})) && ($hreadsrax{$group}{$hpv}{'AK'}{'FAP'}!~/^\s*$/)){
			$fapreadsak=$hreadsrax{$group}{$hpv}{'AK'}{'FAP'};
		}
		if((defined($hreadsrax{$group}{$hpv}{'AK'}{'FAPM1'})) && ($hreadsrax{$group}{$hpv}{'AK'}{'FAPM1'}!~/^\s*$/)){
			$fapm1readsak=$hreadsrax{$group}{$hpv}{'AK'}{'FAPM1'};
		}
		
		print OUT $hpv."\t".$hinforax{$hpv}."\t".$fapreadshs."\t".$fapm1readshs."\t".$fapreadsak."\t".$fapm1readsak."\n";
	}
	close(OUT);
}

################################
################################	HIV tissu sep
################################
#~ for my $group (sort keys %hreadsmega){
	#~ open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_tissu_BlastN.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_tissu_BlastN.txt\n";
	#~ print OUT "Name\tSpecies\tHIVneg\t\t\tHIVpos\t\t\n";
	#~ print OUT "\t\tCUT\tFAP\tFAPM1\tCUT\tFAP\tFAPM1\n";
	#~ for my $hpv (sort keys %{$hreadsmega{$group}}){

		#~ my $fapreadshivneg="-";
		#~ my $fapm1readshivneg="-";
		#~ my $cutreadshivneg="-";
		
		#~ my $fapreadshivpos="-";
		#~ my $fapm1readshivpos="-";
		#~ my $cutreadshivpos="-";
		
		#~ if((defined($hreadsmega{$group}{$hpv}{'HIVneg'}{'FAP'})) && ($hreadsmega{$group}{$hpv}{'HIVneg'}{'FAP'}!~/^\s*$/)){
			#~ $fapreadshivneg=$hreadsmega{$group}{$hpv}{'HIVneg'}{'FAP'};
		#~ }
		#~ if((defined($hreadsmega{$group}{$hpv}{'HIVneg'}{'FAPM1'})) && ($hreadsmega{$group}{$hpv}{'HIVneg'}{'FAPM1'}!~/^\s*$/)){
			#~ $fapm1readshivneg=$hreadsmega{$group}{$hpv}{'HIVneg'}{'FAPM1'};
		#~ }
		#~ if((defined($hreadsmega{$group}{$hpv}{'HIVneg'}{'CUT'})) && ($hreadsmega{$group}{$hpv}{'HIVneg'}{'CUT'}!~/^\s*$/)){
			#~ $cutreadshivneg=$hreadsmega{$group}{$hpv}{'HIVneg'}{'CUT'};
		#~ }
		
		#~ if((defined($hreadsmega{$group}{$hpv}{'HIVpos'}{'FAP'})) && ($hreadsmega{$group}{$hpv}{'HIVpos'}{'FAP'}!~/^\s*$/)){
			#~ $fapreadshivpos=$hreadsmega{$group}{$hpv}{'HIVpos'}{'FAP'};
		#~ }
		#~ if((defined($hreadsmega{$group}{$hpv}{'HIVpos'}{'FAPM1'})) && ($hreadsmega{$group}{$hpv}{'HIVpos'}{'FAPM1'}!~/^\s*$/)){
			#~ $fapm1readshivpos=$hreadsmega{$group}{$hpv}{'HIVpos'}{'FAPM1'};
		#~ }
		#~ if((defined($hreadsmega{$group}{$hpv}{'HIVpos'}{'CUT'})) && ($hreadsmega{$group}{$hpv}{'HIVpos'}{'CUT'}!~/^\s*$/)){
			#~ $cutreadshivpos=$hreadsmega{$group}{$hpv}{'HIVpos'}{'CUT'};
		#~ }
		
		#~ print OUT $hpv."\t".$hinfomega{$hpv}."\t".$cutreadshivneg."\t".$fapreadshivneg."\t".$fapm1readshivneg."\t".$cutreadshivpos."\t".$fapreadshivpos."\t".$fapm1readshivpos."\n";
	#~ }
	#~ close(OUT);
#~ }

#~ for my $group (sort keys %hreadsrax){
	#~ open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_tissu_RaxML.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_tissu_RaxML.txt\n";
	#~ print OUT "Name\tSpecies\tHIVneg\t\t\tHIVpos\t\t\n";
	#~ print OUT "\t\tCUT\tFAP\tFAPM1\tCUT\tFAP\tFAPM1\n";
	#~ for my $hpv (sort keys %{$hreadsrax{$group}}){
		
		#~ my $fapreadshivneg="-";
		#~ my $fapm1readshivneg="-";
		#~ my $cutreadshivneg="-";
		
		#~ my $fapreadshivpos="-";
		#~ my $fapm1readshivpos="-";
		#~ my $cutreadshivpos="-";
		
		#~ if((defined($hreadsrax{$group}{$hpv}{'HIVneg'}{'FAP'})) && ($hreadsrax{$group}{$hpv}{'HIVneg'}{'FAP'}!~/^\s*$/)){
			#~ $fapreadshivneg=$hreadsrax{$group}{$hpv}{'HIVneg'}{'FAP'};
		#~ }
		#~ if((defined($hreadsrax{$group}{$hpv}{'HIVneg'}{'FAPM1'})) && ($hreadsrax{$group}{$hpv}{'HIVneg'}{'FAPM1'}!~/^\s*$/)){
			#~ $fapm1readshivneg=$hreadsrax{$group}{$hpv}{'HIVneg'}{'FAPM1'};
		#~ }
		#~ if((defined($hreadsrax{$group}{$hpv}{'HIVneg'}{'CUT'})) && ($hreadsrax{$group}{$hpv}{'HIVneg'}{'CUT'}!~/^\s*$/)){
			#~ $cutreadshivneg=$hreadsrax{$group}{$hpv}{'HIVneg'}{'CUT'};
		#~ }
		
		#~ if((defined($hreadsrax{$group}{$hpv}{'HIVpos'}{'FAP'})) && ($hreadsrax{$group}{$hpv}{'HIVpos'}{'FAP'}!~/^\s*$/)){
			#~ $fapreadshivpos=$hreadsrax{$group}{$hpv}{'HIVpos'}{'FAP'};
		#~ }
		#~ if((defined($hreadsrax{$group}{$hpv}{'HIVpos'}{'FAPM1'})) && ($hreadsrax{$group}{$hpv}{'HIVpos'}{'FAPM1'}!~/^\s*$/)){
			#~ $fapm1readshivpos=$hreadsrax{$group}{$hpv}{'HIVpos'}{'FAPM1'};
		#~ }
		#~ if((defined($hreadsrax{$group}{$hpv}{'HIVpos'}{'CUT'})) && ($hreadsrax{$group}{$hpv}{'HIVpos'}{'CUT'}!~/^\s*$/)){
			#~ $cutreadshivpos=$hreadsrax{$group}{$hpv}{'HIVpos'}{'CUT'};
		#~ }
		
		#~ print OUT $hpv."\t".$hinforax{$hpv}."\t".$cutreadshivneg."\t".$fapreadshivneg."\t".$fapm1readshivneg."\t".$cutreadshivpos."\t".$fapreadshivpos."\t".$fapm1readshivpos."\n";
		#~ print $hpv."\t".$hinforax{$hpv}."\t".$cutreadshivneg."\t".$fapreadshivneg."\t".$fapm1readshivneg."\t".$cutreadshivpos."\t".$fapreadshivpos."\t".$fapm1readshivpos."\n";
	#~ }
	#~ close(OUT);
#~ }

#~ for my $group (sort keys %hreadsmega){
	#~ open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_BlastN.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_BlastN.txt\n";
	#~ print OUT "Name\tSpecies\tHIVneg\tHIVpos\n";
	#~ for my $hpv (sort keys %{$hreadsmega{$group}}){

		#~ my $hivnegreads="-";
		#~ my $hivposreads="-";
		
		#~ if((defined($hreadsmega{$group}{$hpv}{'HIVneg'})) && ($hreadsmega{$group}{$hpv}{'HIVneg'}!~/^\s*$/)){
			#~ $hivnegreads=$hreadsmega{$group}{$hpv}{'HIVneg'};
		#~ }
		#~ if((defined($hreadsmega{$group}{$hpv}{'HIVpos'})) && ($hreadsmega{$group}{$hpv}{'HIVpos'}!~/^\s*$/)){
			#~ $hivposreads=$hreadsmega{$group}{$hpv}{'HIVpos'};
		#~ }
		#~ print OUT $hpv."\t".$hinfomega{$hpv}."\t".$hivnegreads."\t".$hivposreads."\n";
	#~ }
	#~ close(OUT);
#~ }

#~ for my $group (sort keys %hreadsrax){
	#~ open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_RaxML.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_RaxML.txt\n";
	#~ print OUT "Name\tSpecies\tHIVneg\tHIVpos\n";
	#~ for my $hpv (sort keys %{$hreadsrax{$group}}){
		#~ my $hivnegreads="-";
		#~ my $hivposreads="-";
		#~ if((defined($hreadsrax{$group}{$hpv}{'HIVneg'})) && ($hreadsrax{$group}{$hpv}{'HIVneg'}!~/^\s*$/)){
			#~ $hivnegreads=$hreadsrax{$group}{$hpv}{'HIVneg'};
		#~ }
		#~ if((defined($hreadsrax{$group}{$hpv}{'HIVpos'})) && ($hreadsrax{$group}{$hpv}{'HIVpos'}!~/^\s*$/)){
			#~ $hivposreads=$hreadsrax{$group}{$hpv}{'HIVpos'};
		#~ }
		#~ print OUT $hpv."\t".$hinforax{$hpv}."\t".$hivnegreads."\t".$hivposreads."\n";
	#~ }
	#~ close(OUT);
#~ }

################################
################################	HSAK tissu sep
################################
#~ for my $group (sort keys %hreadsmega){
	#~ open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_BlastN.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_BlastN.txt\n";
	#~ print OUT "Name\tSpecies\tHS\tAK\n";
	#~ for my $hpv (sort keys %{$hreadsmega{$group}}){
		#~ my $hsreads="-";
		#~ my $akreads="-";
		#~ if((defined($hreadsmega{$group}{$hpv}{'HS'})) && ($hreadsmega{$group}{$hpv}{'HS'}!~/^\s*$/)){
			#~ $hsreads=$hreadsmega{$group}{$hpv}{'HS'};
		#~ }
		#~ if((defined($hreadsmega{$group}{$hpv}{'AK'})) && ($hreadsmega{$group}{$hpv}{'AK'}!~/^\s*$/)){
			#~ $akreads=$hreadsmega{$group}{$hpv}{'AK'};
		#~ }
		#~ print OUT $hpv."\t".$hinfomega{$hpv}."\t".$hsreads."\t".$akreads."\n";
	#~ }
	#~ close(OUT);
#~ }

#~ for my $group (sort keys %hreadsrax){
	#~ open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_RaxML.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_RaxML.txt\n";
	#~ print OUT "Name\tSpecies\tHS\tAK\n";
	#~ for my $hpv (sort keys %{$hreadsrax{$group}}){
		#~ my $hsreads="-";
		#~ my $akreads="-";
		#~ if((defined($hreadsrax{$group}{$hpv}{'HS'})) && ($hreadsrax{$group}{$hpv}{'HS'}!~/^\s*$/)){
			#~ $hsreads=$hreadsrax{$group}{$hpv}{'HS'};
		#~ }
		#~ if((defined($hreadsrax{$group}{$hpv}{'AK'})) && ($hreadsrax{$group}{$hpv}{'AK'}!~/^\s*$/)){
			#~ $akreads=$hreadsrax{$group}{$hpv}{'AK'};
		#~ }
		#~ print OUT $hpv."\t".$hinforax{$hpv}."\t".$hsreads."\t".$akreads."\n";
	#~ }
	#~ close(OUT);
#~ }

################################
################################	HIV primer sep
################################
#~ for my $group (sort keys %hreadsmega){
	#~ open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_BlastN.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_BlastN.txt\n";
	#~ print OUT "Name\tSpecies\tCUT\tFAP\tFAPM1\n";
	#~ for my $hpv (sort keys %{$hreadsmega{$group}}){

		#~ my $fapreads="-";
		#~ my $fapm1reads="-";
		#~ my $cutreads="-";
		
		#~ if((defined($hreadsmega{$group}{$hpv}{'FAP'})) && ($hreadsmega{$group}{$hpv}{'FAP'}!~/^\s*$/)){
			#~ $fapreads=$hreadsmega{$group}{$hpv}{'FAP'};
		#~ }
		#~ if((defined($hreadsmega{$group}{$hpv}{'FAPM1'})) && ($hreadsmega{$group}{$hpv}{'FAPM1'}!~/^\s*$/)){
			#~ $fapm1reads=$hreadsmega{$group}{$hpv}{'FAPM1'};
		#~ }
		#~ if((defined($hreadsmega{$group}{$hpv}{'CUT'})) && ($hreadsmega{$group}{$hpv}{'CUT'}!~/^\s*$/)){
			#~ $cutreads=$hreadsmega{$group}{$hpv}{'CUT'};
		#~ }
		#~ print OUT $hpv."\t".$hinfomega{$hpv}."\t".$cutreads."\t".$fapreads."\t".$fapm1reads."\n";
	#~ }
	#~ close(OUT);
#~ }

#~ for my $group (sort keys %hreadsrax){
	#~ open(OUT, ">".$in_folder."/Table_HPV_presence_".$group."_RaxML.txt") or die "$! : $in_folder/Table_HPV_presence_${group}_RaxML.txt\n";
	#~ print OUT "Name\tSpecies\tCUT\tFAP\tFAPM1\n";
	#~ for my $hpv (sort keys %{$hreadsrax{$group}}){
		#~ my $fapreads="-";
		#~ my $fapm1reads="-";
		#~ my $cutreads="-";
		#~ if((defined($hreadsrax{$group}{$hpv}{'FAP'})) && ($hreadsrax{$group}{$hpv}{'FAP'}!~/^\s*$/)){
			#~ $fapreads=$hreadsrax{$group}{$hpv}{'FAP'};
		#~ }
		#~ if((defined($hreadsrax{$group}{$hpv}{'FAPM1'})) && ($hreadsrax{$group}{$hpv}{'FAPM1'}!~/^\s*$/)){
			#~ $fapm1reads=$hreadsrax{$group}{$hpv}{'FAPM1'};
		#~ }
		#~ if((defined($hreadsrax{$group}{$hpv}{'CUT'})) && ($hreadsrax{$group}{$hpv}{'CUT'}!~/^\s*$/)){
			#~ $cutreads=$hreadsrax{$group}{$hpv}{'CUT'};
		#~ }
		#~ print OUT $hpv."\t".$hinforax{$hpv}."\t".$cutreads."\t".$fapreads."\t".$fapm1reads."\n";
		#~ print $hpv."\t".$hinfomega{$hpv}."\t".$cutreads."\t".$fapreads."\t".$fapm1reads."\n";
	#~ }
	#~ close(OUT);
#~ }

sub getindexmax{
	my $myarray=shift;
	my $index = 0;
	my $maxval = @$myarray[$index];
	
	for(my $i=0; $i<=$#{$myarray}; $i++){
		if($maxval < @$myarray[$i]){
			$maxval = @$myarray[$i];
			$index = $i;
		}
	}
	return($index);
}

