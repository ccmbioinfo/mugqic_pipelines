#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;
  

my ($read_files, $data_type)=@ARGV;
#my $read_files="cow";
#my $data_type="assembly";


my $Datapath="~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath="";



my $IDfile;
my $IDfile2;
my $IDfile21;
my $infile;
my $outfile;
my $outfile1;
if ($data_type eq 'assembly') {
	$IDfile=$Workpath.$read_files."_trinity_bwa_pairs.txt";
	$infile=$Workpath.$read_files."_contigs.fasta";

	$outfile=$Workpath.$read_files."_contigs_IDs.txt";
	$outfile1=$Workpath.$read_files."_contigs_IDs_length.txt";
	print "Outputfile: $outfile1\n\n";
}elsif ($data_type eq 'bwa') {
	$IDfile=$Workpath.$read_files."_contigs_IDs.txt";
	$IDfile2=$Workpath.$read_files."_contigs_micro_cds_bwa_IDs.txt"; 
	$IDfile21=$Workpath.$read_files."_singletons_micro_cds_bwa_IDs.txt"; 
}elsif ($data_type eq 'blat') {
	$IDfile=$Workpath.$read_files."_contigs_IDs.txt";
	$IDfile2=$Workpath.$read_files."_contigs_n_micro_cds_blat_IDs.txt"; 
	$IDfile21=$Workpath.$read_files."_singletons_n_micro_cds_blat_IDs.txt"; 
}elsif ($data_type eq 'diamond') {
	$IDfile=$Workpath.$read_files."_contigs_IDs.txt";
	$IDfile2=$Workpath.$read_files."_contigs_nr_diamond_IDs.txt"; 
	$IDfile21=$Workpath.$read_files."_singletons_nr_diamond_IDs.txt"; 
}




if ($data_type eq 'assembly') {
	my %IDs;
	open(INPUT0, $IDfile) or die "Error opening $IDfile : $!\n";
	while  ( my $line = <INPUT0>) {
		chomp($line);
		my @line1 = split(/\t/, $line);
		if (exists $IDs{$line1[1]}) {
			$IDs{$line1[1]} ++;
		}else{
			$IDs{$line1[1]} = 1;
		}
	}
	close INPUT0; 


	my %contigs;
	my $total=0;
	unlink($outfile); 
	open(OUTPUT,'>>', $outfile) or die "Error opening $outfile : $!\n";
	unlink($outfile1); 
	open(OUTPUT1,'>>', $outfile1) or die "Error opening $outfile1 : $!\n";
	my $seqIO=Bio::SeqIO->new(-file=>$infile, -format=>'fasta');
	while (my $seqobj=$seqIO->next_seq) {
		my $head=$seqobj->display_id;
		$contigs{$head}=$seqobj->length;
		
		if (not exists $IDs{$head}) {
			$IDs{$head}=1;
		}
		$total += $IDs{$head};
		print OUTPUT "$head\t$IDs{$head}\n";
		print OUTPUT1 "$head\t$IDs{$head}\t$contigs{$head}\n";
	}
	close OUTPUT1;
	close OUTPUT;
	print "number of all contigs:  " . keys( %contigs ) . ".\n";
	print "total number of reads mapping to contigs = $total.\n\n\n";
}else{
	my %IDs;
	open(INPUT0, $IDfile) or die "Error opening $IDfile : $!\n";
	while  ( my $line = <INPUT0>) {
		chomp($line);
		my @line1 = split(/\t/, $line);
		$IDs{$line1[0]} = $line1[1];
	}
	close INPUT0; 


	my $total=0;
	open(INPUT2, $IDfile2) or die "Error opening $IDfile2 : $!\n";
	while  ( my $line = <INPUT2>) {
		chomp($line);
		if (exists $IDs{$line}) {
			$total += $IDs{$line};
		}else{
			$total += 1;
		}
	}
	close INPUT2; 
	open(INPUT2, $IDfile21) or die "Error opening $IDfile21 : $!\n";
	while  ( my $line = <INPUT2>) {
		chomp($line);
		if (exists $IDs{$line}) {
			$total += $IDs{$line};
		}else{
			$total += 1;
		}
	}
	close INPUT2; 
	print "Total number of mapped-reads = $total.\n\n\n";
}


exit 0;
