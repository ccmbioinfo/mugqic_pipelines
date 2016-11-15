#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;
  

my ($read_files)=@ARGV;


my $Datapath="~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath="";


my $IDfile1=$Workpath.$read_files."_contigs_micro_cds_bwa_hitsID.txt";
my $IDfile2=$Workpath.$read_files."_contigs_n_micro_cds_blat_hitsID.txt";
my $IDfile3=$Workpath.$read_files."_singletons_micro_cds_bwa_hitsID.txt";
my $IDfile4=$Workpath.$read_files."_singletons_n_micro_cds_blat_hitsID.txt";

my $tmpfile=$Workpath.$read_files."_bwablat_hitsID.txt";	
system("cat $IDfile1 $IDfile2 $IDfile3 $IDfile4 > $tmpfile");


my $infile=$Workpath.$read_files."_bwablat_hitsID_unique.txt";	
my $outfile=$Workpath."microbial_cds_sub.fasta";
print "Outputfiles: $outfile\n\n";


my %hits;
unlink($infile); 
open(OUTPUT,'>>', $infile) or die "Error opening $infile : $!\n";
open(INPUT, $tmpfile) or die "Error opening $tmpfile : $!\n";
while  ( my $line = <INPUT>) {
	chomp($line);
	if (not exists $hits{$line}) {
		print OUTPUT "$line\n";
	}
	$hits{$line} = 1;
}
close INPUT; 
close OUTPUT; 
print "number of mapped genes:  " . keys( %hits ) . ".\n\n\n";



system("xargs samtools faidx \$BLASTDB/microbial_all_cds.fasta < $infile > ".$outfile); 



exit 0;