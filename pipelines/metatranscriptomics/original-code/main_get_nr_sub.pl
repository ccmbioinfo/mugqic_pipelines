#!/usr/local/bin/perl -w
  

my ($read_files)=@ARGV;


my $Datapath="~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath="";


my $outfile1=$Workpath."nr_all_sub.fasta";
my $outfile2=$Workpath."nr_all_sub_IDs_giID.txt";
print "Outputfiles: $outfile1\n$outfile2\n\n";


my $IDfile1=$Workpath.$read_files."_contigs_nr_diamond_hitsID_sub.txt";
my $IDfile2=$Workpath.$read_files."_singletons_nr_diamond_hitsID_sub.txt";
my $tmpfile=$Workpath.$read_files."_diamond_histID.txt";	
system("cat $IDfile1 $IDfile2 > $tmpfile");

my %hits;
open(INPUT, $tmpfile) or die "Error opening $tmpfile : $!\n";
while  ( my $line = <INPUT>) {
	chomp($line);
	my $hitgiID=$line;
	$hitgiID =~ s/^gi\|//; 
	$hitgiID =~ s/\|([a-zA-Z]*)\|([a-zA-Z0-9._]*)\|([a-zA-Z0-9._]*)$//; 
	$hits{$hitgiID} = 1;
}
close INPUT; 
print "number of mapped proteins:  " . keys( %hits ) . ".\n";


unlink($outfile2); 
open(OUTPUT2,'>>', $outfile2) or die "Error opening $outfile2 : $!\n";
while ( my ($key, $value) = each(%hits) ) {
	print OUTPUT2 "$key\n";
}
close OUTPUT2;


system("blastdbcmd -db \$BLASTDB/nr -dbtype prot -entry_batch ".$outfile2." -outfmt %f -target_only -out ".$outfile1);


exit 0;