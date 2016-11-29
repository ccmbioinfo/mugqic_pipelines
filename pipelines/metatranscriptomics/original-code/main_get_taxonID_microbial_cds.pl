#!/usr/local/bin/perl -w
  

my ($read_files)=@ARGV;


my $Datapath="~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath="";


my $infile0=$ENV{'BLASTDB'}."/microbial_all_cds_IDs_map_taxid.txt";
my $infile=$Workpath."microbial_cds_sub_IDs_length.txt";
my $outfile=$Workpath."microbial_cds_sub_IDs_map_taxid.txt";
print "Outputfile: $outfile\n\n";


my %refseq;
open(INPUT0, $infile0) or die "Error opening $infile0 : $!\n";
while( my $line1 = <INPUT0> ) {
	chomp($line1);	
	my @line=split(/\t/,$line1);	
	$refseq{$line[0]} = $line1;
}
close INPUT0;
print "number of refseqIDs:  " . keys( %refseq ) . ".\n";


my %geneIDs;
my %taxonIDs;
unlink($outfile); 
open(OUTPUT,'>>', $outfile) or die "Error opening $outfile : $!\n";
open(INPUT, $infile) or die "Error opening $infile : $!\n";
while( my $line1 = <INPUT> ) {
	chomp($line1);	
	my @line=split(/\t/,$line1);

	my $myrefseq1=$line[0];
	$geneIDs{$myrefseq1}=1;

	if (not exists $refseq{$myrefseq1}) {
		print "$line[0]\n$myrefseq1\n";
		exit 1;
	}else{
		my @tmp=split(/\t/,$refseq{$myrefseq1});
		if (not exists $taxonIDs{$tmp[3]}) {
			$taxonIDs{$tmp[3]}=1;
		}
		$taxonIDs{$tmp[3]}+=1;

		print OUTPUT "$refseq{$myrefseq1}\n";
	}
}
close INPUT;
close OUTPUT;
print "number of geneIDs:  " . keys( %geneIDs ) . ".\n";
print "number of taxonIDs:  " . keys( %taxonIDs ) . ".\n\n\n";


exit 0;

