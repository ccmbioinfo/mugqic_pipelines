#!/usr/local/bin/perl -w
  

my ($read_files)=@ARGV;
#my $read_files='cow';


my $Datapath="~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath="";



my $infile1=$Workpath."microbial_cds_sub_table_counts.txt";	##[geneID, length, taxonID, specie, phylum, #reads b#]
my $infile2=$Workpath."nr_all_sub_table_counts.txt";	##[proteinID, length, taxonID, specie, phylum, #reads, b#]
my $outfile1=$Workpath.$read_files."_table_counts_all.txt";	##[geneID/proteinID, length, taxonID, specie, phylum, #reads, b#]
my $outfile2=$Workpath.$read_files."_table_RPKM_all.txt";	##[geneID/proteinID, length, taxonID, specie, phylum, RPKM, b#]
print "Outputfile: $outfile2\n\n";


system("cat $infile1 $infile2 > $outfile1");

my $N = 1;
my @N_total=(0)x $N;


my %genes;
open(INPUT, $outfile1) or die "Error opening $outfile1 : $!\n";
while  ( my $line1 = <INPUT>) {
	chomp($line1);
	my @line=split(/\t/,$line1);
	$genes{$line[0]}=$line1;
	for (my $i=0; $i<$N; $i++) {
		$N_total[$i] += $line[5+$i];
	}
}
close INPUT; 
### genes=[geneID/proteinID, length, taxonID, specie, phylum, #reads, b#]
print "number of genes/proteins:  " . keys( %genes ) . ".\n";

print "total number of reads:  @N_total.\n\n\n";



unlink($outfile2); 
open(OUTPUT2,'>>', $outfile2) or die "Error opening $outfile2 : $!\n";
foreach $key (sort keys %genes){
	my @line=split(/\t/,$genes{$key});

	my @line_rpkm;
	for (my $i=0; $i<$N; $i++) {
		$line_rpkm[$i]=(10**9)*$line[5+$i]/($N_total[$i]*$line[1]);
		$line_rpkm[$i]=sprintf("%.4f", $line_rpkm[$i]);
	}
	print OUTPUT2 join("\t",@line[0..1],@line[5..(5+$N-1)],@line[2..4],@line_rpkm,$line[-1])."\n";
}
close OUTPUT2;
