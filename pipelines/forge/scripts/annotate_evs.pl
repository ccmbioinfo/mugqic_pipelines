#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);

(@ARGV == 0) and usage("");

my ($vcfFile, $evsdb, $help, $verbose);
GetOptions(	'vcf=s' => \$vcfFile,
			'evsdb=s' => \$evsdb,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");
			
open(EBS_DB, "<".$evsdb) or die("Failed to open EVS database file: $evsdb\n");
open(VCF_FILE, "<".$vcfFile) or die("Failed to open variants input file: $vcfFile\n");

my %evsHash;

# Read in entire EVS file
while (<EBS_DB>)
{
	chomp();
	(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
	(/^\s*$/) and next; #ignore white space lines
	my @line = split(/\t/);
	$evsHash{"chr".lc($line[0]).":".$line[1]} = $_;
}
close(EBS_DB);


while (<VCF_FILE>)
{
	if (/^#/)
	{
		# This is a header line, so just print it.
		print;
		next;
	}
	chomp();
	my @vcfLine = split(/\t/, $_, -1);

	my ($MAF, $GTC, $readDepth, $GERP, $grantham);
	
	my $key = lc($vcfLine[0]).":".$vcfLine[1];
	if (exists $evsHash{$key})
	{
		my @evsLine = split(/\t/, $evsHash{$key});
		
		# Scale the MAF to be a fraction rather than a %, to be consistent with thousand genomes
		($evsLine[7] =~ /;MAF=([^,\t]+),([^,\t]+),([^,\t]+);/) and $MAF = $3 / 100.0;
		($evsLine[7] =~ /;DP=([^;\t]+);/) and $readDepth = $1;
		($evsLine[7] =~ /;CG=([^;\t]+);/) and $GERP = $1;
		if ($evsLine[7] =~ /;GS=([^;\t]+);/)
		{
			($grantham) = split(/,/, $1);
			($grantham eq "NA") and $grantham = "";
		}
		if ($evsLine[7] =~ /;GTC=([^;\t]+);/)
		{
			$GTC = $1;
			$GTC =~ s/,/:/g;
		}
	}
	(!defined $MAF) and $MAF = 0;
	$vcfLine[7] .= ";EVSMAF=$MAF";
	($GTC) and $vcfLine[7] .= ";EVSGTC=$GTC";
	($readDepth) and $vcfLine[7] .= ";EVSRD=$readDepth";
	($GERP) and $vcfLine[7] .= ";EVSGERP=$GERP";
	($grantham) and $vcfLine[7] .= ";EVSGRA=$grantham";
	print join("\t", @vcfLine)."\n";
}
close(VCF_FILE);


###############################################################################

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;
This program takes an input vcf file and a VCF file of the exome variant server (EVS)
variant frequencies and adds a number of fields to the VCF INFO field of the variant.
Fields added:
EVSMAF   EVS variant minor allele frequency
EVSGTC   EVS sample genotype counts (#homozygous alt, #het alt, #hom ref)
EVSRD    EVS average sample read depth
EVSGERP  EVS reported GERP score
EVSGRA   EVS reported grantham score

Usage:
$0 --vcf FILE --evsdb FILE > output.vcf

OPTIONS:
--vcf       FILE  VCF file to annotate with EVS info
--evsdb     FILE  VCF file with info on exome variant server variant frequencies
--help|h          Prints usage information

EOF
    exit 1;
}
