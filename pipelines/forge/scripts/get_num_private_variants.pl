#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);

(@ARGV == 0) and usage("");

my ($vcfFile, $name, $help, $verbose);
GetOptions(	'vcf=s' => \$vcfFile,
			'name=s' => \$name,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");
			
open(VCF_FILE, "<".$vcfFile) or die("Failed to open variants input file: $vcfFile\n");

my $numPrivateVariants = 0;
my $numPrivateStopVariants = 0;

while (<VCF_FILE>)
{
	(/^#/) and next; # Ignore header lines
	chomp();
	my @vcfLine = split(/\t/);
	my $info = $vcfLine[7];
	my ($variantType) = ($info =~ /VT=([^;\t]+)/);
	(!defined $variantType) and next;
	
	if ($info =~ /EVSMAF=([^;\t]+)/)
	{
		($1 > 0.0) and next;
	}
	elsif ($info =~ /THGMAF=([^;\t]+)/)
	{
		($1 > 0.0) and next;
	}
	elsif ($vcfLine[2] !~ /\./)
	{
		next;
	}
	elsif ($info =~ /PSN=([^;\t]+)/)
	{
		($1 > 1) and next;
	}

	# We passed all the checks - the variant has not been seen before.
	$numPrivateVariants++;
	if ($variantType =~ /^frameshift|^stopgain|^splicing$/)
	{
		#print STDERR $_."\n";
		$numPrivateStopVariants++;
	}
}
close(VCF_FILE);

$name and print $name."\t";
print $numPrivateVariants."\t".$numPrivateStopVariants."\n";


###############################################################################

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;
This program takes an input vcf file that is annotated with fields PSN (number
of previously seen samples), VT (variant type), and THGMAF and EVS MAF for the
minor allele frequencies from 1000 genomes and the exome variant server, and
returns the number of variants found in the sample but never seen in any other
samples.

Usage:
$0 --vcf FILE

OPTIONS:
--vcf       FILE  VCF file to determine number of private mutations for
--help|h          Prints usage information

EOF
    exit 1;
}
