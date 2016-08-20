#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);

(@ARGV == 0) and usage("");

my ($rohFile, $vcfFile, $help, $verbose);
GetOptions(	'rohFile=s' => \$rohFile,
			'vcf=s', \$vcfFile,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");
			
open(REGIONS_OF_HOMOZYGOSITY_FILE, "<".$rohFile) or die("Failed to open region of homozygosity file: $rohFile\n");
open(INPUT, "<".$vcfFile) or die("Failed to open vcf input file: $vcfFile\n");

my @rohArray; #regions of homozygosity array
my @line;

# Read in entire ROH file
while (<REGIONS_OF_HOMOZYGOSITY_FILE>)
{
	chomp();
	(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
	(/^\s*$/) and next; #ignore white space lines
	@line = split(/\t/);
	
	# Verify that the ROH coordinates are on the same chromosome
	my ($chrm1, $pos1) = split(/:/, $line[0]);
	my ($chrm2, $pos2) = split(/:/, $line[1]);
	if ($chrm1 ne $chrm2)
	{
		print STDERR "Region of homozygosity: chromosome \"$chrm1\" doesn't match \"$chrm2\"";
		next;
	}
	push(@rohArray, lc($chrm1)."\t".$pos1."\t".lc($chrm2)."\t".$pos2);
}
close(REGIONS_OF_HOMOZYGOSITY_FILE);


while (<INPUT>)
{
	if (/^#/)
	{
		# This is a header line, so just print it.
		print;
		next;
	}
	(/^\s*$/) and next; #ignore white space lines
	chomp();
	my ($chrm, $pos, $rsID, $ref, $alt, $qual, $filter, $info, @restOfLine) = split(/\t/, $_, -1);
	
	my $inROH = 0;
	foreach my $roh (@rohArray)
	{
		my ($chrm1, $rohPosStart, $chrm2, $rohPosEnd) = split(/\t/, $roh);
		if ((lc($chrm) eq $chrm1) and ($pos >= $rohPosStart && $pos < $rohPosEnd))
		{
			$inROH = 1;
			last;
		}
	}
	$info .= ";ROH=$inROH";
	print join("\t", $chrm, $pos, $rsID, $ref, $alt, $qual, $filter, $info, @restOfLine)."\n";
}
close(INPUT);


###############################################################################

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;
This program takes a file listing regions of homozygosity (ROH) and an input file of variants
and outputs the variants to STDOUT with an annotation of 1 or 0 for in or out of the ROH.

Usage:
$0 --rohFile FILE --input FILE > output.txt

OPTIONS:
--rohFile  FILE  File listing regions of homozygosity
--input    FILE  File of variants to annotate is in or out of ROH
--help|h         Prints usage information

EOF
    exit 1;
}
