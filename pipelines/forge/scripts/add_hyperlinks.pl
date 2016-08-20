#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);

(@ARGV == 0) and usage("");

my ($inputFile, $forExcel, $help, $verbose);
GetOptions(	'input=s' => \$inputFile,
			'excel' => \$forExcel,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

$help and usage("");
$inputFile or usage("Missing required parameter: exonvars. You must specify the annovar exonic variant function output file");

open(INPUT_FILE, "<".$inputFile) or die("Failed to open input file: $inputFile");

(defined $forExcel) or $forExcel = 1;

while (<INPUT_FILE>)
{
	chomp();
	(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
	(/^\s*$/) and next; #ignore white space lines
	if (/^\s*Position/)
	{
		my @headerLine = split(/\t/, $_);
		print join("\t", $headerLine[0], "UCSC Link", @headerLine[1..$#headerLine])."\n";
		next;
	}
	my @line = split(/\t/, $_, -1);
	
	my $hyperlink = "http://genome.ucsc.edu/cgi-bin/hgTracks?hgt_doJsCommand=&hgt.out3=10x&position=".$line[0];
	my $UCSClink = $forExcel ? "=HYPERLINK(\"".$hyperlink."\",\"UCSC link\")" : $hyperlink;
	
	print join("\t", $line[0], $UCSClink, @line[1..$#line])."\n";
}
close(INPUT_FILE);


###############################################################################

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;
Takes a tab-delimited file of variants as input. Reads the first column assuming
it is in the format chr1:2234567, and adds a new column as column 2 which is a
hyperlink to the UCSC browser at the variant position. If parameter --excel
is specified, the new column is formatted to be an automatic link in Excel.

Usage:
perl add_hyperlinks.pl --input FILE > output.for_excel.txt

OPTIONS:
--input FILE    Tab-delimited input file
--excel         Delimiter to use to separate values (DEFAULT: tab)
--help|h        Prints usage information

EOF
    exit 1;
}
