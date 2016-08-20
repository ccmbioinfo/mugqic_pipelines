#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);

(@ARGV == 0) and usage("");

my ($vcfFile, $inputDelim, $missingVal, $help, $verbose);
my @cols;
my @headers;
my @escapeCols;

GetOptions(	'vcf=s' => \$vcfFile,
			'delim=s' => \$inputDelim,
			'missingVal=s' => \$missingVal,
			'cols=s' => \@cols,
			'escapeCols=s' => \@escapeCols,
			'headers=s' => \@headers,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

(!$inputDelim) and $inputDelim = "\t";
(!$missingVal) and $missingVal = "";

# Allow either multiple --cols parameters, or tab-separated cols (e.g. --cols CHROM	POS	ID)
@cols = split(/,/, join(',', @cols));
if (scalar(@headers) > 0)
{
	@headers = split(/,,/, join(',,', @headers));
}
@escapeCols = split(/,/, join(',', @escapeCols));
my %escapeColsHash = map { $_ => 1 } @escapeCols;

my $numCols = scalar(@cols);
($numCols == 0) and usage("Missing required parameter: cols. You must specify the columns to include in the output.");

my $numHeaders = scalar(@headers);
if ($numHeaders > 0 and $numHeaders != $numCols)
{
	usage("Number of header strings specified ($numHeaders) is different than number of columns specified ($numCols). There should be the same number.\n");
}


open(VCF_FILE, "<".$vcfFile) or die("Failed to open variants input file $vcfFile: $!\n");

#Print the header line
if ($numHeaders > 0)
{
	print join("\t", @headers)."\n";
}

my %basicColsHash = (
	"CHROM" => 0,
	"POS" => 1,
	"ID" => 2,
	"REF" => 3,
	"ALT" => 4,
	"QUAL" => 5,
	"FILTER" => 6
);

while (<VCF_FILE>)
{
	# If this is a header line, skip it. We printed our own header already.
	(/^#/) and next;
	chomp();
	my @vcfLine = split(/\t/, $_, -1);
	my $outputLine = "";
	my $curDelim = ""; #delim starts off empty so that the line doesn't start with a delimiter
	foreach my $colName (@cols)
	{
		# First try to handle standard VCF columns
		if (exists $basicColsHash{$colName})
		{
			my $vcfColIndex = $basicColsHash{$colName};
			$outputLine .= $curDelim.$vcfLine[$vcfColIndex];
		}
		elsif ($colName eq "CHRPOS")
		{
			# Special column name that tells us to combine Chrom and Pos
			$outputLine .= $curDelim.$vcfLine[0].":".$vcfLine[1];
		}
		else
		{
			# This is a custom column - we find the value in the VCF INFO column
			my $val = $missingVal;
			if ($vcfLine[7] =~ /$colName=([^;\t]+)/)
			{
				$val = $1;
			
				# If col is flagged to be "escaped", then add a space before the value.
				# When opening a file in Excel this will prevent it from interpreting certain
				# text strings as numbers.
				(exists $escapeColsHash{$colName}) and $val = " ".$val;
			}
			$outputLine .= $curDelim.$val;
		}
		$curDelim = $inputDelim;
	}
	print $outputLine."\n";
}
close(VCF_FILE);


###############################################################################

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;
This program takes an input VCF file and converts it to a delimited format by extracting
fields from the INFO field of the VCF. Parameters specify what name the custom fields
have and what heading to give them in the output. The basic VCF columns that can be included
are CHROM, POS, ID, REF, ALT, QUAL, FILTER. Any other custom columns can be extracted from
the INFO field assuming the format "COLID=<value>". CHRMPOS is a special column name meaning
that the CHROM and POS cols should be combined together.

Usage:
perl vcf2columns.pl --vcf FILE --cols "colID1,colID2,colID3..." > output.vcf

OPTIONS:
--vcf        FILE  VCF file to split into delimited columns
--delim      CHAR  Delimiter to use to separate values (DEFAULT: tab)
--missingVal STR   Delimiter to use to separate values (DEFAULT: empty)
--cols       STR   Comma-separated list of column IDs to include.
--headers    STR   Double-comma-separated list of column headings to output. (DEFAULT: do not output)
--help|h           Prints usage information

EOF
    exit 1;
}
