#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my ($ccdsGenesPath, $sampleCoveragePath, $printIntervals, $verbose);

GetOptions(	'ccdsGenes=s' => \$ccdsGenesPath,
			'sample=s' => \$sampleCoveragePath,
			'printIntervals' => \$printIntervals,
			'v' => \$verbose) or &usage();

$ccdsGenesPath or &usage("ERROR: Missing ccdsGenes argument, which should be the main CCDS file.\n");
$sampleCoveragePath or &usage("ERROR: Missing sample argument, which is the coverage file for a sample.\n");

open(CCDS_GENES, "<$ccdsGenesPath") or die "Failed to open ccdsGenes file: $!";
open(SAMPLE_COVERAGE, "<$sampleCoveragePath") or die "Failed to open sample coverage file: $!";

# First load entire sample coverage file into a hash table.
my $includeChr = 0;
my %sampleCoverage;
while (<SAMPLE_COVERAGE>)
{
	/^Target/ and next;
	chomp();
	my @line = split(/\t/);
	$sampleCoverage{$line[0]} = $_;
	($line[0] =~ /chr/) and $includeChr = 1;
}
close(SAMPLE_COVERAGE);

# Go through each line in the CCDS genes file and sum coverage for the associated intervals
while (<CCDS_GENES>)
{
	/^#/ and next;
	chomp();
	my ($chrm, $refseqID, $geneSymbol, $geneID, $ccdsID, $ccdsStatus, $strand, $cdsFrom, $cdsTo, $cdsLocs) = split(/\t/);
	($ccdsStatus =~ /Withdrawn/ || $cdsLocs !~ /\[(.*)\]/) and next;
	
	my $matchLocs = $cdsLocs;
	if ($cdsLocs =~ /\[(.*)\]/)
	{
		$matchLocs = $1;
	}
	#print STDERR $cdsLocs."\t".$matchLocs."\n";
	my @cdsIntervals = split(/, /, $matchLocs);
	
	my $totalBases = 0;
	my $totalCoverage = 0;
	my @covgAboveXArray = (0, 0, 0, 0, 0, 0, 0);
	for (my $i=0; $i < scalar(@cdsIntervals); $i++)
	{
		my ($start, $end) = split(/-/, $cdsIntervals[$i]);
		
		# I don't know why by after processing by GATK all intervals have their start
		# indices incremented by one, though the interval end index is not incremented.
		$start = $start + 1;
		$end = $end + 1;
		my $key = "";
		($includeChr) and $key = "chr";
		$key .= $chrm.":".$start."-".$end;
		
		#print STDERR $key."\n";
		my $intervalLine = $sampleCoverage{$key};
		if (!exists $sampleCoverage{$key} or !$intervalLine)
		{
			$verbose and print STDERR "Interval $key not found for CCDS ID: $ccdsID\n";
			next;
		}
		my @intervalData = split(/\t/, $intervalLine);
		
		my $intervalSize = $end - $start + 1;
		$totalBases += $intervalSize;
		$totalCoverage += $intervalData[1];
		# The "%_above_X coverage" columns start at 8 (zero-based)
		for (my $j=8; $j < scalar(@intervalData); $j++)
		{
			$covgAboveXArray[$j-8] += $intervalSize * $intervalData[$j];
		}
	}
	
	if ($totalBases <= 0)
	{
		$verbose and print STDERR "No intervals found for CCDS ID: $ccdsID\n";
		next;
	}
	
	# Change covgAboveXArray to percentages.
	for (my $j=0; $j < scalar(@covgAboveXArray); $j++)
	{
		$covgAboveXArray[$j] = $covgAboveXArray[$j] / $totalBases;
	}
	
	# First print out the totals for the gene
	print join("\t", "", $geneSymbol, $refseqID, $ccdsID, $totalCoverage, $totalCoverage / $totalBases, @covgAboveXArray)."\n";
	
	if ($printIntervals)
	{
		# Go through all intervals again and print out the interval data along with columns identifying the gene
		for (my $i=0; $i < scalar(@cdsIntervals); $i++)
		{
			my ($start, $end) = split(/-/, $cdsIntervals[$i]);
			$start = $start + 1;
			$end = $end + 1;
			my $key = "";
			($includeChr) and $key = "chr";
			$key .= $chrm.":".$start."-".$end;
			my $intervalLine = $sampleCoverage{$key};
			if (!$intervalLine)
			{
				$verbose and print STDERR "Interval $key not found for CCDS ID: $ccdsID\n";
				next;
			}
			my @intervalData = split(/\t/, $intervalLine);
			
			print join("\t", $intervalData[0], $geneSymbol, $refseqID, $ccdsID, $intervalData[1], $intervalData[2], @intervalData[8 .. $#intervalData])."\n";
		}
	}
}
close(CCDS_GENES);



###############################################################################

sub usage() 
{
	($_[0]) and print $_[0]."\n";
    print STDERR <<EOF;

Usage:
$0 --ccdsGenes pathToCCDSFile --ccdsKgMap pathToKgMap --knownToRefSeq pathToKgRefSeqTable --refSeqGenes pathToRefSeqTable [-v] > newCCDSFile.txt

Uses the CCDS ID to find the RefSeq gene name associated with it, and adds
this name to the CCDS genes list in column 13, which is normally blank in
CCDS.

EOF

    exit;
}

