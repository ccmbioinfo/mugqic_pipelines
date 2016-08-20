#!/usr/bin/env perl
# Determine regions of homozygosity from a VCF file with annotated dbSNP and homozygosity.
# The "windowSize" and "maxHetsInWindow" parameters determine where a homozygous regions is called.
use strict;
use warnings;
use Getopt::Long;

sub usage($);
sub getROHSize($$);

(@ARGV == 0) and usage("");
my ($vcfFile, $windowSize, $maxHetsInWindow, $verbose);

GetOptions(	'vcf=s' => \$vcfFile,
			'windowSize=i' => \$windowSize,
			'maxHetsInWindow=i' => \$maxHetsInWindow,
			'v|verbose' => \$verbose) or usage("");

($windowSize) or $windowSize = 25;
(defined $maxHetsInWindow) or $maxHetsInWindow = 2;

my $minReads = 10;
my $minAltReads = 3;
my $minMappingQuality = 50;
my $inROH=0;
my $numHetsInWindow=0;
my $rohStartLoc="";
my $lastChrm="";
my $curChromosome="";
my $chr;
my $start;
my $end;
my $numSNVsInROH=0;
my $totalROHSize=0;
my $totalNonXROHSize=0;
my $candidateROHEndLoc;

# Will hold the last X number of SNV lines where X is the window size
my @SNVArray;

# Will hold "0" or "1" for the last X number of SNVs, indicating whether the SNV is a het
my @homOrHetArray;

# Print a header noting what parameters were used to generate the output
print "#A homozygous region is identified as a genomic window with at least $windowSize total SNVs including no more than $maxHetsInWindow heterozygous SNVs.\n";
print "#Only SNVs with at least $minReads reads and $minMappingQuality mapping quality are considered.\n";
print "#Start pos:\tEnd pos:\t# SNVs\tROH Size\n";


open(VCF_FILE, "<", $vcfFile) or dieScript("Failed to load vcf file: \"${vcfFile}\"");
while (<VCF_FILE>)
{
	(/^#/) and next; #ignore header lines
	(/^\s*$/) and next; #ignore white space lines
	chomp();
	my ($chrm, $pos, $rsID, $ref, $alt, $qual, $filter, $info, @restOfLine) = split(/\t/);
	my ($numAltReads, $numReads, $homozygosity, $mapQ);
	
	# Ignore indels
	(length($ref) != length($alt)) and next;
	
	($info =~ /ALTC=([^;\t]+)/)		and $numAltReads = $1;
	($info =~ /RDC=([^;\t]+)/)		and $numReads = $1;
	($info =~ /HMZ=([^;\t]+)/)		and $homozygosity = $1;
	($info =~ /MQ=([^;\t]+)/)		and $mapQ = $1;

	(!$numReads or $numReads < $minReads) and next;	#only consider SNVs with sufficient coverage
	($mapQ < $minMappingQuality) and next;
	($homozygosity =~ /possibly/) and next;

	# IF annotated as a het, require at least 3 alt reads to count it
	if ($homozygosity =~ /het/)
	{
		(!$numAltReads or $numAltReads < $minAltReads) and next;
	}
	
	if ($lastChrm && $chrm ne $lastChrm)
	{
		# Any existing ROH is ended at the end of the chromosome
		if ($inROH)
		{
			$inROH = 0;
			my @lastSNVLine = split(/\t/, $SNVArray[$#SNVArray]);
			($chr, $start) = @lastSNVLine[0..1];
			my $rohEndLoc = $chr.":".$start;
			#print STDERR "End of chromosome ending a ROH. Last SNV: ".$rohEndLoc."\n";
			my $rohSize = getROHSize($rohStartLoc, $rohEndLoc);
			print $rohStartLoc."\t".$rohEndLoc."\t".$numSNVsInROH."\t".$rohSize."\n";	#print the info on the ROH that just finished
			
			$totalROHSize += $rohSize;
			my ($chrm1, $pos1) = split(/:/, $rohStartLoc);
			if (lc($chrm1) !~ /chrx|chry/)
			{
				$totalNonXROHSize += $rohSize;
			}
		}# else {
			#print STDERR "End of chromosome, not in a ROH. Last SNV: ".$rohEndLoc."\n";
		#}
		
		@SNVArray = ();
		@homOrHetArray = ();
		$rohStartLoc = "";
		$numHetsInWindow = 0;
		$numSNVsInROH = 0;
	}
	
	#Add the new SNV line to the current window
	push(@SNVArray, $_);
	
	my $curSNVIsHet = ($homozygosity =~ /het/) ? 1 : 0;
	push(@homOrHetArray, $curSNVIsHet);
	$numHetsInWindow += $curSNVIsHet;
	
	if (scalar(@SNVArray) > $windowSize)
	{
		my $removedSNVIsHet = $homOrHetArray[0];
		shift(@homOrHetArray);
		shift(@SNVArray);
		$numHetsInWindow -= $removedSNVIsHet;
	}
	
	if ($curSNVIsHet)
	{
		if (!$candidateROHEndLoc)
		{
			$candidateROHEndLoc = $chrm.":".$pos;
		}
		if ($numHetsInWindow > $maxHetsInWindow)
		{
			if ($inROH)
			{
				$inROH = 0;
				# The last SNV in the ROH is the new element that was just pushed on.
				my $rohEndLoc = $candidateROHEndLoc;
				#print STDERR "Het ending a ROH:".$rohEndLoc."\n";
				my $rohSize = getROHSize($rohStartLoc, $rohEndLoc);
				print $rohStartLoc."\t".$rohEndLoc."\t".$numSNVsInROH."\t".$rohSize."\n";	#print the info on the ROH that just finished
				
				$totalROHSize += $rohSize;
				my ($chrm1, $pos1) = split(/:/, $rohStartLoc);
				if (lc($chrm1) !~ /chrx|chry/)
				{
					$totalNonXROHSize += $rohSize;
				}
				@SNVArray = ();
				@homOrHetArray = ();
				$rohStartLoc = "";
				$numHetsInWindow = 0;
				$numSNVsInROH = 0;
			}
		}
	}
	else
	{
		if (scalar(@SNVArray) >= $windowSize && $numHetsInWindow <= $maxHetsInWindow)
		{
			# We meet the criteria to be in a ROH
			if (!$inROH)
			{
				$inROH = 1;
				# We are in a new ROH. The start of the ROH is the first SNV in our list.
				my @startSNVLine = split(/\t/, $SNVArray[0]);
				my ($chr, $start) = @startSNVLine[0..1];
				$rohStartLoc = $chr.":".$start;
				
				$numSNVsInROH = $windowSize-1;
			}
		}
		$candidateROHEndLoc = "";
	}
	$lastChrm = $chrm;
	$numSNVsInROH++;
}
close(VCF_FILE);

print "#Total ROH Size\t$totalROHSize\n";
print "#Total autosomal ROH Size\t$totalNonXROHSize\n";


###############################################################################
sub getROHSize($$)
{
	my @roh1 = split(/:/, $_[0]);
	my @roh2 = split(/:/, $_[1]);
	return $roh2[1] - $roh1[1];
}

###############################################################################

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;
This program takes a VCF file as input and assumes that it has been annotated with dbSNV in
the ID column, and annotated for homozygosity in the INFO field (HMZ=het|multiple het|hom|possibly hom).
It outputs a list of chromosome ranges that are predicted as regions of homozygosity for the sample.

Usage:
$0 --sampleName NAME --sampleConfig FILE [options]

OPTIONS:
--vcf FILE              VCF file to use as input
--windowSize INT        There must be a run of 25 homozygous SNVs with no more than maxHetsInWindow het
--maxHetsInWindow INT   variants to call a region of homozygosity.
--verbose               Writes additional information to STDERR

EOF
    exit;
}
