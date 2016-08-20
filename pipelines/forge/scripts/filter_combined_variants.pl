#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);
sub fitsSSEProfile($$$$$$$$);

(@ARGV == 0) and usage("");

my ($vcfFile, $variantType, $minReadCount, $minAltCount, $minSNVReadRatio, $minIndelReadRatio, $minQuality, $minMapQ, $cnvMaxP, $maxCNVs, $filterVarTypes, $keepOnlyVarTypes, $filterSSE, $filterRandoms, $prevSeenThreshold, $maxMAF, $printFilteredFile, $doNotRemoveFiltered, $help, $verbose);
GetOptions(	'vcf=s' => \$vcfFile,
			'variantType=s', \$variantType,
			'minReadCount=i', \$minReadCount,
			'minAltCount=i', \$minAltCount,
			'minSNVReadRatio=f', \$minSNVReadRatio,
			'minIndelReadRatio=f', \$minIndelReadRatio,
			'minQ=i', \$minQuality,
			'minMapQ=i', \$minMapQ,
			'cnvMaxP=f', \$cnvMaxP,
			'maxCNVs=i', \$maxCNVs,
			'filterVarTypes=s', \$filterVarTypes,
			'keepOnlyVarTypes=s', \$keepOnlyVarTypes,
			'filterSSE', \$filterSSE,
			'filterRandoms', \$filterRandoms,
			'prevSeenThreshold=i', \$prevSeenThreshold,
			'maxMAF=f', \$maxMAF,
			'printFiltered=s' => \$printFilteredFile,
			'doNotRemoveFiltered' => \$doNotRemoveFiltered,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

$help and usage("");
$vcfFile or usage("Must specify input file name with argument --vcf");

if (defined $printFilteredFile)
{
	open(FILTERED, ">$printFilteredFile") || die "Failed to create file to print filtered variants: $printFilteredFile\n";
}

# Set up defaults
my $variantTypeStrict = 0;
(defined $variantType and $variantType eq "strict") and $variantTypeStrict = 1;

(defined $minReadCount) or $minReadCount = 3;
(defined $minAltCount) or $minAltCount = 3;
(defined $minSNVReadRatio) or $minSNVReadRatio = 0.2;
(defined $minIndelReadRatio) or $minIndelReadRatio = 0.15;
(defined $minQuality) or $minQuality = 20;
(defined $minMapQ) or $minMapQ = 15;
(defined $filterVarTypes) or $filterVarTypes = "";
(defined $prevSeenThreshold) or $prevSeenThreshold = 9999999;
(defined $maxMAF) or $maxMAF = 9;
(defined $filterSSE) or $filterSSE = 1;
my $numCNVs = 0;

open(VCF_FILE, "<$vcfFile") || die "Failed to open input file for reading: $vcfFile\n";

my %topCNVs;
if ($maxCNVs)
{
	my @cnvArray;
	# Go through the VCF file once first to make a list of all CNVs, so we can know
	# which the top N CNVs are.
	while(<VCF_FILE>)
	{
		(/^#/) and next; # skip header lines
		my @vcfLine = split(/\t/, $_, -1);
		my $chr = $vcfLine[0];
		my $pos = $vcfLine[1];
		my $qual = $vcfLine[5];
		my $vcfInfo = $vcfLine[7];
		if ($vcfInfo =~ /SVTYPE=CNV/)
		{
			push(@cnvArray, "$chr\t$pos\t$qual");
		}
	}
	seek(VCF_FILE, 0, 0);
	
	# Now output the full list of variants and IDs of samples in which each was seen.
	# Sort the variants by chromosomal position first.
	my %chrmOrder = (
		"1" => 1,
		"2" => 2,
		"3" => 3,
		"4" => 4,
		"5" => 5,
		"6" => 6,
		"7" => 7,
		"8" => 8,
		"9" => 9,
		"10" => 10,
		"11" => 11,
		"12" => 12,
		"13" => 13,
		"14" => 14,
		"15" => 15,
		"16" => 16,
		"17" => 17,
		"18" => 18,
		"19" => 19,
		"20" => 20,
		"21" => 21,
		"22" => 22,
		"X" => 23,
		"Y" => 24,
		"M" => 25
	);	

	# Sort the CNVs by quality
	my @sortedCNVArray = sort {
		my @a1 = split(/\t/, $a);
		my @b1 = split(/\t/, $b);
		return -($a1[2] <=> $b1[2]);
	} @cnvArray;
	for (my $i = 0; $i < $maxCNVs and $i <= $#sortedCNVArray; $i++)
	{
		my ($chr, $pos, $qual) = split(/\t/, $sortedCNVArray[$i]);
		$topCNVs{"$chr\t$pos"} = $qual;
	}
}

while(<VCF_FILE>)
{
	if (/^#/)
	{
		# This is a header line, so just print it.
		print;
		next;
	}
	chomp();
	my @vcfLine = split(/\t/, $_, -1);
	my $chr = $vcfLine[0];
	my $pos = $vcfLine[1];
	my $qual = $vcfLine[5];
	my $vcfInfo = $vcfLine[7];
	
	my ($variation, $altCount, $readCount, $dp4, $pv4, $mapQ, $numPrevSamples);
	my $oldFilter = $vcfLine[6];
	my $newFilter;
	my $thgMaf = 0;
	my $evsMaf = 0;
	my $cnvPval = 0;
	
	($vcfInfo =~ /VT=([^;\t]+)/)		and $variation = $1;
	($vcfInfo =~ /ALTC=([^;\t]+)/)		and $altCount = $1;
	($vcfInfo =~ /RDC=([^;\t]+)/)		and $readCount = $1;
	($vcfInfo =~ /DP4=([^;\t]+)/)		and $dp4 = $1;
	($vcfInfo =~ /PV4=([^;\t]+)/)		and $pv4 = $1;
	($vcfInfo =~ /MQ=([^;\t]+)/)		and $mapQ = $1;
	($vcfInfo =~ /PSN=([^;\t]+)/)		and $numPrevSamples = $1;
	($vcfInfo =~ /THGMAF=([^;\t]+)/)	and $thgMaf = $1;
	($vcfInfo =~ /EVSMAF=([^;\t]+)/)	and $evsMaf = $1;
	
	my $isIndel = ($vcfLine[3] eq "-" or $vcfLine[4] eq "-");
	my $isCNV = 0;
	if ($vcfInfo =~ /SVTYPE=CNV/)
	{
		$isCNV = 1;
		$numCNVs++;
		($vcfInfo =~ /P_BH=([^,;\t]+)/)	and $cnvPval = $1;
	}
	
	if ($filterVarTypes && $variation && $filterVarTypes =~ /^$variation$|^$variation,|,$variation$|,$variation,/)
	{
		$newFilter = "Variant type excluded";
	}
	elsif ($keepOnlyVarTypes && $variation && $keepOnlyVarTypes !~ /^$variation$|^$variation,|,$variation$|,$variation,/)
	{
		$newFilter = "Variant type excluded";
	}
	if ($isCNV)
	{
		if ($isCNV and defined $cnvMaxP and $cnvPval > $cnvMaxP)
		{
			$newFilter = "CNV P val > $cnvMaxP";
		}
		elsif ($isCNV and not exists $topCNVs{"$chr\t$pos"})
		{
			$newFilter = "CNV not in top $maxCNVs CNVs";
		}
	}
	else
	{
		if (defined $qual and $qual < $minQuality)
		{
			$newFilter = "Qual < $minQuality";
		}
		elsif ($readCount < $minReadCount)
		{
			$newFilter = "Read depth < $minReadCount";
		}
		elsif ($altCount < $minAltCount)
		{
			$newFilter = "Alt count < $minAltCount";
		}
		elsif ($mapQ < $minMapQ)
		{
			$newFilter = "MapQ < $minMapQ";
		}
		elsif ($numPrevSamples > $prevSeenThreshold)
		{
			$newFilter = "num prev seen samples > $prevSeenThreshold";
		}
		elsif ($thgMaf > $maxMAF or $evsMaf > $maxMAF)
		{
			$newFilter = "MAF > $maxMAF";
		}
		elsif ($variantTypeStrict and ($variation =~ /extended/))
		{
			$newFilter = "Extended splicing variant";
		}
		elsif (!$isIndel and ($altCount / $readCount) < $minSNVReadRatio)
		{
			$newFilter = "Alt read ratio < $minSNVReadRatio";
		}
		elsif ($isIndel and ($altCount / $readCount) < $minIndelReadRatio)
		{
			$newFilter = "Alt read ratio < $minIndelReadRatio";
		}
		elsif ($filterRandoms and $chr =~ /random|Un_/)
		{
			$newFilter = "random or unassembled chromosome";
		}
		if ($filterSSE and $dp4 and $pv4)
		{
			my ($fwdRef, $revRef, $fwdAlt, $revAlt) = split(/,/, $dp4);
			my ($strandBias, $baseQBias, $mapQBias, $endDistBias) = split(/,/, $pv4);
			if (defined $strandBias and defined $baseQBias and defined $mapQBias and defined $endDistBias
				and fitsSSEProfile($fwdRef, $revRef, $fwdAlt, $revAlt, $strandBias, $baseQBias, $mapQBias, $endDistBias))
			{
				$newFilter = "Probable SSE";
			}
		}
	}
	
	# Right now we will include the old filter in the Filter column if it was present,
	# but we won't REMOVE the variant based on this. Only remove the variant based on 
	# the current ("new") filter.
	my $filterToPrint = $newFilter;
	if ($newFilter)
	{
		if ($oldFilter ne ".")
		{
			$filterToPrint = $oldFilter.",".$newFilter;
		} else {
			$filterToPrint = $newFilter;
		}
	}
	else
	{
		$filterToPrint = $oldFilter;
	}
	
	if ($newFilter and $newFilter ne "." and $printFilteredFile)
	{
		print FILTERED join("\t", @vcfLine[0..5], $filterToPrint, @vcfLine[7..$#vcfLine])."\n";
	}
	
	if (!$newFilter or $newFilter eq "." or $doNotRemoveFiltered)
	{
		print join("\t", @vcfLine[0..5], $filterToPrint, @vcfLine[7..$#vcfLine])."\n";
	}
}
close(VCF_FILE);

if (defined $printFilteredFile)
{
	close(FILTERED);
}



###############################################################################

sub fitsSSEProfile($$$$$$$$)
{
	my ($fwdRef, $revRef, $fwdAlt, $revAlt, $strandBias, $baseQBias, $mapQBias, $endDistBias) = @_;
	return 0;
}

###############################################################################

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;
This program takes as parameters a sample name and a configuration file with information about
where to find sample sequencing reads. It runs a series of other tools to check read qualities,
align the sample to produce a BAM, and collect coverage information. Progress information is
written to STDOUT.

Usage:
$0 --sampleName NAME --sampleConfig FILE [options]

OPTIONS:
--vcf               FILE  File that is output by combine_annovar_files.pl
--minReadCount      int   Filters variants with read count < minReadCount
--minAltCount       int   Filters variants with alt count < minAltCount
--minSNVReadRatio   float Filters SNV variants with read ratio < minSNVReadRatio
--minIndelReadRatio float Filters INDEL variants with read count < minIndelReadRatio
--minQ              int   Filters variants with variant quality < minQ
--minMapQ           int   Filters variants with RMS mapping quality < minMapQ
--filterSSE               Filters variants with characteristics of sequence-specific error
--prevSeenThreshold int   Filters variants seen in > prevSeenThreshold total samples
--maxMAF            float Filters variants with MAF > maxMAF
--printFiltered     FILE  Removes filtered variants from output, but writes them to FILE
--doNotRemoveFiltered     Puts a note in the Filter column of the output but does not remove variants
--help|h                  Prints usage information

EOF
    exit 1;
}
