#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);
sub getChrShort($);

my ($help, $verbose, $variantsDB, $vcfFile, $detailsThreshold);

GetOptions(	'variantsDB=s' => \$variantsDB,
			'vcfFile=s' => \$vcfFile,
			'detailsThreshold=i' => \$detailsThreshold,
			'v|verbose' => \$verbose,
			'h|help' => \$help) or usage("");

$help and usage("");
$variantsDB or usage("ERROR: Missing variantsDB argument.\n");
$vcfFile or usage("ERROR: Missing vcfFile argument.\n");
$detailsThreshold or $detailsThreshold = 50;

open(VARDB, "<$variantsDB") or die "Failed to open variants DB file: $!";
open(VCF_FILE, "<".$vcfFile) or die "Failed to open $vcfFile: $!";

# First read in the whole VCF file. Since the VCF is much smaller than the DB, it is
# better to store the VCF in memory and go through the DB a line at a time.
my $t = localtime();
$verbose and print STDERR "$t\tReading VCF file into memory...\n";
my %vcfVariants;
my @vcfVariantKeys; # Keep sorted copy of keys to avoid having to sort all keys later
while (<VCF_FILE>)
{
	(/^\s*#/) and print and next; #ignore comment lines (first non-ws character is #)
	chomp();
	my $numErrors = 0;

	my ($chrmStr, $pos, $id, $ref, $alt, $qual, $filter, $info, @restOfLine) = split(/\t/);
	my $chrm = getChrShort($chrmStr);
	if (!$chrm)
	{
		$numErrors++;
		if ($numErrors <= 10)
		{
			print STDERR "Warning: invalid chromosome '$chrmStr' in file $vcfFile line $.\n";
			($numErrors == 10) and print STDERR "10 chromosome errors printed... no more warnings will be shown for this file.\n\n";
		}
		next;
	}

	my @altAlleles = split(/,/, uc($alt));
	foreach my $altAllele (@altAlleles)
	{
		my $key = join("\t", $chrm, $pos, uc($ref), $altAllele);
		$vcfVariants{$key} = join("\t", $chrm, $pos, $id, $ref, $altAllele, $qual, $filter, $info.";PSN=0;PS=-", @restOfLine);
		push(@vcfVariantKeys, $key);
	}
}
close(VCF_FILE);
my $numVariants = scalar(keys %vcfVariants);
$t = localtime();
print STDERR "$t\tFinished reading VCF file: $numVariants variants.\n";


# Now read the variants DB.
# First read the table of sample info. All lines beginning with ## are sample info
# in the format:
# sampleName\tsampleID\tsampleIsCancer\tsamplePhenotype
#
# sampleName is a text string
# sampleID is an integer
# sampleIsCancer is 0 or 1 (false or true)
# samplePhenotype is a text string
my @sampleNames;
my @sampleIsCancer;
my @samplePhenotypes;
$verbose and print STDERR "Reading variants db: $variantsDB\n";
my $numSamplesRead = 0;
while (<VARDB>)
{
	if ($_ !~ /^##/)
	{
		last;
	}
	chomp();

	if ($_ =~ /^##ID_CUTOFF/)
	{
		my ($junk, $cutoff) = split(/=/);
		if ($cutoff)
		{
			if ($detailsThreshold > $cutoff)
			{
				print STDERR "Warning: detailsThreshold specified was $detailsThreshold, but this is higher than the variant db sample count cutoff of $cutoff. Details threshold will be reduced to $cutoff. You must rebuild the DB from scratch if you want a higher cutoff.\n";
				$detailsThreshold = $cutoff;
				next;
			}
		}
	}
	if ($_ =~ /^##SAMPLE=(<.*>)/)
	{
		my ($sampleID, $sampleName, $isCancer, $phenotype);
		my $str = $1;
		($str =~ /ID=([^\t]*)[\t>]/) and $sampleID = $1;
		($str =~ /Name=([^\t]*)[\t>]/) and $sampleName = $1;
		($str =~ /Cancer=([^\t]*)[\t>]/) and $isCancer = $1;
		($str =~ /Phenotype=["]?(.*)["]?[\t>]/) and $phenotype = $1;
		if (!$sampleID or !$sampleID or !$sampleID or !$sampleID)
		{
			die "ERROR: Sample line format not recognized: line $.:\n".$_."\n";
		}

		$verbose and print STDERR "Read from existingdb sample $sampleID, $sampleName, $isCancer, \"$phenotype\"\n";
		# sampleID should be an integer
		($sampleID =~ /^\d+$/) or die "ERROR: sample ID '$sampleID' is not an integer in database $variantsDB.";

		$sampleNames[$sampleID] = $sampleName;
		$sampleIsCancer[$sampleID] = $isCancer;
		$samplePhenotypes[$sampleID] = $phenotype;
		$numSamplesRead++;
	}
}
$t = localtime();
$verbose and print STDERR "$t\tRead sample table: $numSamplesRead samples\n";

my $numVariantsFound = 0;
while (<VARDB>)
{
	/^\s*#/ and next;
	chomp();
	my @line = split(/\t/);
	my $key = join("\t", @line[0 .. 3]);

	unless (exists $vcfVariants{$key}) { next;}
	$numVariantsFound++;

	my ($chrm, $pos, $id, $ref, $alt, $qual, $filter, $info, @restOfLine) = split(/\t/, $vcfVariants{$key});

	my $prevSampleNames = "-";
	if ($line[4] <= $detailsThreshold)
	{
		my @prevSamples = split(/,/, $line[5]);
		my $first = 1;
		for (my $i = 0; $i <= $#prevSamples; $i++)
		{
			if ($first)
			{
				$first = 0;
				$prevSampleNames = $sampleNames[$prevSamples[$i]];
			} else {
				$prevSampleNames .= ",".$sampleNames[$prevSamples[$i]];
			}
		}
	}
	my $addInfoString = "PSN=".$line[4].";PS=".$prevSampleNames;
	$info =~ s/PSN=0;PS=-/$addInfoString/;

	if ($line[6])
	{
		# The variant DB had a filter flag on the variant. Put this in the filter
		# field of the VCF line.
		if ($filter and $filter ne ".")
		{
			$filter .= ",".$line[6];
		} else {
			$filter = $line[6];
		}
	}

	$vcfVariants{$key} = join("\t", $chrm, $pos, $id, $ref, $alt, $qual, $filter, $info, @restOfLine);
}
close(VARDB);
$t = localtime();
$verbose and print STDERR "$t\tFinished going through variants DB: $numVariantsFound variants found.\n";


for (my $i = 0; $i <= $#vcfVariantKeys; $i++)
{
	print "chr".$vcfVariants{$vcfVariantKeys[$i]}."\n";
}
$t = localtime();
$verbose and print STDERR "$t\tFinished outputting variants.\n";


###############################################################################
sub getChrShort($)
{
	($_[0] =~ /^(chr)?(X|Y|M|[1-9][0-9]?)$/i) and return uc($2);
	return $_[0];
}

###############################################################################
sub usage($)
{
	($_[0]) and print $_[0]."\n";
    print STDERR <<EOF;
Usage:
annotate_vardb.pl --vcfFile sample.vcf --variantsDB FILE [options] > sample.annotated.vcf

Adds entries to the INFO field of each variant in the input VCF specified with --vcfFile.
The new field "PSN=" indicates the number of times the variant was seen before, while the
new field "PS=" indicates the names of samples in which it was seen.

--help                   print usage info
--vcfFile                VCF file to annotate
--variantsDB             DB of variants to use for annotation, built with build_variants_db.pl
--detailsThreshold INT   If variant seen more than this many times, sample IDs will be reported as '-'

EOF

    exit;
}

