#!/usr/bin/env perl
# This script creates a master list of previously seen variants (BOTH SNVs and indels). 
# As input it takes a list of VCF filenames and IDs to associate with each file. The
# ID for each file must be unique.
# Optionally, it can also apply a filter to the entire dataset to determine likely 
# false variants. This option takes more time and memory. If specified, an additional
# column is added with the filter that failed.
use strict;
use warnings;
use Getopt::Long;

sub usage($);
sub trim($);
sub getChrShort($);
sub getSampleNames($);
sub doCalculateFilter($);
sub approxBinomialPValue($$$);

my ($vcflist, $existingdb, $calculateFilter, $samplecountcutoff, $minAltReads, $sortOutput, $help, $verbose);
GetOptions(	'vcflist=s' => \$vcflist,
			'existingdb=s' => \$existingdb,
			'calcFilter' => \$calculateFilter,
			'samplecountcutoff=i' => \$samplecountcutoff,
			'minAltReads=i' => \$minAltReads,
			'sort' => \$sortOutput,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

$help and usage("");
$vcflist or usage("ERROR: Missing vcflist argument.\n");

($existingdb and $calculateFilter) and usage("ERROR: Cannot specify --calcFilter and --existingdb. The filter can only be calculated when the DB is built from scratch.\n");

my $dp4StrandBiasPhredFilter = 5;
my $dp4StrandBiasEffectSizeFilter = 0.2;

my $pv4StrandBiasFilter = 30;
my $pv4StrandBiasEffectSizeFilter = 0.5;

my $pv4BaseQBiasFilter = 50;
my $pv4BaseQBiasEffectSizeFilter = 5.0;

my $pv4MapQBiasFilter = 30;
my $pv4MapQBiasEffectSizeFilter = 4.0;

my $pv4EndDistanceBiasFilter = 50;
my $pv4EndDistanceBiasEffectSizeFilter = 2.0;

($minAltReads) or $minAltReads = 0;

# If an existing database was specified with --existingdb, then we first load
# that DB's variants. This is important since the IDs for samples in the existing DB
# will go from 1...numSamples. We want to keep the same ordering, and new VCF files
# in the vcflist will get subsequent sample ID numbers.
my $lastSampleID = 0;
my %sampleIDs;
my @sampleNames;
my %sampleIsCancer;
my %samplePhenotypes;
if (defined $existingdb)
{
	print STDERR "Reading table of samples from existing variants DB: $existingdb.\n";
	open(EXISTING_DB, "<", $existingdb) or die "Failed to open database file: $existingdb\n";

	my $existingDBSampleCountCutoff = 0;
	my ($existingDBNumSNPs, $existingDBNumIndels) = (0, 0);
	while (<EXISTING_DB>)
	{
		if ($_ !~ /^##/)
		{
			last;
		}
		chomp();
		if ($_ =~ /^##ID_CUTOFF/)
		{
			my ($junk, $existingDBSampleCountCutoff) = split(/=/);
			if ($existingDBSampleCountCutoff)
			{
				if ($samplecountcutoff and $samplecountcutoff > $existingDBSampleCountCutoff)
				{
					die "ERROR: samplecountcutoff specified was $samplecountcutoff, but this is higher than the existing db sample count cutoff of $existingDBSampleCountCutoff. You must rebuild the DB from scratch if you want a higher cutoff.\n";
				}
				$samplecountcutoff = $existingDBSampleCountCutoff;
			}
			next;
		}
		if ($_ =~ /^##SAMPLE=(<.*>)/)
		{
			my ($sampleID, $sampleName, $isCancer, $phenotype);
			my $str = $1;
			($str =~ /ID=([^\t]*)[\t>]/) and $sampleID = $1;
			($str =~ /Name=([^\t]*)[\t>]/) and $sampleName = $1;
			($str =~ /Cancer=([^\t]*)[\t>]/) and $isCancer = $1;
			($str =~ /Phenotype=([^\t]*)[\t>]/) and $phenotype = $1;
			if (!$sampleID or !$sampleID or !$sampleID or !$sampleID)
			{
				die "ERROR: Sample line format not recognized: line $.:\n".$_."\n";
			}
			$sampleName = trim($sampleName);
			
			$verbose and print STDERR "Read from existingdb sample $sampleID, $sampleName, $isCancer, \"$phenotype\"\n";
			# sampleID should be an integer
			($sampleID =~ /^\d+$/) or die "ERROR: sample ID '$sampleID' is not an integer in database $existingdb.";
			
			$sampleIDs{$sampleName} = $sampleID;
			$sampleNames[$sampleID] = $sampleName;
			$sampleIsCancer{$sampleName} = $isCancer;
			$samplePhenotypes{$sampleName} = $phenotype;
			$lastSampleID = $sampleID;
			next;
		}
	}
}
$samplecountcutoff or $samplecountcutoff = 1000;


open(VCF_LIST, "<", $vcflist) or die "Failed to open file listing samples to put in DB: $vcflist\n";
print STDERR "Reading list of new VCF files from file: $vcflist.\n";

# Read in entire list of variant filenames, and map them to their IDs
my @varFilesArray;
while (<VCF_LIST>)
{
	chomp();
	(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
	(/^\s*$/) and next; #ignore white space lines

	my ($includeSample, $sampleName, $file, $isCancer, $phenotype) = split(/\t/);
	
	(! $includeSample) and next;
	$sampleName = trim($sampleName);

	if (!$file)
	{
		$file = "${sampleName}.flt.vcf";
	}
	unless (-e $file)
	{
		print STDERR "File $file DOESN'T EXIST. Variants from this sample will not be in the master list.\n";
		next;
	}
	unless ($sampleName)
	{
		print STDERR "No sample name specified in column 2 for file \"$file\". Variants from this sample will not be in the master list.\n";
		next;
	}
	unless (defined $isCancer)
	{
		die "No true/false value specified for column 4, 'is cancer', for file \"$file\".\n";
	}
	if ($isCancer =~ /true/i or $isCancer eq "1")
	{
		$isCancer = 1;
	}
	elsif ($isCancer =~ /false/i or $isCancer eq "0")
	{
		$isCancer = 0;
	}
	else
	{
		die "Invalid true/false value specified for column 4, 'is cancer', for file \"$file\".\n";
	}
	
	if (exists $sampleIDs{$sampleName})
	{
		(!defined $existingdb) and die "ERROR: Sample name '$sampleName' was present more than once in the VCF list.\n";
		$verbose and print STDERR "Sample $sampleName was present in the existing database (or was present more than once in the VCF list).\n";
	}
	else
	{
		$verbose and print STDERR "Adding file: $file\n";
		push(@varFilesArray, join("\t", $includeSample, $sampleName, $file, $isCancer, $phenotype));
		# Assign the next integer as the ID for this sample
		$lastSampleID++;
		$sampleIDs{$sampleName} = $lastSampleID;
		$sampleNames[$lastSampleID] = $sampleName;
		$sampleIsCancer{$sampleName} = $isCancer;
		$samplePhenotypes{$sampleName} = $phenotype;
	}
}
close(VCF_LIST);


#Read in all variants in the existing variants DB
my %variantSamples; # For each variant, contains first the number of samples it was seen in, followed by the sample IDs
if (defined $existingdb)
{
	my $t = localtime();
	print STDERR "\n$t\tReading all variants from existing variants DB: $existingdb.\n";
	my ($existingDBNumSNPs, $existingDBNumIndels) = (0,0);
	while (<EXISTING_DB>)
	{
		chomp();
		(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
		
		my ($chrm, $pos, $ref, $alt, $numSamples, $sampleIDs) = split(/\t/);
		
		my $key = join("\t", $chrm, $pos, $ref, $alt);
		(exists $variantSamples{$key}) and die "Error in existing DB: identical variant key read more than once. Key=$key\n";
		if ($numSamples > $samplecountcutoff)
		{
			$sampleIDs = "-";
		}
		$variantSamples{$key} = "$numSamples\t$sampleIDs";
		
		my $isSNV = (length($ref) == 1 and length($alt) == 1);
		if (!$isSNV) #indels - right now they are handled the same as SNVs
		{
			$existingDBNumIndels++;
		} else {
			$existingDBNumSNPs++;
		}
	}
	close(EXISTING_DB);
	print STDERR "Read variants from existing DB: #indels=$existingDBNumIndels; #snps=$existingDBNumSNPs\n";
}


my $numSamplesToAdd = scalar(@varFilesArray);
my $t = localtime();
print STDERR "\n$t\tReady to add $numSamplesToAdd vcf files to create a database.\n\n";

my %pv4Values;
my %dp4Values;

for (my $fileIndex = 0; $fileIndex <= $#varFilesArray; $fileIndex++)
{
	my $varFile = $varFilesArray[$fileIndex];

	my ($curFileNumSNPs, $curFileNumIndels) = (0, 0);
	my ($includeSample, $curSampleName, $filename, $isCancer, $phenotype) = split(/\t/, $varFile);
	my $curSampleID = $sampleIDs{$curSampleName};
	($curSampleID) or die "Sample ID not assigned for sample $curSampleName, file $filename.\n";
	my $numErrors = 0;

	open(VARIANT_FILE, "<", $filename) or die("Failed to open variant file $filename.\n");
	while (<VARIANT_FILE>)
	{
		chomp();
		(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
		
		my ($chrmStr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split(/\t/);
		my $chrm = getChrShort($chrmStr);
		if (!$chrm)
		{
			$numErrors++;
			if ($numErrors <= 10)
			{
				print STDERR "Warning: invalid chromosome '$chrmStr' in file $filename line $.\n";
				($numErrors == 10) and print STDERR "10 chromosome errors printed... no more warnings will be shown for this file.\n\n";
			}
			next;
		}
		
		my ($dp4, $pv4);
		($info =~ /DP4=([^;\t]*)/) and $dp4 = $1;
		if ($dp4 and $minAltReads)
		{
			my @dp4Vals = split(/[;,]/, $dp4, -1);
			my $altReads = $dp4Vals[2] + $dp4Vals[3];
			if ($altReads < $minAltReads)
			{
				next;
			}
		}
		
		# Determine the key for this variant -- note that there can be multiple ALT alleles
		# specified. We add each of these separately. Also note that the same variant could be
		# described with different REF and ALT fields; however, we just keep the position and
		# whatever REF/ALT was specified in the VCF -- we make no adjustments to the position
		# of indels based on different REF/ALT alleles.
		my @altAlleles = split(/,/, $alt);
		foreach my $altAllele (@altAlleles)
		{
			my $isSNV = (length($ref) == 1 and length($altAllele) == 1);
			my $key = $chrm."\t".$pos."\t".uc($ref)."\t".uc($altAllele);

			if (!$isSNV) #indels - right now they are handled the same as SNVs
			{
				if (exists $variantSamples{$key})
				{
					my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
					++$numSamples;
					if ($numSamples > $samplecountcutoff)
					{
						$variantSamples{$key} = "$numSamples\t-";
					} else {
						$variantSamples{$key} = "$numSamples\t$sampleIDs,$curSampleID";
					}
				} else {
					$variantSamples{$key} = "1\t".$curSampleID;
				}
				$curFileNumIndels++;
			} else {
				if (exists $variantSamples{$key})
				{
					my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
					++$numSamples;
					if ($numSamples > $samplecountcutoff)
					{
						$variantSamples{$key} = "$numSamples\t-";
					} else {
						$variantSamples{$key} = "$numSamples\t$sampleIDs,$curSampleID";
					}
				} else {
					$variantSamples{$key} = "1\t".$curSampleID;
				}
				$curFileNumSNPs++;
			}
			
			if ($calculateFilter)
			{
				# If DP4 or PV4 fields are specified then add them to a list of these fields for
				# each sample this variant was seen in. Later we will use these to determine
				# filter status for the variant.
				($info =~ /DP4=([^;\t]*)/) and $dp4 = $1;
				if ($dp4) # Only try to include this sample variant if DP4 is specified.
				{
					($info =~ /PV4=([^;\t]*)/) and $pv4 = $1;
					(!$pv4) and $pv4 = ",,,";
					if (exists $pv4Values{$key})
					{
						$pv4Values{$key} .= ";$pv4";
						$dp4Values{$key} .= ";$dp4";
					} else {
						$pv4Values{$key} = "$pv4";
						$dp4Values{$key} = "$dp4";
					}
				}
			}
		}
	}
	
	my $fileNum = $fileIndex+1;
	my $numFiles = scalar(@varFilesArray);
	print STDERR "File $fileNum/$numFiles\t$filename: #indels=$curFileNumIndels; #snps=$curFileNumSNPs\n";
	close(VARIANT_FILE);
}

# First output a table of sample information -- sample IDs and names.
my @sortedNames = sort {
	return $sampleIDs{$a} <=> $sampleIDs{$b};
} keys %sampleIDs;

print "##fileformat=VARDBv0.1\n";
if ($samplecountcutoff)
{
	print "##ID_CUTOFF=".$samplecountcutoff."\n";
}
foreach my $sampleName (@sortedNames)
{
	print "##SAMPLE=<ID=$sampleIDs{$sampleName}\tName=$sampleName\tCancer=$sampleIsCancer{$sampleName}\tPhenotype=$samplePhenotypes{$sampleName}>\n";
}

# print "\n";

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

print STDERR "\n";
my @sortedKeys;
if ($sortOutput)
{
	$t = localtime();
	print STDERR "$t\tSorting variant keys...\n";
	# Sort the VCF file variant keys by chrm and position
	@sortedKeys = sort {
		my @a1 = split(/\t/, $a);
		my @b1 = split(/\t/, $b);
		
		# Get ordering, if chrm is 1-22 or X, Y, M
		my $chrAorder = $chrmOrder{uc($a1[0])};
		my $chrBorder = $chrmOrder{uc($b1[0])};
		if (!defined $chrAorder)
		{
			if (!defined $chrBorder)
			{
				# Both chrs are unknown
				my $val = (uc($a1[0]) cmp uc($b1[0]));
				if ($val == 0)
				{
					$val = ($a1[1] <=> $b1[1]);
				}
				return $val;
			}
			return 1; # a is not defined but b is, so a > b
		}
		if (!defined $chrBorder)
		{
			return -1; # a is defined by b is not, so a < b
		}
		if ($chrAorder == $chrBorder)
		{
			return $a1[1] <=> $b1[1];
		}
		return $chrAorder <=> $chrBorder;
	} keys %variantSamples;
}
else
{
	@sortedKeys = keys %variantSamples;
}

$t = localtime();
print STDERR "$t\tOutputting variants\n";

my $numVariantsSBFilter = 0;
my $numVariantsMQFilter = 0;
my $numVariantsBQFilter = 0;
my $numVariantsEDFilter = 0;

my $totalVariants = 0;
foreach my $key (@sortedKeys)
{
	my $filter = "";
	if ($calculateFilter)
	{
		$filter = "\t".doCalculateFilter($key);
	}
	print $key."\t".$variantSamples{$key}.$filter."\n";
	my ($numSamples) = split(/\t/, $variantSamples{$key});
	$totalVariants += $numSamples;
}

my $uniqueVariants = scalar keys %variantSamples;
$t = localtime();
print STDERR "$t\nTotal variants:  ".$totalVariants."\tUnique variants: ".$uniqueVariants."\n";

if ($calculateFilter)
{
	print STDERR "Variants failing strand bias filter: $numVariantsSBFilter\n";
	print STDERR "Variants failing baseQ bias filter:  $numVariantsBQFilter\n";
	print STDERR "Variants failing mapQ bias filter:   $numVariantsMQFilter\n";
	print STDERR "Variants failing end distance bias filter: $numVariantsEDFilter\n";
}

###############################################################################
sub getChrShort($)
{
	($_[0] =~ /^(chr)?(X|Y|M|[1-9][0-9]?)/i) and return uc($2);
	return $_[0];
}

sub getSampleNames($)
{
	my @sampleIDs = split(/,/, $_[0]);
	my $namesStr = "";
	foreach my $id (@sampleIDs)
	{
		(length($namesStr) > 0) and $namesStr .= ",";
		$namesStr .= $sampleNames[$id];
	}
	return $namesStr;
}

sub doCalculateFilter($)
{
	my $key = $_[0];
	my $filter = "";
	# Split the PV4 and DP4 items on either ; or , so that each consecutive 4 array items
	# corresponds to one sample. This avoids having to call split once for each sample.
	if (exists $dp4Values{$key})
	{
		my ($chr, $pos, $ref, $alt) = split('\t', $key);
		my $isSNV = (length($ref) == 1 and length($alt) == 1);
		
		my @dp4Vals = split(/[;,]/, $dp4Values{$key}, -1);
		my @pv4Vals = split(/[;,]/, $pv4Values{$key}, -1);
		# $verbose and print STDERR "key($key)   DP4=$dp4Values{$key}    VALS=(@dp4Vals)\n";
		# $verbose and print STDERR "key($key)   PV4=$pv4Values{$key}    VALS=(@pv4Vals)\n";
		(scalar(@pv4Vals) == scalar(@dp4Vals)) or die("ERROR: Assertion failed: Calculating filter - DP4 and PV4 arrays should have the same size. $key; DP4=(@dp4Vals); PV4=(@pv4Vals)\n");
		# Put each PV4 p-value type into a separate array (strand bias, baseQ bias, mapQ bias, end dist bias)
		# my (@pv4SB, @pv4BaseQB, @pv4MapQB, @pv4EndDB);
		my ($dp4RefFwd, $dp4RefRev, $dp4AltFwd, $dp4AltRev) = (0, 0, 0, 0);
		my ($sbPV4Product, $bqPV4Product, $mqPV4Product, $edPV4Product) = (1, 1, 1, 1);
		for (my $i=0; $i <= $#pv4Vals; $i+=4)
		{
			my $readCount = $dp4Vals[$i] + $dp4Vals[$i+1] + $dp4Vals[$i+2] + $dp4Vals[$i+3];
			$dp4RefFwd += $dp4Vals[$i];
			$dp4RefRev += $dp4Vals[$i+1];
			$dp4AltFwd += $dp4Vals[$i+2];
			$dp4AltRev += $dp4Vals[$i+3];
			if (length($pv4Vals[$i]) > 0)
			{
				# Transform the P-values to get natural log(p) over readcount
				# According to presentations on genome.gov this gives a better estimate of
				# likelihood the variant is real
				# ($pv4Vals[$i] > 0)   and push(@pv4SB, -log($pv4Vals[$i]) / $readCount) and $sbPV4Product += log($pv4Vals[$i]);
				# ($pv4Vals[$i+1] > 0) and push(@pv4BaseQB, -log($pv4Vals[$i+1]) / $readCount) and $bqPV4Product += log($pv4Vals[$i+1]);
				# ($pv4Vals[$i+2] > 0) and push(@pv4MapQB,  -log($pv4Vals[$i+2])  / $readCount) and $mqPV4Product += log($pv4Vals[$i+2]);
				# ($pv4Vals[$i+3] > 0) and push(@pv4EndDB,  -log($pv4Vals[$i+3]) / $readCount) and $edPV4Product += log($pv4Vals[$i+3]);
				($pv4Vals[$i] > 0)   and $sbPV4Product += log($pv4Vals[$i]);
				($pv4Vals[$i+1] > 0) and $bqPV4Product += log($pv4Vals[$i+1]);
				($pv4Vals[$i+2] > 0) and $mqPV4Product += log($pv4Vals[$i+2]);
				($pv4Vals[$i+3] > 0) and $edPV4Product += log($pv4Vals[$i+3]);
			}
		}
		
		my $totalReadCount = $dp4RefFwd + $dp4RefRev + $dp4AltFwd + $dp4AltRev;
		my $refCount = $dp4RefFwd + $dp4RefRev;
		my $altCount = $dp4AltFwd + $dp4AltRev;
		my $minRefAltCount = ($refCount < $altCount) ? $refCount : $altCount;
		my ($sbPV4stat, $bqPV4stat, $mqPV4stat, $edPV4stat) = (0,0,0,0);
		my ($sbPV4effectSize, $bqPV4effectSize, $mqPV4effectSize, $edPV4effectSize) = (0,0,0,0);
		if ($minRefAltCount > 5)
		{
			$sbPV4stat = -10*$sbPV4Product / log(10);
			$bqPV4stat = -10*$bqPV4Product / log(10);
			$mqPV4stat = -10*$mqPV4Product / log(10);
			$edPV4stat = -10*$edPV4Product / log(10);
			$sbPV4effectSize = $sbPV4stat / $minRefAltCount;
			$bqPV4effectSize = $bqPV4stat / $minRefAltCount;
			$mqPV4effectSize = $mqPV4stat / $minRefAltCount;
			$edPV4effectSize = $edPV4stat / $minRefAltCount;
		}
		
		if ($sbPV4stat > $pv4StrandBiasFilter and $sbPV4effectSize > $pv4StrandBiasEffectSizeFilter)
		{
			$numVariantsSBFilter++;
			(length($filter) > 1) and $filter .= ",";
			$filter .= "SB";
			my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
			my $sampleNames = "-";
			($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
			my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
			if ($verbose) { print STDERR "SB PV4\t$key2\t$sbPV4stat\t$sbPV4effectSize\t($dp4RefFwd, $dp4RefRev, $dp4AltFwd, $dp4AltRev)\t$numSamples\t$sampleNames\n"; }
		}
		if ($isSNV and $bqPV4stat > $pv4BaseQBiasFilter and $bqPV4effectSize > $pv4BaseQBiasEffectSizeFilter)
		{
			$numVariantsBQFilter++;
			(length($filter) > 1) and $filter .= ",";
			$filter .= "BQB";
			my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
			my $sampleNames = "-";
			($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
			my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
			if ($verbose) { print STDERR "BQ PV4\t$key2\t$bqPV4stat\t$bqPV4effectSize\t($dp4RefFwd, $dp4RefRev, $dp4AltFwd, $dp4AltRev)\t$numSamples\t$sampleNames\n"; }
		}
		if ($mqPV4stat > $pv4MapQBiasFilter and $mqPV4effectSize > $pv4MapQBiasEffectSizeFilter)
		{
			$numVariantsMQFilter++;
			(length($filter) > 1) and $filter .= ",";
			$filter .= "MQB";
			my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
			my $sampleNames = "-";
			($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
			my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
			if ($verbose) { print STDERR "MQ PV4\t$key2\t$mqPV4stat\t$mqPV4effectSize\t($dp4RefFwd, $dp4RefRev, $dp4AltFwd, $dp4AltRev)\t$numSamples\t$sampleNames\n"; }
		}
		if ($edPV4stat > $pv4EndDistanceBiasFilter and $edPV4effectSize > $pv4EndDistanceBiasEffectSizeFilter)
		{
			$numVariantsEDFilter++;
			(length($filter) > 1) and $filter .= ",";
			$filter .= "EDB";
			my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
			my $sampleNames = "-";
			($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
			my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
			if ($verbose) { print STDERR "ED PV4\t$key2\t$edPV4stat\t$edPV4effectSize\t($dp4RefFwd, $dp4RefRev, $dp4AltFwd, $dp4AltRev)\t$numSamples\t$sampleNames\n"; }
		}
		# else
		# {
			# Do a strand bias filter on all reads combined.
			# my $fwdCount = $dp4RefFwd + $dp4AltFwd;
			# my $revCount = $dp4RefRev + $dp4AltRev;
			# if ($fwdCount > 5 and $revCount > 5 and $refCount > 5 and $altCount > 5)
			# {
				# my $nullHypBinomialP = $fwdCount / $totalReadCount;
				# my $altBinomialP = $dp4AltFwd / $altCount;
				# my $refBinomialP = $dp4RefFwd / $refCount;
				# # Restrict possible binomial P values so that we don't get any zero P values
				# ($nullHypBinomialP > 0.98) and $nullHypBinomialP = 0.98;
				# ($altBinomialP > 0.98) and $altBinomialP = 0.98;
				# ($refBinomialP > 0.98) and $refBinomialP = 0.98;
				# ($nullHypBinomialP < 0.02) and $nullHypBinomialP = 0.02;
				# ($altBinomialP < 0.02) and $altBinomialP = 0.02;
				# ($refBinomialP < 0.02) and $refBinomialP = 0.02;
				# # print STDERR "nullp=$nullHypBinomialP;altp=$altBinomialP;refp=$refBinomialP\n";
				# my $logNullHypLikelihood = $fwdCount*log($nullHypBinomialP) + $revCount*log(1-$nullHypBinomialP);
				# my $logAltHypLikelihood = $dp4RefFwd*log($refBinomialP) + $dp4RefRev*log(1-$refBinomialP) + $dp4AltFwd*log($altBinomialP) + $dp4AltRev*log(1-$altBinomialP);
				# my $logLikelihoodRatio = $logNullHypLikelihood - $logAltHypLikelihood;
				# # Convert to a phred score
				# my $phredScore = -10*$logLikelihoodRatio / log(10);
				# # Determine the "effect size" - strand bias PER read
				# my $effectSize = $phredScore / $minRefAltCount;
				# if ($phredScore > $dp4StrandBiasPhredFilter) #&& $effectSize > $dp4StrandBiasEffectSizeFilter)
				# {
					# $filter .= "SB";
					# my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
					# my $sampleNames = "-";
					# my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
					# ($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
					# if ($verbose) { print STDERR "SB DP4\t$key2\t$phredScore\t$effectSize\t($dp4RefFwd, $dp4RefRev, $dp4AltFwd, $dp4AltRev)\t$numSamples\t$sampleNames\n"; }
				# }
			# }
		# }
		
		# if (scalar(@pv4SB) > 0 and &median(\@pv4SB) > $thresholdStrandBiasFilter)
		# {
			# $filter .= "SB";
			# my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
			# my $sampleNames = "-";
			# ($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
			# my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
			# if ($verbose) { my $med = &median(\@pv4SB); print STDERR "SB PV4\t$key2\t$med\t$numSamples\t$sampleNames\n"; }
		# }
		# elsif ($isSNV and scalar(@pv4BaseQB) > 0 and &median(\@pv4BaseQB) > $thresholdBaseQBiasFilter)
		# {
			# $filter .= "BQB";
			# my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
			# my $sampleNames = "-";
			# ($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
			# my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
			# if ($verbose) { my $med = &median(\@pv4BaseQB); print STDERR "BQ PV4\t$key2\t$med\t$numSamples\t$sampleNames\n"; }
		# }
		# elsif (scalar(@pv4MapQB) > 0 and &median(\@pv4MapQB) > $thresholdMapQBiasFilter)
		# {
			# $filter .= "MQB";
			# my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
			# my $sampleNames = "-";
			# ($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
			# my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
			# if ($verbose) { my $med = &median(\@pv4MapQB); print STDERR "MQ PV4\t$key2\t$med\t$numSamples\t$sampleNames\n"; }
		# }
		# elsif (scalar(@pv4EndDB) > 0 and &median(\@pv4EndDB) > $thresholdEndDistanceBiasFilter)
		# {
			# $filter .= "EDB";
			# my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
			# my $sampleNames = "-";
			# ($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
			# my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
			# if ($verbose) { my $med = &median(\@pv4EndDB); print STDERR "ED PV4\t$key2\t$med\t$numSamples\t$sampleNames\n"; }
		# }
		# elsif (scalar(@pv4SB) < 5 and (scalar(@dp4Vals)/4) > 5)
		# {
			# # Fewer than 5 samples had a PV4 value for the variant. Try to determine
			# # a filter based on strand bias across all variant calls.

			# # Estimate a binomial probability based on the ref base fwd/rev read counts.
			# # If there are too few ref reads, there is no point in this.
			# if (($dp4RefFwd + $dp4RefRev) > 5)
			# {
				# my $binomialP = $dp4RefFwd / ($dp4RefFwd + $dp4RefRev);
				# if (approxBinomialPValue($binomialP, ($dp4AltFwd + $dp4AltRev), $dp4AltFwd) < 0.001 or
					# approxBinomialPValue(1-$binomialP, ($dp4AltFwd + $dp4AltRev), $dp4AltRev) < 0.001)
				# {
					# $filter .= "SB";
					# if ($verbose)
					# {
						# my $binP1 = approxBinomialPValue($binomialP, ($dp4AltFwd + $dp4AltRev), $dp4AltFwd);
						# my $binP2 = approxBinomialPValue(1-$binomialP, ($dp4AltFwd + $dp4AltRev), $dp4AltRev);
						# my $key2 = $key; $key2 =~ s/\t/:/; ($key2) = split(/\t/, $key2);
						# my $sampleNames = "-";
						# my ($numSamples, $sampleIDs) = split(/\t/, $variantSamples{$key});
						# ($sampleIDs ne "-") and $sampleNames = getSampleNames($sampleIDs);
						# print STDERR "SB DP4\t$key2\tp1=$binP1,p2=$binP2,($dp4RefFwd, $dp4RefRev, $dp4AltFwd, $dp4AltRev)\t$numSamples\t$sampleNames\n";
					# }
				# }
			# }
		# }
		
	}
	return $filter;
}

sub median()
{
	my $refArray = $_[0];
	my @sortedArray = sort(@{$refArray});
	my $count = scalar(@sortedArray);
	if ($count % 2)
	{ 
		return $sortedArray[int($count/2)]; 
	}
	return ($sortedArray[$count/2] + $sortedArray[$count/2 - 1]) / 2; 
}

sub approxBinomialPValue($$$)
{
	my ($p, $n, $k) = @_;
	($k > $n*$p) and return 1; #Definitely in the expected range

	# Return easy to calculate upper bound for p value of observed k
	return exp(-2 * ($n*$p - $k)**2 / $n) / 2;
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

###############################################################################

sub usage($) 
{
	print STDERR $_[0]."\n";
    print STDERR <<EOF;
Usage:
build_variants_db.pl --vcflist FILE [options] > allVariantsDB.txt

Creates a text database of previously seen variants by combining the variant calls from
a set of VCF files. This database can then be used to annotate VCF files with new fields
PSN and PS that indicate the number of times a variant was seen, and in which samples.

The vcflist text file must have the following format (TAB-separated):
#Include in DB? Sample ID           VCF file       Is Cancer?    Phenotype
[0|1]           sampleID_for_file1  path_to_file1  [TRUE|FALSE]  text string for sample phenotype
[0|1]           sampleID_for_file1  path_to_file1  [TRUE|FALSE]  text string for sample phenotype

Additional options that can be passed:
--existingdb FILE         Specifies a DB previously built with build_variants_db.pl.
                          Variants in new samples will be added to this DB.
--samplecountcutoff INT   If a variant is seen in more than this many samples, sample names
                          will not be stored (saves memory)
--sort                    Variants in the output DB will be coordinate sorted. This is not really
                          important and takes at least 20 minutes extra for a large DB.

EOF

    exit;
}
