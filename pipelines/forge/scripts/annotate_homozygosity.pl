#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub homHetOrder($);
sub usage($);

(@ARGV == 0) and usage("");

my ($vcfFile, $help, $verbose);
GetOptions(	'vcf=s' => \$vcfFile,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

open(VCF_FILE, "<".$vcfFile) or die("Failed to open input file: $vcfFile");
my @vcfLines = <VCF_FILE>;
close(VCF_FILE);

my %timesGeneSeen;
my @variantClassArray;

for (my $i=0; $i <= $#vcfLines; $i++)
{
	my $inputLine = $vcfLines[$i];
	if ($inputLine =~ /^#/)
	{
		# This is a header line, skip it for now. But add an empty variant "class".
		push(@variantClassArray, 0);
		next;
	}
	chomp($inputLine);

	my @vcfFields = split(/\t/, $inputLine);
	my $isIndel = ($vcfFields[3] eq "-" or $vcfFields[4] eq "-");

	my $altCount = 0;
	my $readCount = 1;
	($vcfFields[7] =~ /ALTC=([^;\t]+)/)		and $altCount = $1;
	($vcfFields[7] =~ /RDC=([^;\t]+)/)		and $readCount = $1;

	my $refBaseCount = $readCount - $altCount;
	my $pctAlt = $altCount / $readCount;
	my $class = "het";
	if (($readCount >= 8 && $altCount == $readCount) or
	    ($readCount >= 15 && $refBaseCount <= 1) or
		($isIndel and $readCount >= 20 && $pctAlt >= 0.90))
	{
		$class = "hom";
	}
	elsif (($altCount == $readCount) or
	       ($readCount >= 4 && $refBaseCount <= 1) or
	       ($readCount >= 15 && $refBaseCount <= 2) or
	       ($readCount >= 30 && $pctAlt >= 0.85) or
		   ($isIndel and $readCount >= 10 && $pctAlt >= 0.70))
	{
		# Possibly a homozygote
		$class = "possibly hom";
	}
	push(@variantClassArray, $class);

	my $gene;
	($vcfFields[7] =~ /GENE=([^;\t]*)/) and $gene = $1;
	if ($gene)
	{
		(! exists $timesGeneSeen{$gene}) and $timesGeneSeen{$gene} = 0;
		$timesGeneSeen{$gene}++;
	}
}

# Change multiple "het" calls in the same gene to "multiple het"
for (my $i=0; $i <= $#vcfLines; $i++)
{
	($vcfLines[$i] =~ /^#/) and next; # header line, skip for now
	my @vcfFields = split(/\t/, $vcfLines[$i]);
	my $gene;
	($vcfFields[7] =~ /GENE=([^;\t]*)/) and $gene = $1;
	if ($gene and $variantClassArray[$i] eq "het" and $timesGeneSeen{$gene} >= 2)
	{
		$variantClassArray[$i] = "multiple het";
	}

	chomp($vcfFields[7]);
	$vcfFields[7] .= ";HMZ=".$variantClassArray[$i]."\n";
	$vcfLines[$i] = join("\t", @vcfFields);
}


# my @sortedVariants = sort {
	# my @a1 = split("\t", $a);
	# my @b1 = split("\t", $b);
	# if (homHetOrder($a1[6]) < homHetOrder($b1[6])) {
		# return -1;
	# }
	# elsif (homHetOrder($b1[6]) < homHetOrder($a1[6])) {
		# return 1;
	# }
	# # sort by gene symbol
	# return ($a1[11] cmp $b1[11]);
# } @variantLines;


for (my $i=0; $i <= $#vcfLines; $i++)
{
	print $vcfLines[$i];
}


###############################################################################
sub homHetOrder($)
{
	($_[0] eq "hom") and return 0;
	($_[0] eq "possibly hom") and return 1;
	($_[0] eq "multiple het") and return 2;
	($_[0] eq "het") and return 3;
	return 4;
}

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;

Usage:
$0 --input variantfile.txt > variantfile.hom.txt

Adds annotation of homozygosity (hom, possibly hom, multiple het, het) as column 6,
based on read counts of each variant.

--input FILE    File output by combine_annovar_files.pl

EOF

    exit;
}
