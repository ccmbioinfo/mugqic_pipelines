#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);
sub trim($);

(@ARGV == 0) and usage("");

my ($mutationCountsFile, $vcfFile, $help, $verbose);
GetOptions(	'mutationCounts=s' => \$mutationCountsFile,
			'vcf=s', \$vcfFile,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

open(MUTATION_COUNTS_FILE, "<".$mutationCountsFile) or die("Failed to open region of mutation counts file: $mutationCountsFile\n");
open(VCF_FILE, "<".$vcfFile) or die("Failed to open variants input file: $vcfFile\n");


# Read in entire Mutation Counts file
my %geneDetailsHash;
while (<MUTATION_COUNTS_FILE>)
{
	chomp();
	(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
	(/^\s*$/) and next; #ignore white space lines
	my @line = split(/\t/);

	$geneDetailsHash{$line[0]} = $line[1].":".$line[2].":".$line[4];
}
close(MUTATION_COUNTS_FILE);


while (<VCF_FILE>)
{
	if (/^#/)
	{
		# This is a header line, so just print it.
		print;
		next;
	}
	chomp();
	my @line = split(/\t/, $_, -1);
	my $gene;
	($line[7] =~ /GENE=([^;\t]*)/) and $gene = trim($1);

	my $annotation = "-";
	if (defined $gene)
	{
		($verbose and !exists $geneDetailsHash{$gene}) and print STDERR "Note: Gene $gene not found when annotating gene mutation frequency.\n";
		(exists $geneDetailsHash{$gene}) and $annotation = $geneDetailsHash{$gene};
	}
	$line[7] .= ";GMF=".$annotation;
	print join("\t", @line)."\n";
}
close(VCF_FILE);


###############################################################################

# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub usage($)
{
	print $_[0]."\n";
    print STDERR <<EOF;
This program takes a vcf file as input and a file listing mutation counts for all genes, and
adds a field "GMF=" to the VCF INFO column with details of the gene mutation frequency in
the mutation counts file passed in.

Usage:
$0 --vcf FILE --mutationCounts FILE > output.vcf

OPTIONS:
--vcf             FILE  VCF file of variants to annotate
--mutationCounts  FILE  File listing the counts of rare mutations in each gene
--help|h                Prints usage information

EOF
    exit 1;
}
