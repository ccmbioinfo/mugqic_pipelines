#!/usr/bin/env perl
# Sept 20, 2010
# Created: Feb. 14, 2011, KH
# Modified: 2012/3/7, JS - changed to work on VCF
#
# Query the OMIM database for disease-gene associations.
#
# See usage() subroutine at the bottom of this script for more info or
# call script with "-h" option.
use strict;
use warnings;
use Getopt::Long;

sub usage($);
sub trim($);

(@ARGV == 0) and usage("");

my ($help, $vcfFile, $omimdb, $verbose);
GetOptions(	'vcf=s' => \$vcfFile,
			'omimdb=s', \$omimdb,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

$help and usage("");
$vcfFile or usage("Must specify input file name with argument --vcf");
$omimdb or usage("Must specify omim db file with argument --db");


open OMIM_TSV, "<$omimdb" or die "Failed to open $omimdb: $!\n";

my %gene_desc;
my %gene_omim;

while (<OMIM_TSV>)
{
    chomp;
    my @line = split("\t");

    # store values as a "hash of a hash" to keep the annotation list unique
    $gene_desc{$line[0]}{$line[1]}++;
    $gene_omim{$line[0]}{$line[2]}++ if defined $line[2];
}
close OMIM_TSV;


###############################################################################

# keep track of number of gene associations found
my ($totalVariants, $omim, $desc) = (0,0,0);

# open input file
open VCF_FILE, "<$vcfFile" or die "Failed to open input file $vcfFile";
while (<VCF_FILE>)
{
	if (/^#/)
	{
		# This is a header line, so just print it.
		print;
		next;
	}
	chomp();
	my @vcfLine = split(/\t/, $_, -1);
	my $gene;
	($vcfLine[7] =~ /GENE=([^;\t]*)/) and $gene = trim($1);
	
	my ($geneName, $disorder);
	if (defined $gene)
	{
		my (%disorders, %gene_names);
        $gene =~ s/^"//;
        $gene =~ s/"$//;

        # find OMIM disease association
        if (defined $gene_omim{$gene})
		{
            my $dis = &trim( join("; ", keys %{$gene_omim{$gene}}) );
            $disorder = &verify_fields($dis);
            $omim++;
        }

        # get gene description
        if (defined $gene_desc{$gene})
		{
            my $name = &trim( join("; ", keys %{$gene_desc{$gene}}) );
            $geneName = &verify_fields($name);
            $desc++;
        }
	}
	
	my $geneAnnotation;
	my $disorderAnnotation;
	if (defined $geneName)
	{
		$vcfLine[7] .= ";GN=${geneName}";
	}
	if (defined $disorder)
	{
		$vcfLine[7] .= ";OMIM=${disorder}";
	}

    print join("\t", @vcfLine)."\n";

    $totalVariants++;
}
close(VCF_FILE);


# print summary
printf STDERR "# Found %d (%.1f%%) gene descriptions out of %d entries.\n",
    ($desc, $desc/$totalVariants*100, $totalVariants);
printf STDERR "# Found %d (%.1f%%) OMIM associations.\n",
    ($omim, $omim/$totalVariants*100);


###############################################################################

sub usage($) 
{
	print $_[0]."\n";
	print STDERR <<EOF;
Query the OMIM database for disease-gene associations and the DAVID database
for gene annotation.

The OMIM and DAVID database is stored as a tsv flat file that is specified with
the --omimdb option. (Stored at /data/GRID/exome_sequencing/omim/david_omim_genemap.tsv.)

The VCF file specified with --vcf is annotated by adding to the INFO column the
fields "GN=" (descriptive gene name) and "OMIM=" (omim disease association) if
such annotations are found. No fields are added if no annotation is found.

usage: $0 --vcf FILE --omimdb FILE > annotated.vcf

-h           : print usage info
-vcf FILE    : VCF file to annotate
-omimdb FILE : specifies the location of the omim TSV file

EOF

    exit 1;
} 

# trim whitespace
sub trim($)
{
    my $text = shift;

    $text =~ s/^\s+//;
    $text =~ s/\s+$//;

    return $text;
}

# check if the field contains delimiters, which may cause problems for
# future parsing
sub verify_fields()
{
    my $text = shift;
	# Replace any occurrences of ";" or "\t", which are delimiters in VCF
	$text =~ s/[;\t]/\/\//;
    return $text;
}
