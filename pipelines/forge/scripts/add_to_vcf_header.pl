#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

sub usage($);

(@ARGV == 0) and usage("");

my ($inputFile, @linesToAdd, $help, $verbose);
GetOptions(	'input=s' => \$inputFile,
			'add=s' => \@linesToAdd,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

$help and usage("");
$inputFile or $inputFile="-";

open(INPUT_FILE, "<".$inputFile) or die("Failed to open input file: $inputFile");

my $linesAlreadyAdded = 0;

while (<INPUT_FILE>)
{
	(/^##/) and print and next; #skip past all header lines
	
	if (!$linesAlreadyAdded)
	{
		for (my $i=0; $i <= $#linesToAdd; $i++)
		{
			print $linesToAdd[$i]."\n";
		}
		$linesAlreadyAdded = 1;
	}
	print;
}
close(INPUT_FILE);


###############################################################################

sub usage($) 
{
	($_[0]) and print $_[0]."\n";
    print STDERR <<EOF;

Usage:
$0 --input FILE.vcf --add "lineToAdd" --add "secondLineToAdd" > FILE2.vcf

EOF

    exit;
}

