#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;
use Sys::Hostname;
use Cwd 'abs_path';

sub runCommand($);
sub dieScript($);


# Set so that files created by the pipeline will be writable by group members.
# Otherwise it becomes impossible for more than one person to work on a sample,
# since all files will be writable only by them.
umask 0002;

# Flush any output as soon as it is written
$| = 1;

## OPTIONS
my ($inputFile, $outputFile, $pgFile, $headerFile, $ref_genome, $bwa_ver);

GetOptions(
    'inputFile=s'              => \$inputFile,
    'outputFile=s'             => \$outputFile,
    'pgFile=s'                 => \$pgFile,
    'headerFile=s'             => \$headerFile,
    'ref_genome=s'             => \$ref_genome,
    'bwa_ver=s'                => \$bwa_ver
);

open(PG_FILE, "<", $pgFile) or dieScript("Failed to open file for recording programs run to produce BAM (for BAM header): ${pgFile}");

# Read whole file and build up an ordered list of the PG commands contained.
# When there are duplicate commands with the same ID, the latest one is used.
# To build the list, we first store PG commands by their IDs in a hash.
my $retval = 0;
my $pgCommandNumber = 1;
my %pgCommandsHash;
my $genomeRefFilePGLine = "\@PG\tID:ref_genome_indexing\tVN:${bwa_ver}\tCL:bwa index -a bwtsw ${ref_genome}\n";
$pgCommandsHash{"ref_genome_indexing"} = $genomeRefFilePGLine."\t0";

while(<PG_FILE>)
{
    chomp();
    my @line = split(/\t/);
    if ($line[0] ne "\@PG" or $_ !~ /VN:[ -~]+/ or $_ !~ /ID:[ -~]+/ or $_ !~ /CL:[ -~]+/)
    {
	dieScript("ERROR:  The PG file should only contain lines of the format: \@PG     ID:<cmdID>      VN:<version>    [PP:<prevPG>]   CL:<cmdline>. Found line:\n$_\n");
    }
    # Store the PG line as well as the index (order position) of the command
    ($line[1] =~ /ID:([ -~]+)/) and $pgCommandsHash{$1} = $_."\t".$pgCommandNumber;
    $pgCommandNumber++;
}
close(PG_FILE);

# Get keys to the PG lines, ordered by command number
my @sortedKeys = sort {
    my @a1 = split(/\t/, $pgCommandsHash{$a});
    my @b1 = split(/\t/, $pgCommandsHash{$b});
                ($a1[$#a1] eq $b1[$#b1]) and dieScript("FAILED ASSERT: Two PG command lines have the same command number - this should be impossible.\n");
    return ($a1[$#a1] <=> $b1[$#b1]);
} keys %pgCommandsHash;

# For each command that has a PP tag, verify that the command number of the
# previous program run is less than the current one. This ensures that we don't
# put an invalid order of PG tags in the header. For example, if the pipeline were
# restarted from a particular command, this would ensure that it was run to
# completion after that command in order to produce a valid BAM.
for (my $i = 0; $i <= $#sortedKeys; $i++)
{
    my $curLine = $pgCommandsHash{$sortedKeys[$i]};
    if ($curLine =~ /PP:([ -~]+)/)
    {
	my @curItems = split(/\t/, $curLine);
	my $curCmdNum = $curItems[$#curItems];
	my $ppID = $1;
	(exists $pgCommandsHash{$ppID}) or dieScript("ERROR: The PG file has a PP tag referencing ID \"$ppID\", but no PG line with this ID was found.\n");

	my @ppLine = split(/\t/, $pgCommandsHash{$ppID});
	my $ppCmdNum = $ppLine[$#ppLine];
	if ($curCmdNum <= $ppCmdNum)
	{
	    dieScript("ERROR: The PG line with ID \"$sortedKeys[$i]\" has a command number <= its 'previous' command with ID \"$ppID\".\n");
	}
    }
}

# Get the current BAM file header and remove all existing PG lines
#my $headerFile = "${noArchiveSampleFolder}/${sampleName}.header.sam";
runCommand("samtools view -H ${inputFile} | grep -v \"\@PG\" - > ${headerFile}");

# Write out the PG lines in descending order, most recent PG first, as required by FORGE
open(BAM_HEADER_FILE, ">>", $headerFile) or dieScript("Failed to open temporary header file for recording programs run to produce BAM (for BAM header): ${headerFile}");
for (my $i = $#sortedKeys; $i >= 0; $i--)
{
    my @line = split(/\t/, $pgCommandsHash{$sortedKeys[$i]});
    # Strip off the command number and print the line
    print BAM_HEADER_FILE join ("\t", @line[0..$#line-1])."\n";
}
close(BAM_HEADER_FILE);

#my $finalOutputBamName = "${sampleName}.t30l30.realigned.markdup.sorted.fixmate";
runCommand("samtools reheader ${headerFile} ${inputFile} > ${outputFile}");
runCommand("samtools index ${outputFile}");

# Calculate MD5 checksum for BAM file
runCommand("md5sum --binary ${outputFile} > ${outputFile}.md5");


###############################################################################

sub runCommand($)
{
    ($retval = system($_[0])) and dieScript("\nSCRIPT DIED: Error code $retval running system command: $_[0]\n");
    return 0;
}

sub dieScript($)
{
    print "This script failed with the error below. Note that more details (error output from the command) may be available in the log file/output stream.\n";
    print STDERR $_[0];
    die $_[0];
}

