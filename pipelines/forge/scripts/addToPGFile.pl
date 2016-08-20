#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;
use Sys::Hostname;
use Cwd 'abs_path';

sub notifyError($);
sub dieScript($);


# Set so that files created by the pipeline will be writable by group members.
# Otherwise it becomes impossible for more than one person to work on a sample,
# since all files will be writable only by them.
umask 0002;

# Flush any output as soon as it is written
$| = 1;

## OPTIONS
my ($cmdID, $cmdVer, $cmdLine, $cmdPP, $pgFile);

GetOptions(
    'cmdID=s'                  => \$cmdID,
    'cmdVer=s'                 => \$cmdVer,
    'cmdLine=s'                => \$cmdLine,
    'cmdPP=s'                  => \$cmdPP,
    'pgFile=s'                 => \$pgFile
);

my $ppSection = "";
($cmdPP) and $ppSection = "PP:${cmdPP}\t";
my $newPgLine = "\@PG\tID:${cmdID}\tVN:${cmdVer}\t${ppSection}CL:${cmdLine}\n";
if (! -e $pgFile)
{
    open(PG_FILE, ">", $pgFile) or dieScript("Failed to open file for recording programs run to produce BAM (for BAM header): ${pgFile}");
    print PG_FILE $newPgLine;
    close(PG_FILE);
} else {
    # Just add the new line to the file, regardless whether one with the same ID was
    # already there.
    open(PG_FILE, ">>", $pgFile) or dieScript("Failed to open file for recording programs run to produce BAM (for BAM header): ${pgFile}");
    print PG_FILE $newPgLine;
    close(PG_FILE);
}


###############################################################################

sub dieScript($)
{
    print "This script failed with the error below. Note that more details (error output from the command) may be available in the log file/output stream.\n";
    print STDERR $_[0];
    die $_[0];
}
