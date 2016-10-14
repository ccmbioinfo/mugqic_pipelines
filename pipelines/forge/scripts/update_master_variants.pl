#!/usr/bin/env perl
#

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;
use Cwd;
use Cwd 'abs_path';

sub dieScript($);
sub runCommand($);
sub updateVCFList($$$);
sub annovarPipeline($);
sub trim($);


my $scriptPath = abs_path($0);
my ($junkFilename, $SCRIPTS_DIR, $junkSuffix) = fileparse($scriptPath);
chop($SCRIPTS_DIR);  # remove trailing "/"

my $scriptName = "update_master_variants.pl";
my $sampleName;

my ($sampleNameList, $sampleConfig, $IN_VCF_FOLDER, $vardb, $vcflist, $VCF_FOLDER, $skipBuildingFile, $rebuildVardb, $verbose);
GetOptions(	'sampleNames=s' => \$sampleNameList,
			'sampleConfig=s' => \$sampleConfig,
		        'in_vcf_fol=s' => \$IN_VCF_FOLDER,
			'vardb=s' => \$vardb,
			'vcflist=s' => \$vcflist,
		        'vcffolder=s' => \$VCF_FOLDER,
			'skipBuildingFile' => \$skipBuildingFile,
			'rebuildVardb' => \$rebuildVardb,
			'v|verbose' => \$verbose) or &usage("");
$sampleNameList or &usage("ERROR: Missing sampleNames argument.\n");
$sampleConfig or &usage("ERROR: Missing sampleConfig argument.\n");


my @sampleNames = split(/,/, $sampleNameList);
# Trim whitespace from sample names
for (my $i = 0; $i <= $#sampleNames; $i++)
{
	$sampleNames[$i] = trim($sampleNames[$i]);
}

my %samplePhenotypes;
my %sampleIsCancer;

my $retval = 0;
my $sc = "";
my $replaceID;

my $commandLine = "$scriptName --sampleNames $sampleNameList --sampleConfig $sampleConfig";
$vardb and $commandLine .= " --vardb $vardb";
$vcflist and $commandLine .= " --vcflist $vcflist";
$skipBuildingFile and $commandLine .= " --skipBuildingFile";
$rebuildVardb and $commandLine .= " --rebuildVardb";
$verbose and $commandLine .= " -v";
print "\n######################### UPDATING MASTER LIST OF VARIANTS ########################\n";
print "Command line:\n$commandLine\n";



open(SAMPLE_INFO, "<", $sampleConfig) or die "Failed to open sample configuration file $sampleConfig: $!\n";
while (<SAMPLE_INFO>)
{
	chomp();
	(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
	(/^\s*$/) and next; #ignore white space lines
	my ($sampleNameConfig, $keyValuePair) = split(/\t/, $_, -1);
	$sampleName = $sampleNameConfig;

	my ($arg, $value) = split(/=/, $keyValuePair);
	if ($arg eq "PHENOTYPE")
	{
		$samplePhenotypes{$sampleNameConfig} = $value;
	}
	elsif ($arg eq "CANCER")
	{
		$sampleIsCancer{$sampleNameConfig} = $value;
	}
	elsif ($arg eq "VARDB")
	{
		$vardb = $value;
		(! -e $vardb) and dieScript("Missing variants DB: ${vardb}\n");
	}
}
close(SAMPLE_INFO);


foreach my $sampleName (@sampleNames)
{
	if (!exists $sampleIsCancer{$sampleName} and !exists $sampleIsCancer{"ALL"})
	{
		dieScript("No sample CANCER field specified for sample $sampleName in config file $sampleConfig.\n");
	}
}

chdir $VCF_FOLDER;

foreach my $sampleName (@sampleNames)
{
	$verbose and print "\n##### Sample ${sampleName}: updating master list of variants\n";
	runCommand("date");

	(-e "${VCF_FOLDER}/${sampleName}.flt.vcf") and print STDERR "Variants file ${VCF_FOLDER}/${sampleName}.flt.vcf already exists. OVERWRITING this file.\n";
	runCommand("cp ${IN_VCF_FOLDER}/${sampleName}.flt.vcf ${VCF_FOLDER}/${sampleName}.flt.vcf");
	updateVCFList("${sampleName}.flt.vcf", $sampleName, $vcflist);

}

if (!$skipBuildingFile)
{
	# First make a backup of the all variants file
	runCommand("cp $vardb $vardb.backup");

	my $options = "";
	(!$rebuildVardb) and $options .= " --existingdb $vardb.backup";
	my $cmd = "${SCRIPTS_DIR}/build_variants_db.pl $options --vcflist $vcflist 2>&1 > $vardb";

	print "###### Building variants database:\n";
	runCommand($cmd);
}

runCommand("date");
print "##### DONE updating master list of variants.\n";



###############################################################################

sub runCommand($)
{
	$verbose and !($_[0] eq "date") and print "\n$_[0]\n";
	($retval = system($_[0])) and dieScript("\nERROR code $retval running system command: $_[0]\n");
}

sub dieScript($)
{
	my $errorMsg = $_[0];
	die $errorMsg;
}

sub updateVCFList($$$)
{
	my $newSampleVariantsFile = $_[0];
	my $newSampleID = $_[1];
	my $allVariantsListFile = $_[2];
	my $addVariantFile = 1;

	# First read in the whole file to make sure this sample isn't already there.
	open(VARIANT_FILES_LIST, "+<", $allVariantsListFile);
	while (<VARIANT_FILES_LIST>)
	{
		chomp();
		(/^\s*#/) and next; #ignore comment lines (first non-ws character is #)
		(/^\s*$/) and next; #ignore white space lines
		my ($includeSample, $sampleID) = split(/\t/);

		if (lc($newSampleID) eq lc(trim($sampleID)))
		{
			# If the sample is already in the list then don't change anything.
			print "SampleID \"$newSampleID\" is already in variants file $allVariantsListFile. No need to update file.\n";
			$addVariantFile = 0;
			last;
		}
	}

	if ($addVariantFile)
	{
		my $isCancer = $sampleIsCancer{$newSampleID};
		(!$isCancer) and $isCancer = $sampleIsCancer{"ALL"};

		my $phenotype = $samplePhenotypes{$newSampleID};
		(!$phenotype) and $phenotype = $samplePhenotypes{"ALL"};

		# We didn't find this sampleID in the file, so let's add it.
		print "Adding sample to variants file $allVariantsListFile: $newSampleVariantsFile\t$newSampleID\n";
		print VARIANT_FILES_LIST "\n1\t$newSampleID\t$newSampleVariantsFile\t$isCancer\t$phenotype";
	}
	close(VARIANT_FILES_LIST);
}

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}


###############################################################################

sub usage
{
	print $_[0]."\n";
    print STDERR <<EOF;
Usage:
update_master_variants.pl --sampleNames NAME --sampleConfig FILE [options] > output.txt

This script takes a list of sample names to be added to the standard exome variants database and
does a few things: copies the sample.flt.vcf file to the all_variants/vcf folder; updates the
sample list used to build the variants database; makes a backup and then rebuilds the variants database.

Parameters:
--sampleNames Name1,Name2,... Comma-separated list of sample names; identify a sample's info in the config file.
--sampleConfig FILE           Path to the config file specifying the sample output dir and phenotype.
--vardb FILE                  The existing variants DB to update.
--vcflist FILE                The VCFlist to update, which is used to build the variants DB.
--vcffolder FOLDER            The folder where the VCF files are located (and will be placed).
--skipBuildingFile            Only copy files and add to the VCFlist, do not rebuild the variants DB.
--rebuildVardb                Rebuild the variants DB from scratch, rather than updating just with the new sample(s).
EOF

    exit;
}
