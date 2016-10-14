#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

sub handleVariantFunctionLine($$$$);
sub replaceDelimiters($);

# Note: Left the ability to deal with the old way (not table_annovar) in case an older version of annovar needs to be used

(@ARGV == 0) and usage("");

my ($vcfFile, $exonVarsFile, $allVarsFile, $allVarsExtendedFile, $dbsnpFile, $thgFile, $consFile, $tableFile, $siftFile, $polyphenFile, $mtFile, $lrtFile, $gerpFile, $dgvFile, $segdupFile, $exacFile, $clinvarFile, $includeAll, $includeSynonymous, $includeIntronic, $includeUTR, $help, $verbose);
GetOptions(	'vcf=s' => \$vcfFile,			#the original VCF file that was annotated - needed only to retain the VCF header
			'exonvars=s' => \$exonVarsFile, #annovar_out.exonic_variant_function
			'allvars=s' => \$allVarsFile,   #annovar_out.variant_function (use to check if there are any splice site SNPs)
			'allvarsExtended=s' => \$allVarsExtendedFile,   #annovar_out.variant_function.splicing_extended
		        'table=s' => \$tableFile,       #annovar_out.hg19_multianno.txt
			'dbsnp=s' => \$dbsnpFile,       #annovar_out.hg19_snp131_dropped
			'thg=s' => \$thgFile,           #annovar_out.hg19_ALL.sites.2010_11_dropped
			'cons=s' => \$consFile,         #annovar_out.hg19_mce46way
			'sift=s' => \$siftFile,         #annovar_out.hg19_ljb_sift_dropped
			'polyphen=s' => \$polyphenFile, #annovar_out.asdf
			'mt=s' => \$mtFile,             #annovar_out.hg19_ljb_mt_dropped
			'lrt=s' => \$lrtFile,           #annovar_out.hg19_ljb_lrt_dropped
			'gerp=s' => \$gerpFile,         #annovar_out.hg19_ljb_gerp+++_dropped
			'dgv=s' => \$dgvFile,           #annovar_out.hg19_dgv
			'segdup=s' => \$segdupFile,     #annovar_out.hg19_genomicSuperDups
			'exac=s' => \$exacFile,		#annovar_out.hg19_exac03
			'clinvar=s' => \$clinvarFile,	#annovar_out.hg19_clinvar
			'includeAll' => \$includeAll,
			'includeSynonymous' => \$includeSynonymous,
			'includeIntronic' => \$includeIntronic,
			'includeUTR' => \$includeUTR,
			'help|h' => \$help,
			'v|verbose' => \$verbose) or usage("");

$help and usage("");
$vcfFile or usage("Missing required parameter: vcf. You must specify the vcf file that was annotated (its header is needed).");
$exonVarsFile or usage("Missing required parameter: exonvars. You must specify the annovar exonic variant function output file.");
$allVarsFile or usage("Missing required parameter: allvars. You must specify the annovar variant function output file.");
# $allVarsExtendedFile or usage("Missing required parameter: allvarsExtended. You must specify the annovar variant function output file.");

open(VCF_FILE, "<".$vcfFile) or die("Failed to open vcf file: ${vcfFile}");
open(CODING_FILE, "<".$exonVarsFile) or die("Failed to open exonic_variant_function file: ${exonVarsFile}");
open(ALL_FILE, "<".$allVarsFile) or die("Failed to open variant_function file: ${allVarsFile}");

$tableFile and (open(TABLE_FILE, "<".$tableFile) or die("Failed to open table file: ${tableFile}"));
$dbsnpFile and (open(DBSNP_FILE, "<".$dbsnpFile) or die("Failed to open dbSNP file: ${dbsnpFile}"));
$thgFile and (open(THG_FILE, "<".$thgFile) or die("Failed to open 1000 genomes file: ${thgFile}"));
$consFile and (open(CONS_FILE, "<".$consFile) or die("Failed to open conservation file: ${consFile}"));
$siftFile and (open(SIFT_FILE, "<".$siftFile) or die("Failed to open SIFT file: ${siftFile}"));
$polyphenFile and (open(POLYPHEN_FILE, "<".$polyphenFile) or die("Failed to open Polyphen file: ${polyphenFile}"));
$mtFile and (open(MT_FILE, "<".$mtFile) or die("Failed to open SIFT file: ${mtFile}"));
$lrtFile and (open(LRT_FILE, "<".$lrtFile) or die("Failed to open SIFT file: ${lrtFile}"));
$gerpFile and (open(GERP_FILE, "<".$gerpFile) or die("Failed to open SIFT file: ${gerpFile}"));
$dgvFile and (open(DGV_FILE, "<".$dgvFile) or die("Failed to open DGV file: ${dgvFile}"));
$segdupFile and (open(SEGDUP_FILE, "<".$segdupFile) or die("Failed to open Segdup file: ${segdupFile}"));


my %variantLines; #text line with VCF details
my %variantInfoToAdd; #text to add to Info field for variant
my %proteinChanges;
my $lineNum = 0;
my $numCodingVarsSkipped = 0;
while (<CODING_FILE>)
{
	$lineNum++;
	chomp();
	my @line = split(/\t/);
	if (scalar(@line) > 50)
	{
		$verbose and print STDERR "ERROR: line $. does not have expected number of fields in file $exonVarsFile. Skipping this variant.\n";
		$numCodingVarsSkipped++;
		next;
	}
	if ($includeAll or $includeSynonymous or ($line[1] !~ /^synonymous SNV$/ and $line[1] !~ /^exonic$/))
	{
		my @details = split(/:/, $line[2]);
		(my $geneName) = split(/[;:(,]/, $details[0]);
		my $varKey = "$line[3]:$line[4]-$line[5];$line[6]>$line[7]";
		if ($line[8] !~ /[chr]?[\dXYM]+/)
		{
			print STDERR "WARNING: Likely Annovar annotation error. Empty or bad value ($line[8]) found for column 9 of line $. in file ${exonVarsFile}.\n";
			next;
		}
		# The VCF portion includes everything from col 4 onward
		$variantLines{$varKey} = join("\t", @line[8..$#line]);

		$variantInfoToAdd{$varKey} = ";VT=".replaceDelimiters($line[1]).";GENE=".$geneName.";DTLS=".replaceDelimiters($line[2]);
		my @detailBits = split(":", $line[2]);
		if (@detailBits > 4 && substr($detailBits[4], 0, 1) eq "p")
		{
			#NOTE - the location of the protein change info differs between Annovar versions!!
			my @pChange = split(",",$detailBits[4]);
			$variantInfoToAdd{$varKey} .= ";PC=".$pChange[0];
		}

		my $dp4 = "";
		my ($altCount, $readCount);
		my $vcfInfo = $line[15];
		($vcfInfo =~ /DP4=([^;\t]*)/) and $dp4 = $1;
		if ($dp4)
		{
			my ($fwdRef, $revRef, $fwdAlt, $revAlt) = split(/,/, $dp4);
			$altCount = $fwdAlt + $revAlt;
			$readCount = $altCount + $fwdRef + $revRef;
		} else {
			# See if there is an AD (allelic depth) field in the per-sample genotype info,
			# which is provided if variants are called by GATK.
			if ($line[16])
			{
				my @fieldCodes = split(/:/, $line[16]);
				my $adFieldPos = -1;
				for (my $i = 0; $i <= $#fieldCodes; $i++)
				{
					if ($fieldCodes[$i] =~ /^AD$/i)
					{
						$adFieldPos = $i;
						last;
					}
				}
				if ($adFieldPos >= 0)
				{
					my @fields = split(/:/, $line[17]);
					my $refCount;
					($refCount, $altCount) = split(/,/, $fields[$adFieldPos]);
					$readCount = $altCount + $refCount;
				}
			}
		}
		if (!defined $altCount or !defined $readCount)
		{
			($verbose) and print STDERR "Warning: no DP4 field found for line ${lineNum} of file ${exonVarsFile}\n";
			next;
		}
		$variantInfoToAdd{$varKey} .= ";ALTC=".$altCount;
		$variantInfoToAdd{$varKey} .= ";RDC=".$readCount;
	}
}
close(CODING_FILE);
($numCodingVarsSkipped > 0) and print STDERR "WARNING: $numCodingVarsSkipped variant lines skipped which did not have expected number of fields in file $exonVarsFile.\n";

my $numExtendedVarsSkipped = 0;
$lineNum = 0;
while (<ALL_FILE>)
{
	$lineNum++;
	chomp();
	handleVariantFunctionLine($_, $lineNum, $allVarsFile, 0);
}
close(ALL_FILE);


if ($allVarsExtendedFile)
{
	open(EXTENDED_FILE, "<".$allVarsExtendedFile) or die("Failed to open splicing_extended.variant_function file: ${allVarsExtendedFile}");
	$lineNum = 0;
	while (<EXTENDED_FILE>)
	{
		$lineNum++;
		chomp();
		handleVariantFunctionLine($_, $lineNum, $allVarsExtendedFile, 1);
	}
	close(EXTENDED_FILE);
	($numExtendedVarsSkipped > 0) and print STDERR "WARNING: $numExtendedVarsSkipped variant lines skipped which did not have expected number of fields in two 'extended' variants files.\n";
}

sub handleVariantFunctionLine($$$$)
{
	my @line = split(/\t/, $_[0]);
	if (scalar(@line) > 50)
	{
		$verbose and print STDERR "ERROR: line $. does not have expected number of fields in file $_[2]. Skipping this variant.\n";
		$numExtendedVarsSkipped++;
		return;
	}
	my $lineNum = $_[1];
	my $isExtended = $_[3];
	my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
	my $addVariant = 0;

	# We only overwrite variants from the exonic vars file in special cases
	if ($includeAll or
		($line[0] =~ /^splicing/) or
		($line[0] =~ /exonic;splicing/) or
		($includeIntronic and $line[0] =~ /intronic/) or
		($includeUTR and $line[0] =~ /^UTR/))
	{
		$addVariant = 1;
	}
	if (exists $variantLines{$varKey})
	{
		# We only overwrite an existing variant if the new annotation appears more deleterious.
		if (! ($line[0] =~ /^splicing/ or ($line[0] =~ /exonic;splicing/ and $variantInfoToAdd{$varKey} =~ /VT=(intronic|^synonymous)/)))
		{
			$addVariant = 0;
		}
	}
	if ($addVariant)
	{
		if ($line[7] !~ /[chr]?[\dXYM]+/)
		{
			print STDERR "WARNING: Likely Annovar annotation error. Empty or bad value ($line[7]) found for column 8 of line $. in extended variants file $_[2].\n";
			return;
		}
		my $variantType = replaceDelimiters($line[0]);
		if ($isExtended)
		{
			# The only case where we want to add a variant from the "extended" file is when it is
			# a "splicing" variant, and it wasn't already present (or it was annotated previously
			# as intronic or synonymous), so it is a "splicing-extended" variant.
			if ($variantType =~ /splicing/ and (!exists $variantLines{$varKey} or $variantInfoToAdd{$varKey} =~ /VT=(intronic|synonymous)/))
			{
				# The new annotation is stronger, so use it; note that it is an "extended" splicing annotation
				$variantType .= "-extended";
			} else {
				return; # Variant was already processed and the new one is not splicing, so skip it
			}
		}

		$variantLines{$varKey} = join("\t", @line[7..$#line]);
		my @info = split(/:/, $line[1]);
		(my $geneName) = split(/[;:(,]/, $info[0]);
		$variantInfoToAdd{$varKey} = ";VT=${variantType};GENE=${geneName};DTLS=$line[1]";

		my $dp4 = "";
		my ($altCount, $readCount);
		my $vcfInfo = $line[14];
		($vcfInfo =~ /DP4=([^;\t]*)/) and $dp4 = $1;
		if ($dp4)
		{
			my ($fwdRef, $revRef, $fwdAlt, $revAlt) = split(/,/, $dp4);
			$altCount = $fwdAlt + $revAlt;
			$readCount = $altCount + $fwdRef + $revRef;
		} else {
			# See if there is an AD (allelic depth) field in the per-sample genotype info,
			# which is provided if variants are called by GATK.
			if ($line[15])
			{
				my @fieldCodes = split(/:/, $line[15]);
				my $adFieldPos = -1;
				for (my $i = 0; $i <= $#fieldCodes; $i++)
				{
					if ($fieldCodes[$i] =~ /^AD$/i)
					{
						$adFieldPos = $i;
						last;
					}
				}
				if ($adFieldPos >= 0)
				{
					my @fields = split(/:/, $line[16]);
					my $refCount;
					($refCount, $altCount) = split(/,/, $fields[$adFieldPos]);
					$readCount = $altCount + $refCount;
				}
			}
		}
		if (!defined $altCount or !defined $readCount)
		{
			($verbose) and print STDERR "Warning: no DP4 field found for line ${lineNum} of file ${exonVarsFile}\n";
			next;
		}
		$variantInfoToAdd{$varKey} .= ";ALTC=".$altCount;
		$variantInfoToAdd{$varKey} .= ";RDC=".$readCount;
	}
}
# Now all mutations that we are going to process are in the %soi hash table

if ($dbsnpFile)
{
	while (<DBSNP_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantLines{$varKey})
		{
			# dbSNP is a bit special because VCF has a field ID (col 3) specially for it
			# So rather than putting it in the INFO field we put it there.
			my @varLine = split(/\t/, $variantLines{$varKey});
			$varLine[2] = $line[1];
			$variantLines{$varKey} = join("\t", @varLine);

			# Possibly a more efficient way of putting in the dbSNP ID?
			#$variantLines{$varKey} =~ s/([^\t]+\t[^\t]+\t)\./$1$line[1]/;
		}
	}
	close(DBSNP_FILE);
}

if ($thgFile)
{
	while (<THG_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			$variantInfoToAdd{$varKey} .= ";THGMAF=$line[1]";
		}
	}
	close(THG_FILE);
}


if ($tableFile)
{
        while(<TABLE_FILE>)
	{
	        chomp();
   	        my @line = split(/\t/);
	        my $varKey = "$line[0]:$line[1]-$line[2];$line[3]>$line[4]";
 	        if (exists $variantInfoToAdd{$varKey})
	        {
		    # dbSNP is a bit special because VCF has a field ID specially for it
		    # So rather than putting it in the INFO field we put it there.
		    my @varLine = split(/\t/, $variantLines{$varKey});
		    $varLine[2] = $line[6];
		    $variantLines{$varKey} = join("\t", @varLine);

		    $variantInfoToAdd{$varKey} .= ";THGMAF=$line[5]";

		    my $numScore = $line[7];
		    if (looks_like_number($numScore)) {
			$variantInfoToAdd{$varKey} .= ";SIFT=$numScore";
		    }

		    my $numScore = $line[9];
		    if (looks_like_number($numScore)) {
			$variantInfoToAdd{$varKey} .= ";PP2_HDIV=$numScore"; # Originally, just "PP2" -- let's se if this is okay
		    }

		    my $numScore = $line[11];
		    if (looks_like_number($numScore)) {
			$variantInfoToAdd{$varKey} .= ";PP2_HVAR=$numScore"; 
		    }

		    my $numScore = $line[13];
		    if (looks_like_number($numScore)) {
			$variantInfoToAdd{$varKey} .= ";LRT=$numScore";
		    }

		    my $numScore = $line[15];
		    if (looks_like_number($numScore)) {
			$variantInfoToAdd{$varKey} .= ";MT=$numScore";
		    }

		    my $numScore = $line[27];
		    if (looks_like_number($numScore)) {
			$variantInfoToAdd{$varKey} .= ";CADD_Phred=$numScore";
		    }

		    my $numScore = $line[28];
		    if (looks_like_number($numScore)) {
			$variantInfoToAdd{$varKey} .= ";GERP=$numScore";
		    }

		    my $numScore = $line[32];
		    if (looks_like_number($numScore)) {
			$variantInfoToAdd{$varKey} .= ";EXAC=$numScore";
		    }

		    # Modified:
		    # 	CLINVAR output had ";"-separated fields which could have caused problems when trying to separate other field, so here
		    # 	the separator of fields have been changed from ";" to ":" 
		    if ($line[33] ne ".") {
		        $line[33] =~ tr/;/:/;
		    	$variantInfoToAdd{$varKey} .= ";CLINVAR=$line[33]";
		    }

		    if ($line[34] ne ".") {
		        $line[34] =~ tr/;/:/;
		        $variantInfoToAdd{$varKey} .= ";PHC=$line[34]";
		    }
	        }
	}
	close(TABLE_FILE);
}

if ($consFile)
{
	while (<CONS_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			my @info = split(/;/, $line[1]);
			my ($junk1, $numScore) = split("=", $info[0]);
			$variantInfoToAdd{$varKey} .= ";PHC=$numScore";
		}
	}
	close(CONS_FILE);
}

if ($siftFile)
{
	while (<SIFT_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			my $numScore = $line[1];
			$variantInfoToAdd{$varKey} .= ";SIFT=$numScore";
		}
	}
	close(SIFT_FILE);
}

if ($polyphenFile)
{
	while (<POLYPHEN_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			my $numScore = $line[1];
			$variantInfoToAdd{$varKey} .= ";PP2=$numScore";
		}
	}
	close(POLYPHEN_FILE);
}

if ($mtFile)
{
	while (<MT_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			my $numScore = $line[1];
			$variantInfoToAdd{$varKey} .= ";MT=$numScore";
		}
	}
	close(MT_FILE);
}

if ($lrtFile)
{
	while (<LRT_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			my $numScore = $line[1];
			$variantInfoToAdd{$varKey} .= ";LRT=$numScore";
		}
	}
	close(LRT_FILE);
}

if ($gerpFile)
{
	while (<GERP_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			my $numScore = $line[1];
			$variantInfoToAdd{$varKey} .= ";GERP=$numScore";
		}
	}
	close(GERP_FILE);
}

if ($dgvFile)
{
	while (<DGV_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			my ($junk, $names) = split("=", $line[1]);
			$variantInfoToAdd{$varKey} .= ";DGV=$names";
		}
	}
	close(DGV_FILE);
}

if ($segdupFile)
{
	while (<SEGDUP_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey})
		{
			my @info = split(/;/, $line[1]);
			my ($junk1, $numScore) = split("=", $info[0]);
			$variantInfoToAdd{$varKey} .= ";SDS=$numScore";
		}
	}

	seek(SEGDUP_FILE, 0, 0);
	while (<SEGDUP_FILE>)
	{
		chomp();
		my @line = split(/\t/);
		my $varKey = "$line[2]:$line[3]-$line[4];$line[5]>$line[6]";
		if (exists $variantInfoToAdd{$varKey}) {
			my @info = split(/;/, $line[1]);
			my ($junk2, $name) = split("=", $info[1]);
			$variantInfoToAdd{$varKey} .= ";SDL=$name";
		}
	}
	close(SEGDUP_FILE);
}


my $varTypeHeader		= "##INFO=<ID=VT,Number=1,Type=String,Description=\"Type of variant (intronic, synonymous, nonsynonymous, etc)\">\n";
my $geneHeader			= "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene symbol\">\n";
my $detailsHeader		= "##INFO=<ID=DTLS,Number=1,Type=String,Description=\"Variant details from Annovar\">\n";
my $altCountHeader		= "##INFO=<ID=ALTC,Number=1,Type=String,Description=\"Alt read count\">\n";
my $readCountHeader		= "##INFO=<ID=RDC,Number=1,Type=String,Description=\"Total read count\">\n";
my $proteinChangeHeader	= "##INFO=<ID=PC,Number=1,Type=String,Description=\"Protein change\">\n";
my $thgHeader			= "##INFO=<ID=THGMAF,Number=1,Type=String,Description=\"Minor allele frequency of variant in thousand genomes\">\n";
my $phastConsHeader		= "##INFO=<ID=PHC,Number=1,Type=String,Description=\"Phastcons conservation score\">\n";
my $tableHeader                 = "##INFO=<ID=TABLE,Number=1,Type=String,Description=\"TABLE Scores - SIFT, Polyphen2 (HDIV + HVAR), MutationTaster, LRT, GERP, Phastcons\">\n";
my $siftHeader			= "##INFO=<ID=SIFT,Number=1,Type=String,Description=\"SIFT score\">\n";
my $polyphenHeader		= "##INFO=<ID=PP2,Number=1,Type=String,Description=\"Polyphen2 score\">\n";
my $mtHeader			= "##INFO=<ID=MT,Number=1,Type=String,Description=\"MutationTaster score\">\n";
my $lrtHeader			= "##INFO=<ID=LRT,Number=1,Type=String,Description=\"LRT score\">\n";
my $snvavgHeader		= "##INFO=<ID=SNVAVG,Number=1,Type=String,Description=\"SNV score average\">\n";
my $gerpHeader			= "##INFO=<ID=GERP,Number=1,Type=String,Description=\"GERP score\">\n";
my $dgvHeader			= "##INFO=<ID=DGV,Number=1,Type=String,Description=\"Database of genomic variants IDs\">\n";
my $segdupScoreHeader	= "##INFO=<ID=SDS,Number=1,Type=String,Description=\"Segmental duplication percent identity\">\n";
my $segdupLocHeader		= "##INFO=<ID=SDL,Number=1,Type=String,Description=\"Segmental duplication locus\">\n";


# Get VCF header from original file
my $infoHeaderSeen = 0;
my $newHeaderPrinted = 0;
while (<VCF_FILE>)
{
	if ($_ !~ /^#/)
	{
		# We are at the end of the header. We don't need any more of the file.
		last;
	}

	if (/^##INFO=/)
	{
		$infoHeaderSeen = 1;
	}
	elsif ($infoHeaderSeen && !$newHeaderPrinted)
	{
		# This is not an INFO header line but one has been seen, so we are just at
		# the end of the INFO header block. Add our INFO lines here.
		print $varTypeHeader;
		print $geneHeader;
		print $detailsHeader;
		print $altCountHeader;
		print $readCountHeader;
		print $proteinChangeHeader;
		($thgFile) and print $thgHeader;
		($consFile) and print $phastConsHeader;
		($tableFile) and print $tableHeader;
		($siftFile) and print $siftHeader;
		($polyphenFile) and print $polyphenHeader;
		($mtFile) and print $mtHeader;
		($lrtFile) and print $lrtHeader;
		($gerpFile) and print $gerpHeader;
		($dgvFile) and print $dgvHeader;
		($dgvFile) and print $dgvHeader;
		print $snvavgHeader;
		($segdupFile) and print $segdupLocHeader;
		$newHeaderPrinted = 1;
	}
	print $_;
}
close(VCF_FILE);

# sort keys by position
my @sortedVariantKeys = sort {
	my ($achr, $aStart, $bchr, $bStart);
	if ($a =~ /^(chr)?(X|Y|\d+):(\d*)-\d*;/)
	{
		$achr = $2;
		$aStart = $3;
		if ($b =~ /^(chr)?(X|Y|\d+):(\d*)-\d*;/)
		{
			$bchr = $2;
			$bStart = $3;
			($achr eq "X") and $achr = 23;
			($achr eq "Y") and $achr = 24;
			($achr eq "M") and $achr = 25;
			($bchr eq "X") and $bchr = 23;
			($bchr eq "Y") and $bchr = 24;
			($bchr eq "M") and $bchr = 25;

			($achr < $bchr) and return -1;
			($achr > $bchr) and return 1;
			return ($aStart <=> $bStart);
		} else {
			# chrB doesn't match so is probably a random chr
			return -1;
		}
	}
	elsif ($b =~ /^(chr)?(X|Y|\d+):(\d*)-\d*;/)
	{
		# chrA doesn't match but chrB does
		return 1;
	}
	#print STDERR "ERROR: comparing keys \"$a\" and \"$b\"\n";
	return $a cmp $b;
} keys %variantLines;


for (my $i=0; $i <= $#sortedVariantKeys; $i++)
{
	my $key = $sortedVariantKeys[$i];
	if (!exists $variantInfoToAdd{$key})
	{
		print STDERR "ERROR: Assertion failed: No variant info to add for variant with key $key\n";
		next;
	}
	my $infoToAdd = $variantInfoToAdd{$key};
	($infoToAdd !~ /THGMAF=/) and $infoToAdd .= ";THGMAF=0";

	my $snvAvg = 0;
	my $numSNVScores = 0;
	($infoToAdd =~ /SIFT=([^;\t]+)/) and $snvAvg += $1 and $numSNVScores++;
	($infoToAdd =~ /PP2=([^;\t]+)/)  and $snvAvg += $1 and $numSNVScores++;
	($infoToAdd =~ /MT=([^;\t]+)/)   and $snvAvg += $1 and $numSNVScores++;
	($infoToAdd =~ /LRT=([^;\t]+)/)  and $snvAvg += $1 and $numSNVScores++;
	if ($numSNVScores > 0)
	{
		$snvAvg = $snvAvg / $numSNVScores;
		$infoToAdd .= ";SNVAVG=$snvAvg";
	}

	my @line = split(/\t/, $variantLines{$key});
	$line[7] .= $infoToAdd;

	print join("\t", @line)."\n";
}


###############################################################################

sub usage
{
	print $_[0]."\n";
    print STDERR <<EOF;
This program is used to combine a set of Annovar annotation files together after
Annovar has been run. It produces a VCF file as output with new entries in the
INFO column. You must specify the original vcf file, as well as the required
Annovar output files specified below.

Usage:
$0 --sampleName NAME --sampleConfig FILE [options]

OPTIONS:
--vcf FILE              REQUIRED. The VCF file that was annotated
--exonvars FILE         REQUIRED. sample.annovar_out.exonic_variant_function
--allvars FILE          REQUIRED. sample.annovar_out.variant_function
--allvarsExtended FILE  REQUIRED. sample.annovar_out.splicing_extended.variant_function
--dbsnp FILE            sample.annovar_out.hg19_snp132_dropped
--thg FILE              sample.annovar_out.hg19_ALL.sites.2010_11_dropped
--cons FILE             sample.annovar_out.hg19_mce46way
--sift FILE             sample.annovar_out.hg19_avsift_dropped
--dgv FILE              sample.annovar_out.hg19_dgv
--segdup FILE           sample.annovar_out.hg19_genomicSuperDups
--includeSynonymous     synonymous variants will be included in the output
--includeIntronic       intronic variants will be included in the output
EOF
    exit;
}

###############################################################################

sub replaceDelimiters($)
{
    my $text = shift;
	# Replace any occurrences of ";" or "\t", which are delimiters in VCF
	$text =~ s/[;\t]/-/;
    return $text;
}
