#!/usr/local/bin/perl -w


my ($infile, $infile2, $outfile, $outfile2, $BLASTDB) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

#my $infile = $Workpath."nr_all_sub_IDs_length.txt";
#my $infile2 = $Workpath."nr_all_sub_IDs_giID.txt";
#my $outfile = $Workpath."nr_all_sub_IDs_map_taxid.txt";
#my $outfile2 = $Workpath."nr_all_sub_IDs_taxonID.txt";
#print "Outputfile: $outfile\n\n";

my %proteinsID;
unlink($infile2);
open(OUTPUT, '>>', $infile2) or die "Error opening $infile2 : $!\n";
open(INPUT, $infile) or die "Error opening $infile : $!\n";
while( my $line1 = <INPUT> ) {
    chomp($line1);

    next if $line1 =~ /^\s*@/;

    my @line = split(/\t/, $line1);

    my $giID1 = $line[0];
    $giID1 =~ s/^gi\|//;
    $giID1 =~ s/\|.*?$//;

#    my @tmp = split(/\|/, $line[0]);
#    my $accession = $tmp[3];

    my $spName1 = $line[1];
    $spName1 =~ s/.*?\[//;
    $spName1 =~ s/\]//;

    $proteinsID{$.} = join("\t", $line[0], $giID1, $spName1);

    print OUTPUT "$giID1\n";
}
close INPUT;
close OUTPUT;
print "number of proteinsID:  ".keys( %proteinsID ).".\n\n\n";

system("blastdbcmd -db $BLASTDB -dbtype prot -entry_batch ".$infile2." -outfmt %T -target_only -out ".$outfile2);

unlink($outfile);
open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
open(INPUT, $outfile2) or die "Error opening $outfile2 : $!\n";
while( my $taxID1 = <INPUT> ) {
    chomp($taxID1);
    print OUTPUT "$proteinsID{$.+ 1}\t$taxID1\n";
}
close INPUT;
close OUTPUT;

exit 0;
