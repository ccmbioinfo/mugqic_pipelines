#!/usr/local/bin/perl -w

my ($read_files) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

my $infile1 = $Workpath."microbial_cds_sub_ecoli_ppi_pairs.txt";
my $infile2 = $Workpath."nr_all_sub_ecoli_ppi_pairs.txt";
my $outfile = $Workpath."PPI_pairs.txt";
print "Outputfile: $outfile\n\n";

my %genes;
open(INPUT1, $infile1) or die "Error opening $infile1 : $!\n";
while( my $line = <INPUT1> ) {
    chomp($line);
    my @values = split(/\t/, $line);
    if (exists $genes{$values[0]}) {
        my $tmp = join(';', $genes{$values[0]}, $values[1]);
        $genes{$values[0]} = $tmp;
    } else {
        $genes{$values[0]} = $values[1];
    }
}
close INPUT1;
print "number of mapped genes:  ".keys( %genes ).".\n";

open(INPUT2, $infile2) or die "Error opening $infile2 : $!\n";
while( my $line = <INPUT2> ) {
    chomp($line);
    my @values = split(/\t/, $line);
    if (exists $genes{$values[0]}) {
        my $tmp = join(';', $genes{$values[0]}, $values[1]);
        $genes{$values[0]} = $tmp;
    } else {
        $genes{$values[0]} = $values[1];
    }
}
close INPUT2;
print "number of mapped genes+proteins:  ".keys( %genes ).".\n\n\n";

unlink($outfile);
open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
foreach $key (sort keys %genes) {
    print OUTPUT $key."\t".$genes{$key}."\n";
}
close OUTPUT;

exit 0;

