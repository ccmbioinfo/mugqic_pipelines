#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;

my ($read_files, $read_type) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

my $IDfile;
my $infile;
my $infile2;
my $IDfile2;
if ($read_type eq 'rRNA') {
    $IDfile = $Workpath.$read_files."1_qual_all_unique_IDs.txt";
    $infile = $Workpath.$read_files."1_qual_all_unique.fastq";
    $IDfile2 = $Workpath.$read_files."2_qual_all_unique_IDs.txt";
    $infile2 = $Workpath.$read_files."2_qual_all_unique.fastq";
} elsif ($read_type eq 'singletons') {
    $infile = $Workpath.$read_files."1_singletons.fastq";
    $infile2 = $Workpath.$read_files."2_singletons.fastq";
} elsif ($read_type eq 'micro_cds_sub') {
    $infile = $Workpath."microbial_cds_sub.fasta";
} elsif ($read_type eq 'nr_sub') {
    $infile = $Workpath."nr_all_sub.fasta";
}

my ($root, $suff) = $infile =~ m/(.*)\.(.*)/;

my $outfile = $root."_IDs_length.txt";
print "Outputfile: $outfile\n\n";

my %reads;
if (defined $IDfile) {
    open(INPUT, $IDfile) or die "Error opening $IDfile : $!\n";
    while( my $line = <INPUT> ) {
        chomp($line);
        my @values = split(/\t/, $line);
        $reads{$values[0]} = $values[1];
    }
    close INPUT;
    print "number of unique reads:  ".keys( %reads ).".\n\n\n";
}

unlink($outfile);
open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
my $seqIO = Bio::SeqIO->new( -file => $infile, -format => $suff );
while (my $seqobj = $seqIO->next_seq) {
    my $headf = $seqobj->description;
    my $head = $seqobj->display_id;
    my $seq_length = $seqobj->length;
    if ($read_type eq 'rRNA') {
        print OUTPUT "$head\t$reads{$head}\t$seq_length\n";
    } elsif ($read_type eq 'micro_cds_sub') {
        print OUTPUT "$head\t$seq_length\n";
    } elsif ($read_type eq 'nr_sub') {
        print OUTPUT "$head$headf\t$head\t$seq_length\n";
    } else {
        print OUTPUT "$head\t1\t$seq_length\n";
    }
}
close OUTPUT;

if (defined $infile2) {
    my ($root2, $suff2) = $infile2 =~ m/(.*)\.(.*)/;
    my $outfile2 = $root2."_IDs_length.txt";
    print "Outputfile: $outfile2\n\n";

    my %reads;
    if (defined $IDfile) {
        open(INPUT2, $IDfile2) or die "Error opening $IDfile2 : $!\n";
        while( my $line = <INPUT2> ) {
            chomp($line);
            my @values = split(/\t/, $line);
            $reads{$values[0]} = $values[1];
        }
        close INPUT2;
        print "number of unique reads:  ".keys( %reads ).".\n\n\n";
    }

    unlink($outfile2);
    open(OUTPUT2, '>>', $outfile2) or die "Error opening $outfile2 : $!\n";
    my $seqIO = Bio::SeqIO->new( -file => $infile2, -format => $suff2 );
    while (my $seqobj = $seqIO->next_seq) {
        my $head = $seqobj->display_id;
        my $seq_length = $seqobj->length;
        if ($read_type eq 'rRNA') {
            print OUTPUT2 "$head\t$reads{$head}\t$seq_length\n";
        } else {
            print OUTPUT2 "$head\t1\t$seq_length\n";
        }
    }
    close OUTPUT2;
}

exit 0;
