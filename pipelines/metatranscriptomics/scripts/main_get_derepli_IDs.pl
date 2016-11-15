#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;

my ($read_files) = @ARGV;
#my $read_files='cow';

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

for (my $i = 1; $i <= 2; $i++) {

    my $IDfile = $Workpath.$read_files.$i."_qual_all_unique.fasta";
    my $infile = $Workpath.$read_files.$i."_qual_all.fastq";

    my $outfile0 = $Workpath.$read_files.$i."_qual_all_unique_IDs.txt";
    my $outfile = $Workpath.$read_files.$i."_qual_all_unique.fastq";
    print "Outputfiles: $outfile0\n $outfile\n\n";

    my %IDs;
    unlink($outfile0);
    open(OUTPUT0, '>>', $outfile0) or die "Error opening $outfile0 : $!\n";
    open(INPUT0, $IDfile) or die "Error opening $IDfile : $!\n";
    while  ( my $line = <INPUT0>) {
        chomp($line);
        if ($line =~ m/^>/) {
            my $line1 = $line;
            $line1 =~ s/;.*//;
            my $tmp = substr($line1, 1);
            $IDs{$tmp} = 1;

            my $line2 = $line;
            $line2 =~ s/.*?=//;
            my $tmp2 = substr($line2, 0, -1);

            print OUTPUT0 "$tmp\t$tmp2\n";
        }
    }
    close INPUT0;
    close OUTPUT0;

    my $seq_out = Bio::SeqIO->new( -file => ">$outfile", -format => 'fastq' );
    my $seq_in = Bio::SeqIO->new( -file => $infile, -format => 'fastq' );
    while (my $seqobj = $seq_in->next_seq) {
        my $head = $seqobj->display_id;
        if (exists $IDs{$head}) {
            $seq_out->write_seq( $seqobj );
        }
    }
}

exit 0;
