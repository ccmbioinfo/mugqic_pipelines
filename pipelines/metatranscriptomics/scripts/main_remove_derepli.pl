#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;

my ($read_files) = @ARGV;
#my $read_files="cow";

my ($IDfile, $infile, $outfile) = @ARGV;

#my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
#my $Workpath = "";

#for (my $i = 1; $i <= 2; $i++) {
#    my $IDfile = "remove_duplicates/cow".$i."_qual_all_unique_IDs.txt";
#    my $IDfile = "remove_duplicates/cow".$i."_qual_all_unique_IDs_length.txt";
#    my $infile = "remove_host_reads/cow".$i."_qual_unique_n_rRNA_n_host.fastq";
#    my $outfile = "add_duplicates/cow".$i."_mRNA.fastq";

#    my $IDfile = $Workpath.$read_files.$i."_qual_all_unique_IDs_length.txt";
#    my $infile = $Workpath.$read_files.$i."_qual_unique_n_rRNA_n_host.fastq";
#    my $outfile = $Workpath.$read_files.$i."_mRNA.fastq";
    print "Outputfile: $outfile\n\n";

    my $total;
    my %IDs;
    open(INPUT0, $IDfile) or die "Error opening $IDfile : $!\n";
    while  (my $line = <INPUT0>) {
        chomp($line);
        my @line1 = split(/\t/, $line);
        $IDs{$line1[0]} = $line1[1];
        $total += $line1[1];
    }
    close INPUT0;
    print "number of unqiue reads:  ".keys( %IDs ).".\n";
    print "number of reads:  ".$total.".\n\n\n";

    unlink($outfile);
    open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
    open(INPUT, $infile) or die "Error opening $infile : $!\n";
    while  (my $head = <INPUT>) {
        chomp($head);
        my $tmp = substr($head, 1);

        my $seq = <INPUT>;
        chomp($seq);

        my $head2 = <INPUT>;
        chomp($head2);

        my $qual = <INPUT>;
        chomp($qual);

        print OUTPUT "$head\n$seq\n$head2\n$qual\n";

        for (my $j = 2; $j <= $IDs{$tmp}; $j++) {
            my $head_m;
            if ($head =~ m/\//) {
                my @tmp2 = split(/\//, $head);
                $head_m = $tmp2[0].'_'.$j.'/'.$tmp2[1];
            } else {
                $head_m = $head.'_'.$j;
            }
            print OUTPUT "$head_m\n$seq\n$head2\n$qual\n";
        }
    }
    close INPUT;
    close OUTPUT;
#}

exit 0;
