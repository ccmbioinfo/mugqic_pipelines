#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;

my ($read_files, $data_type) = @ARGV;
#my $read_files="cow";
#my $data_type="assembly";


my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $workpath = "";

my $IDfile;
my $IDfile2;
my $IDfile21;
my $contigs;
my $out_ids;
my $out_length;
if ($data_type eq 'assembly') {
    my $contig_dir = "index_contigs";
    my $output_dir = "map_reads";

    $IDfile = $output_dir."/".$read_files."_trinity_bwa_pairs.txt";
    $contigs = $contig_dir."/".$read_files."_contigs.fasta";

    $out_ids = $output_dir."/".$read_files."_contigs_IDs.txt";
    $out_length = $output_dir."/".$read_files."_contigs_IDs_length.txt";
    print "Outputfile: $out_length\n\n";
}
# Not necessary (see below)
#} elsif ($data_type eq 'bwa') {
#    $IDfile = $workpath.$read_files."_contigs_IDs.txt";
#    $IDfile2 = $workpath.$read_files."_contigs_micro_cds_bwa_IDs.txt";
#    $IDfile21 = $workpath.$read_files."_singletons_micro_cds_bwa_IDs.txt";
#} elsif ($data_type eq 'blat') {
#    $IDfile = $workpath.$read_files."_contigs_IDs.txt";
#    $IDfile2 = $workpath.$read_files."_contigs_n_micro_cds_blat_IDs.txt";
#    $IDfile21 = $workpath.$read_files."_singletons_n_micro_cds_blat_IDs.txt";
#} elsif ($data_type eq 'diamond') {
#    $IDfile = $workpath.$read_files."_contigs_IDs.txt";
#    $IDfile2 = $workpath.$read_files."_contigs_nr_diamond_IDs.txt";
#    $IDfile21 = $workpath.$read_files."_singletons_nr_diamond_IDs.txt";
#}

if ($data_type eq 'assembly') {
    my %id_to_num_reads;
    open(INPUT0, $IDfile) or die "Error opening $IDfile : $!\n";
    while  ( my $line = <INPUT0>) {
        chomp($line);
        my @line1 = split(/\t/, $line);
        if (exists $id_to_num_reads{$line1[1]}) {
            $id_to_num_reads{$line1[1]}++;
        } else {
            $id_to_num_reads{$line1[1]} = 1;
        }
    }
    close INPUT0;

    my %id_to_length;
    my $total = 0;
    unlink($out_ids);
    open(OUTPUT, '>>', $out_ids) or die "Error opening $out_ids : $!\n";
    unlink($out_length);
    open(OUTPUT1, '>>', $out_length) or die "Error opening $out_length : $!\n";
    my $seqIO = Bio::SeqIO->new( -file => $contigs, -format => 'fasta' );
    while (my $seqobj = $seqIO->next_seq) {
        my $id = $seqobj->display_id;
        $id_to_length{$id} = $seqobj->length;

        if (not exists $id_to_num_reads{$id}) {
            $id_to_num_reads{$id} = 1;
        }
        $total += $id_to_num_reads{$id};
        print OUTPUT "$id\t$id_to_num_reads{$id}\n";
        print OUTPUT1 "$id\t$id_to_num_reads{$id}\t$id_to_length{$id}\n";
    }
    close OUTPUT1;
    close OUTPUT;
    print "number of all contigs:  ".keys( %id_to_length ).".\n";
    print "total number of reads mapping to contigs = $total.\n\n\n";
}
# This code does not accomplish a purpose other than printing some information
#} else {
#    my %IDs;
#    open(INPUT0, $IDfile) or die "Error opening $IDfile : $!\n";
#    while  ( my $line = <INPUT0>) {
#        chomp($line);
#        my @line1 = split(/\t/, $line);
#        $IDs{$line1[0]} = $line1[1];
#    }
#    close INPUT0;
#
#    my $total = 0;
#    open(INPUT2, $IDfile2) or die "Error opening $IDfile2 : $!\n";
#    while  ( my $line = <INPUT2>) {
#        chomp($line);
#        if (exists $IDs{$line}) {
#            $total += $IDs{$line};
#        } else {
#            $total += 1;
#        }
#    }
#    close INPUT2;
#    open(INPUT2, $IDfile21) or die "Error opening $IDfile21 : $!\n";
#    while  ( my $line = <INPUT2>) {
#        chomp($line);
#        if (exists $IDs{$line}) {
#            $total += $IDs{$line};
#        } else {
#            $total += 1;
#        }
#    }
#    close INPUT2;
#    print "Total number of mapped-reads = $total.\n\n\n";
#}

exit 0;
