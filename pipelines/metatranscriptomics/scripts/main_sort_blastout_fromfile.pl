#!/usr/local/bin/perl -w
use List::Util qw[min max];

my ($read_files, $db_type, $map_type, $read_type, $cutoff_num) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

my $infile = $Workpath.$read_files."_".$read_type."_".$db_type.".".$map_type."out";
if ($read_type eq 'singletons') {
    my $read1 = $Workpath.$read_files."1_singletons_".$db_type.".".$map_type."out";
    my $read2 = $Workpath.$read_files."2_singletons_".$db_type.".".$map_type."out";
    system("cat $read1 $read2 > $infile");
}
my $outfile = $Workpath.$read_files."_".$read_type."_".$db_type."_sorted.".$map_type."out";
print "Outputfile: $outfile\n\n";

my %IDs;
unlink($outfile);
open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
open(INPUT, $infile) or die "Error opening $infile : $!\n";
while  (my $line1 = <INPUT>) {
    chomp($line1);
    my @line = split(/\t/, $line1);
    my @head = split('\/', $line[0]);

    if (not exists $IDs{$head[0]}) {
        $IDs{$head[0]} = 1;

        my $cmd = 'grep "'.$head[0].'" '.$infile;
        my @hit_pairs = qx($cmd);
        chomp(@hit_pairs);
        my $N1 = @hit_pairs;

        if ($N1 == 1) {
            print OUTPUT "$hit_pairs[0]\n";
        } else {
            my %pairs;
            my %scores;
            for (my $k = 0; $k < $N1; $k++) {
                my @pair = split(/\t/, $hit_pairs[$k]);
                $pairs{$k} = $hit_pairs[$k];
                $scores{$k} = $pair[11];
            }

            my @hits = sort { $scores{$b} <=> $scores{$a} } keys(%scores);
            for (my $k = 0; $k < min($N1, $cutoff_num); $k++) {
                print OUTPUT "$pairs{$hits[$k]}\n";
            }
        }
    }
}
close INPUT;
close OUTPUT;

exit 0;

