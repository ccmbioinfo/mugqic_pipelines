#!/usr/bin/env perl
use List::Util qw[min max];

my ($read_files, $read_type, $cutoff_type, $cutoff1, $cutoff2) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

if ($cutoff_type == 0) {
    $cutoff1 = 0;
    $cutoff2 = 0;
}

my $II = 1;
if ($read_type eq 'pairs') {
    $II = 2;
}

for (my $i = 1; $i <= $II; $i++) {
    my $lengthfile = "remove_rrna/cow" . $i . "_IDs_length.txt";
#    my $lengthfile = "remove_rrna/cow" . $i . "_qual_all_unique_IDs_length.txt";
    my $infile = "remove_rrna/cow" . $i . "_rRNA.infernalout";

    my $outfile1 = "remove_rrna/cow" . $i . "_rRNA_infernal_IDs.txt";
    my $outfile2 = "remove_rrna/cow" . $i . "_rRNA_infernal_pairs.txt";
    my $outfile3 = "remove_rrna/cow" . $i . "_rRNA_infernal_hitsID.txt";
    print "Outputfile: $outfile1\n\n";

    my %reads;
    my %reads2;
    open(INPUT0, $lengthfile) or die "Error opening $lengthfile : $!\n";
    while( my $line = <INPUT0> ) {
        chomp($line);
        my @values = split(/\t/, $line);
        $reads{$values[0]} = $values[2];
        $reads2{$values[0]} = $values[1];
    }
    close INPUT0;
    print "number of all unique reads:  ".keys( %reads ).".\n";

    my %IDs;
    my %pairs;
    my %hits;
    my $max_evalue;
    open(INPUT, $infile) or die "Error opening $infile : $!\n";
    while  (my $line0 = <INPUT>) {
        chomp($line0);
        if ($line0 !~ m/^\#/) {
            my $tmp1 = $line0;
            $tmp1 =~ s/ .*?$//;

            $line0 =~ s/^.*? RF//;

            my $tmp2 = substr($line0, 8);
            $tmp2 =~ s/ .*?$//;
            my @tmp222 = split(';', $tmp2);
            $tmp2 = $tmp222[0];

            if ($line0 =~ m/ -          cm /g) {
                $line0 =~ s/^.*?( -          cm )//;
            } elsif ($line0 =~ m/ -         hmm /g) {
                $line0 =~ s/^.*?( -         hmm )//;
            }

            my $tmp3 = substr($line0, 0, 26);
            $tmp3 =~ s/.*? //g;

            my $tmp4 = substr($line0, 0, 35);
            $tmp4 =~ s/.*? //g;

            my $tmp5 = substr($line0, 0, 71);
            $tmp5 =~ s/.*? //g;

            my $tmp6 = substr($line0, 0, 81);
            $tmp6 =~ s/.*? //g;

            my @line = ($tmp1, $tmp2, $tmp3, $tmp4, $tmp5, $tmp6);

            if (not exists $IDs{$line[1]}) {
                $max_evalue = $line[5];

                my $mydbID = $line[1];

                my $hitgiID = $line[0];

                my $pident;
                if (not exists $reads{$mydbID}) {
                    $reads{$mydbID} = 100;
                    $reads2{$mydbID} = 2;
                }
                $pident = 100 * abs($line[3] - $line[2] + 1) / $reads{$mydbID};
                $pident = sprintf("%.2f", $pident);

                my $score = $line[4];

                if ($cutoff_type == 0) {
                    $IDs{$mydbID} = $reads2{$mydbID};
                    $pairs{$mydbID} = join("\t", $mydbID, $hitgiID, $pident, $score, $max_evalue);
                    if (exists $hits{$hitgiID}) {
                        $hits{$hitgiID} += $IDs{$mydbID};
                    } else {
                        $hits{$hitgiID} = $IDs{$mydbID};
                    }
                } else {
                    if (($max_evalue <= $cutoff1) and ($pident >= $cutoff2)) {
                        $IDs{$mydbID} = $reads2{$mydbID};
                        $pairs{$mydbID} = join("\t", $mydbID, $hitgiID, $pident, $score, $max_evalue);
                        if (exists $hits{$hitgiID}) {
                            $hits{$hitgiID} += $IDs{$mydbID};
                        } else {
                            $hits{$hitgiID} = $IDs{$mydbID};
                        }
                    }
                }
            }
        }
    }
    close INPUT;

    print "number of mapped rRNA unique reads:  ".keys( %IDs ).".\n";
    my $N_mapped = 0;
    unlink($outfile1);
    open(OUTPUT1, '>>', $outfile1) or die "Error opening $outfile1 : $!\n";
    unlink($outfile2);
    open(OUTPUT2, '>>', $outfile2) or die "Error opening $outfile2 : $!\n";
    foreach my $key (keys %IDs) {
        $N_mapped += $IDs{$key};
        print OUTPUT1 "$key\t$IDs{$key}\n";
        print OUTPUT2 "$pairs{$key}\n";
    }
    close OUTPUT1;
    close OUTPUT2;
    print "number of mapped rRNA reads:  ".$N_mapped.".\n";

    print "number of mapped hits:  ".keys( %hits ).".\n\n\n";
    my $N = 0;
    unlink($outfile3);
    open(OUTPUT3, '>>', $outfile3) or die "Error opening $outfile3 : $!\n";
    foreach my $key (keys %hits) {
        $N = $N + $hits{$key};
        print OUTPUT3 "$key\t$hits{$key}\n";
    }
    close(OUTPUT3);
    #print "number of mapped rRNA reads:  " . $N . ".\n";
}

if ($II == 2) {
    my $output_dir = "remove_rrna/";
    my $infile11 = $output_dir . "cow1_rRNA_infernal_IDs.txt";
    my $infile12 = $output_dir . "cow1_rRNA_infernal_pairs.txt";
    my $infile13 = $output_dir . "cow1_rRNA_infernal_hitsID.txt";

    my $infile21 = $output_dir . "cow2_rRNA_infernal_IDs.txt";
    my $infile22 = $output_dir . "cow2_rRNA_infernal_pairs.txt";
    my $infile23 = $output_dir . "cow2_rRNA_infernal_hitsID.txt";

    my $outfile1 = $output_dir . "cow_rRNA_infernal_IDs.txt";
    my $outfile2 = $output_dir . "cow_rRNA_infernal_pairs.txt";
    my $outfile3 = $output_dir . "cow_rRNA_infernal_hitsID.txt";
    print "Outputfiles: $outfile1\n$outfile2\n$outfile3\n\n";

    system("cat $infile12 $infile22 > $outfile2");

    my %IDs;
    my %IDs0;
    open(INPUT1, $infile11) or die "Error opening $infile11 : $!\n";
    while  (my $line1 = <INPUT1>) {
        chomp($line1);
        my @values = split(/\t/, $line1);
        my @line11 = split('\/', $values[0]); ## remove the suspended "/1" or "/2" of each ID
        $IDs{$values[0]} = $values[1];
        $IDs0{$line11[0]} = $values[1];
    }
    close INPUT1;
    #print "number of unique rRAN reads1:  " . keys( %IDs0 ) . ".\n";

    open(INPUT2, $infile21) or die "Error opening $infile21 : $!\n";
    while  (my $line1 = <INPUT2>) {
        chomp($line1);
        my @values = split(/\t/, $line1);
        my @line11 = split('\/', $values[0]); ## remove the suspended "/1" or "/2" of each ID
        $IDs{$values[0]} = $values[1];
        if (exists $IDs0{$line11[0]}) {
            my $tmpp = $IDs0{$line11[0]};
            $IDs0{$line11[0]} = max($tmpp, $values[1]);
        } else {
            $IDs0{$line11[0]} = $values[1];
        }
    }
    close INPUT2;
    #print "number of unique rRAN reads1+2:  " . keys( %IDs0 ) . ".\n";


    my $N_mapped = 0;
    unlink($outfile1);
    open(OUTPUT1, '>>', $outfile1) or die "Error opening $outfile1 : $!\n";
    foreach  my $key (keys %IDs0) {
        $N_mapped += $IDs0{$key};
        print OUTPUT1 "$key\t$IDs0{$key}\n";
    }
    close(OUTPUT1);
    print "number of mapped rRAN reads:  ".$N_mapped.".\n";

    my %hits;
    open(INPUT3, $infile13) or die "Error opening $infile13 : $!\n";
    while  (my $line1 = <INPUT3>) {
        chomp($line1);
        my @line11 = split('\t', $line1);
        $hits{$line11[0]} = $line11[1];
    }
    close INPUT3;
    #print "number of rRAN hits1:  " . keys( %hits ) . ".\n";
    open(INPUT4, $infile23) or die "Error opening $infile23 : $!\n";
    while  (my $line1 = <INPUT4>) {
        chomp($line1);
        my @line11 = split('\t', $line1);
        if (exists $hits{$line11[0]}) {
            $hits{$line11[0]} += $line11[1];
        } else {
            $hits{$line11[0]} = $line11[1];
        }
    }
    close INPUT4;
    print "number of all rRAN hits:  ".keys( %hits ).".\n";

    my $N = 0;
    unlink($outfile3);
    open(OUTPUT3, '>>', $outfile3) or die "Error opening $outfile3 : $!\n";
    foreach my $key (keys %hits) {
        $N = $N + $hits{$key};
        print OUTPUT3 "$key\t$hits{$key}\n";
    }
    close(OUTPUT3);
    #print "number of all mapped rRAN reads:  " . $N . ".\n";
}

exit 0;

