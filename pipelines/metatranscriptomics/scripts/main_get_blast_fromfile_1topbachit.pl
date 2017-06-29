#!/usr/local/bin/perl -w

#my ($read_files, $db_type, $map_type, $read_type) = @ARGV;
my ($infile1, $infile2, $outfile1, $outfile2, $outfile3, $BLASTDB) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

#my $readfile0 = $read_files."_".$read_type."_".$db_type;

#my $infile1 = $Workpath.$readfile0."_".$map_type."_IDs.txt";
#my $infile2 = $Workpath.$readfile0."_".$map_type."_pairs.txt";
#my $outfile1 = $Workpath.$readfile0."_".$map_type."_hitsID_sub.txt";
#my $outfile2 = $Workpath.$readfile0."_".$map_type."_pairs_sub.txt";
#my $outfile3 = $Workpath.$readfile0."_".$map_type."_hitsID_bacsub.txt";
#print "Outputfiles: $outfile2\n$outfile1\n\n";
#print "infile1: $infile1\n";
#print "infile2: $infile2\n";
#print "outfile1: $outfile1\n";
#print "outfile2: $outfile2\n";
#print "outfile3: $outfile3\n";
#die;

my %IDs;
my %hits;
my %bachits;
unlink($outfile2);
open(OUTPUT2, '>>', $outfile2) or die "Error opening $outfile2 : $!\n";
open(INPUT1, $infile1) or die "Error opening $infile1 : $!\n";
while( my $line = <INPUT1> ) {
    chomp($line);

    next if $line =~ /^\s*@/;

    $id = (split /\t/, $line)[0];

    my $mycmd = "grep \"".$id."\" ".$infile2;
    my @hits_sub = `$mycmd`;
    my $n = @hits_sub;

    $IDs{$id} = $hits_sub[0];
    my @tmp0 = split(/\t/, $hits_sub[0]);
    $hits{$tmp0[1]} += 1;

    #print "$.\n";

    if ($n > 1) {
        for (my $i = 0; $i < $n; $i++) {
            my @tmp = split(/\t/, $hits_sub[$i]);
            my $mycmd = "blastdbcmd -db $BLASTDB -dbtype 'prot' -entry '".$tmp[1]."' -outfmt %l";
            my $tmp2 = `$mycmd`;
            if (scalar $tmp2) {
                $IDs{$id} = $hits_sub[$i];
                $hits{$tmp[1]} += 1;
                $hits{$tmp0[1]} -= 1;
                $bachits{$tmp[1]} = 1;
                last;
            }
        }
    } elsif ($n == 1) {
        my $mycmd = "blastdbcmd -db $BLASTDB -dbtype 'prot' -entry '".$tmp0[1]."' -outfmt %l";
        my $tmp2 = `$mycmd`;
        if (scalar $tmp2) {
            $bachits{$tmp0[1]} = 1;
        }
    }
    print OUTPUT2 "$IDs{$id}";
}
close INPUT1;
close OUTPUT2;
print "number of reads:  ".keys( %IDs ).".\n";
print "number of bachits:  " . keys( %bachits ) . ".\n";


unlink($outfile1);
open(OUTPUT1, '>>', $outfile1) or die "Error opening $outfile1 : $!\n";
foreach $key (keys %hits) {
    if ($hits{$key} > 0) {
        print OUTPUT1 $key."\n";
    } else {
        delete $hits{$key};
    }
}
close(OUTPUT1);
print "number of hits:  ".keys( %hits ).".\n";

unlink($outfile3);
open(OUTPUT3, '>>', $outfile3) or die "Error opening $outfile3 : $!\n";
foreach $key (keys %bachits) {
    print OUTPUT3 $key."\n";
}
close(OUTPUT3);

exit 0;

