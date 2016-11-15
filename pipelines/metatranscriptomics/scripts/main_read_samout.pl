#!/usr/local/bin/perl -w


my ($read_files, $db_type, $map_type, $read_type, $strain_type) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

my $DBfile;
my $reads;
my $infile;
my $outfile1;
my $outfile2;
my $outfile3;
if ($db_type eq 'host') {
    $infile = $Workpath.$read_files."_host.".$map_type."out";
    $outfile1 = $Workpath.$read_files."_host_".$map_type."_IDs.txt";
    $outfile2 = $Workpath.$read_files."_host_".$map_type."_pairs.txt";
    $outfile3 = $Workpath.$read_files."_host_".$map_type."_hitsID.txt";
}
elsif ($db_type eq 'assembly') {
    $infile = $Workpath.$read_files."_trinity.".$map_type."out";
    $outfile1 = $Workpath.$read_files."_trinity_".$map_type."_IDs.txt";
    $outfile2 = $Workpath.$read_files."_trinity_".$map_type."_pairs.txt";
    $outfile3 = $Workpath.$read_files."_trinity_".$map_type."_hitsID.txt";
}
elsif ($db_type eq 'microgenes') {
    $DBfile = "microbial_all_cds.fasta";

    $infile = $Workpath.$read_files."_".$read_type."_".$strain_type.".".$map_type."out";

    $outfile1 = $Workpath.$read_files."_".$read_type."_".$strain_type."_".$map_type."_IDs.txt";
    $outfile2 = $Workpath.$read_files."_".$read_type."_".$strain_type."_".$map_type."_pairs.txt";
    $outfile3 = $Workpath.$read_files."_".$read_type."_".$strain_type."_".$map_type."_hitsID.txt";
}
print "Outputfiles: $outfile1\n$outfile2\n$outfile3\n\n";

my %IDs;
my %pairs;
my %hits;
my $line_past = 0;
my $line_now;
open(INPUT, $infile) or die "Error opening $infile : $!\n";
while(my $line = <INPUT>) {
    chomp($line);
    my @line1;
    @line1 = split(/\t/, $line);

    if ($read_type eq 'pairs') {
        $line_now = $line1[0];
        if ($line_now ne $line_past) {
            $IDs{$line_now} = 1;
            $pairs{$line_now} = join("\t", $line1[0], $line1[2]);
            $hits{$line1[2]} = 1;
            $line_past = $line_now;
        }
    } else {
        $IDs{$line1[0]} = 1;
        $pairs{$line1[0]} = join("\t", $line1[0], $line1[2]);
        $hits{$line1[2]} = 1;
    }
}
close(INPUT);

print "number of mapped reads:  ".keys( %IDs ).".\n";
print "number of mapped pairs:  ".keys( %pairs ).".\n";
unlink($outfile1);
open(OUTPUT1, '>>', $outfile1) or die "Error opening $outfile1 : $!\n";
unlink($outfile2);
open(OUTPUT2, '>>', $outfile2) or die "Error opening $outfile2 : $!\n";
foreach $key (sort keys %IDs) {
    print OUTPUT1 $key."\n";
    print OUTPUT2 $pairs{$key}."\n";
}
close(OUTPUT1);
close(OUTPUT2);

print "number of mapped hits:  ".keys( %hits ).".\n\n\n";
unlink($outfile3);
open(OUTPUT3, '>>', $outfile3) or die "Error opening $outfile3 : $!\n";
foreach $key (sort keys %hits) {
    print OUTPUT3 "$key\n";
}
close(OUTPUT3);

exit 0;
