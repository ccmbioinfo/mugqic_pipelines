#!/usr/bin/env perl


#my ($read_files, $db_type, $map_type, $read_type, $strain_type) = @ARGV;
my ($infile, $out_ids, $out_pairs, $out_hits, $read_files, $db_type, $map_type, $read_type, $strain_type) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";


#my $DBfile;
#my $reads;
#my $infile;
#my $outfile1;
#my $outfile2;
#my $outfile3;
#if ($db_type eq 'host') {
#    my $input_dir = "remove_host_reads";
#    my $output_dir = "remove_host_reads";
#
#    $infile = $input_dir."/".$Workpath.$read_files."_host.".$map_type."out";
#    $outfile1 = $output_dir."/".$Workpath.$read_files."_host_".$map_type."_IDs.txt";
#    $outfile2 = $output_dir."/".$Workpath.$read_files."_host_".$map_type."_pairs.txt";
#    $outfile3 = $output_dir."/".$Workpath.$read_files."_host_".$map_type."_hitsID.txt";
#    $infile = $Workpath.$read_files."_host.".$map_type."out";
#    $outfile1 = $Workpath.$read_files."_host_".$map_type."_IDs.txt";
#    $outfile2 = $Workpath.$read_files."_host_".$map_type."_pairs.txt";
#    $outfile3 = $Workpath.$read_files."_host_".$map_type."_hitsID.txt";
#}
#elsif ($db_type eq 'assembly') {
#    my $input_dir = "map_reads";
#    my $output_dir = "map_reads";
#    $infile = $input_dir."/".$read_files."_trinity.".$map_type."out";
#    $outfile1 = $output_dir."/".$Workpath.$read_files."_trinity_".$map_type."_IDs.txt";
#    $outfile2 = $output_dir."/".$Workpath.$read_files."_trinity_".$map_type."_pairs.txt";
#    $outfile3 = $output_dir."/".$Workpath.$read_files."_trinity_".$map_type."_hitsID.txt";
#}
#elsif ($db_type eq 'microgenes') {
#    $DBfile = "microbial_all_cds.fasta";
#
#    $infile = $Workpath.$read_files."_".$read_type."_".$strain_type.".".$map_type."out";
#
#    $outfile1 = $Workpath.$read_files."_".$read_type."_".$strain_type."_".$map_type."_IDs.txt";
#    $outfile2 = $Workpath.$read_files."_".$read_type."_".$strain_type."_".$map_type."_pairs.txt";
#    $outfile3 = $Workpath.$read_files."_".$read_type."_".$strain_type."_".$map_type."_hitsID.txt";
#}
print "Outputfiles: $out_ids\n$out_pairs\n$out_hits\n\n";

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
unlink($out_ids);
open(OUTPUT1, '>>', $out_ids) or die "Error opening $out_ids : $!\n";
unlink($out_pairs);
open(OUTPUT2, '>>', $out_pairs) or die "Error opening $out_pairs : $!\n";
foreach $key (sort keys %IDs) {
    print OUTPUT1 $key."\n";
    print OUTPUT2 $pairs{$key}."\n";
}
close(OUTPUT1);
close(OUTPUT2);

print "number of mapped hits:  ".keys( %hits ).".\n\n\n";
unlink($out_hits);
open(OUTPUT3, '>>', $out_hits) or die "Error opening $out_hits : $!\n";
foreach $key (sort keys %hits) {
    print OUTPUT3 "$key\n";
}
close(OUTPUT3);

exit 0;
