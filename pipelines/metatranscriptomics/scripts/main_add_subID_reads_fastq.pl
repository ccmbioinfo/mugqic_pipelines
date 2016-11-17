#!/usr/bin/env perl

#my ($read_files) = @ARGV;
#my $read_files='cow';

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

for (my $i = 1; $i <= 2; $i++) {
#    $infile = $Workpath.$read_files.$i.".fastq";
#    $outfile = $Workpath.$read_files.$i."_new.fastq";
    $infile = shift;
    $outfile = shift;
    print "Outputfile: $outfile\n";

    unlink($outfile);
    open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
    open(INPUT, $infile) or die "Error opening $infile : $!\n";
    while  ( my $line = <INPUT>) {
        chomp($line);
        my @tmp1 = split(' ', $line);
        my $read1 = join('/', $tmp1[0], $i);

        $line = <INPUT>;
        chomp($line);
        my $seq1 = $line;

        $line = <INPUT>;
        chomp($line);
        my @tmp2 = split(' ', $line);
        my $read2 = join('/', $tmp2[0], $i);

        $line = <INPUT>;
        chomp($line);
        my $qual1 = $line;

        print OUTPUT "$read1\n$seq1\n$read2\n$qual1\n";
    }
    close INPUT;
    close OUTPUT;
}

exit 0;
