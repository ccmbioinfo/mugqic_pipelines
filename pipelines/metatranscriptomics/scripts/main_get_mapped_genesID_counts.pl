#!/usr/local/bin/perl -w


my ($IDfile, $numfile, $tmpfile, $outfile) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

#my %db_strain = ('micro_cds' => "microbial_cds_sub", 'nr' => "nr_all_sub");
#
#my $IDfile = $Workpath.$db_strain{$strain_type}."_IDs_length.txt";
#my $outfile = $Workpath.$db_strain{$strain_type}."_IDs_counts.txt";
#print "Outputfile: $outfile\n\n";

my $N = 1;
my @Zo = (0)x $N;

my %genes;
my %acc2hit;
if (defined $IDfile) {
    open(INPUT0, $IDfile) or die "Error opening $IDfile : $!\n";
    while  ( my $line = <INPUT0>) {
        chomp($line);

    	next if $line =~ /^\s*@/;

        my @line1 = split(/\t/, $line);

        # Extract accession from fullID (fullID format: gi|gi#|ref/gb|acc#|~)
        # Used Regex instead of split since some files only have accession# instead of fullID
	if(index($IDfile, 'nr_all_sub') != -1) {
#            $accession = $line1[0];
#            $accession =~ s/^gi\|//;
#            $accession =~ s/.*?\|//;
#            $accession =~ s/.*?\|//;
#            $accession =~ s/\|.*?$//;
	     if($line1[0] =~ /([A-Z].*?)(?:\||$)/){
		$accession = $1;
	     } else {
		$accession = $line1[0];
		# Print errorneous line
		print "$line1[0]\n";
	    }
	    $genes{$accession} = join("\t", @Zo);
	}
	else {
	    $genes{$line1[0]} = join("\t", @Zo);
	}

    }
    close INPUT0;
}
print "number of genes/proteins:  ".keys( %genes ).".\n\n\n";

#my $numfile = $Workpath.$read_files."_contigs_IDs.txt";
my %contigs;
open(INPUT1, $numfile) or die "Error opening $numfile : $!\n";
while  ( my $line = <INPUT1>) {
    chomp($line);

    next if $line =~ /^\s*@/;

    my @line1 = split(/\t/, $line);
    $contigs{$line1[0]} = $line1[1];
}
close INPUT1;

#my $tmpfile;
#if ($strain_type eq "micro_cds") {
#    my $IDfile1 = $Workpath.$read_files."_contigs_micro_cds_bwa_pairs.txt";
#    my $IDfile2 = $Workpath.$read_files."_contigs_n_micro_cds_blat_pairs.txt";
#    my $IDfile3 = $Workpath.$read_files."_singletons_micro_cds_bwa_pairs.txt";
#    my $IDfile4 = $Workpath.$read_files."_singletons_n_micro_cds_blat_pairs.txt";
#
#    $tmpfile = $Workpath.$read_files."_bwablat_pairs_tmp.txt";
#    system("cat $IDfile1 $IDfile2 $IDfile3 $IDfile4 > $tmpfile");
#} elsif ($strain_type eq "nr") {
#    my $IDfile1 = $Workpath.$read_files."_contigs_nr_diamond_pairs_sub.txt";
#    my $IDfile2 = $Workpath.$read_files."_singletons_nr_diamond_pairs_sub.txt";
#    $tmpfile = $Workpath.$read_files."_diamond_pairs_tmp.txt";
#    system("cat $IDfile1 $IDfile2 > $tmpfile");
#}

my %hits;
open(INPUT2, $tmpfile) or die "Error opening $tmpfile : $!\n";
while  ( my $line = <INPUT2>) {
    chomp($line);

    next if $line =~ /^\s*@/;

    my @line1 = split(/\t/, $line);

    # Extract accession from fullID (fullID format: gi|gi#|ref/gb|acc#|~)
    # Used Regex instead of split since some files only have accession# instead of fullID
    
    if(index($tmpfile, 'nr_all_sub') != -1) {
#        $accession = $line1[1];
#        $accession =~ s/^gi\|//;
#        $accession =~ s/.*?\|//;
#        $accession =~ s/.*?\|//;
#        $accession =~ s/\|.*?$//;
	if($line1[1] =~ /([A-Z].*?)(?:\||$)/){
	    $accession = $1;
	} else {
	    $accession = $line1[1];
	}

        if (not exists $contigs{$line1[0]}) {
            $hits{$accession} += 1;
        } else {
            $hits{$accession} += $contigs{$line1[0]};
        }
    } else {
        if (not exists $contigs{$line1[0]}) {
            $hits{$line1[1]} += 1;
        } else {
            $hits{$line1[1]} += $contigs{$line1[0]};
        }	
    }

}
close INPUT2;

foreach my $key (keys %hits) {
	if(exists $genes{$key}) {
    	    $genes{$key} = join("\t", $hits{$key});
	}
}

unlink($outfile);
open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
foreach my $key (keys %genes) {
    print OUTPUT $key."\t".$genes{$key}."\n";
}
close OUTPUT;

exit 0;
