#!/usr/local/bin/perl -w


my ($read_files, $strain_type) = @ARGV;

my $Datapath = "~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath = "";

my %db_strain = ('micro_cds' => "microbial_cds_sub", 'nr' => "nr_all_sub");

my $infile1 = $Workpath.$db_strain{$strain_type}."_IDs_length.txt";    ##[geneID length] or [proteinID_full proteinID length]
my $infile2 = $Workpath.$db_strain{$strain_type}."_IDs_map_taxid_phylum.txt";    ##[geneID/proteinID refID/giID specie taxonID phylum]
my $infile3 = $Workpath.$db_strain{$strain_type}."_IDs_counts.txt";    ##[geneID/proteinID #reads]
my $infile4 = $Workpath."PPI_pairs.txt";    ##[geneID/proteinID b#]
my $outfile = $Workpath.$db_strain{$strain_type}."_table_counts.txt";    ##[geneID/proteinID, length, taxonID, specie, #reads]
print "Outputfile: $outfile\n\n";



### for mapped genes/proteins
my %genes;
open(INPUT1, $infile1) or die "Error opening $infile1 : $!\n";
while  ( my $line = <INPUT1>) {
    chomp($line);
    my @tmp = split(/\t/, $line);
    if ($strain_type eq 'nr') {
        $genes{$tmp[1]} = join("\t", $tmp[1], $tmp[2]);
    } else {
        $genes{$tmp[0]} = $line;
    }
}
close INPUT1;
### genes=[geneID, length]
print "number of genes/proteins:  ".keys( %genes ).".\n";

open(INPUT2, $infile2) or die "Error opening $infile2 : $!\n";
while  ( my $line = <INPUT2>) {
    chomp($line);
    my @tmp = split(/\t/, $line);
    my $tmp2 = $genes{$tmp[0]};
    $genes{$tmp[0]} = join("\t", $tmp2, $tmp[3], $tmp[2], $tmp[4]);
}
close INPUT2;
### genes=[geneID, length, taxonID, specie, phylum]
#print "number of genes/proteins:  " . keys( %genes ) . ".\n";


my $N = 1;
my @N_total = (0)x $N;

open(INPUT3, $infile3) or die "Error opening $infile3 : $!\n";
while  ( my $line = <INPUT3>) {
    chomp($line);
    my @tmp = split(/\t/, $line);
    my $ID = $tmp[0];
    shift @tmp;

    for (my $i = 0; $i < $N; $i++) {
        $N_total[$i] += $tmp[$i];
    }

    my @tmp3 = split(/\t/, $genes{$ID});
    if (scalar(@tmp3) == 2) {
        $genes{$ID} = join("\t", @tmp3, '-', '-', '-', @tmp);
    } else {
        $genes{$ID} = join("\t", @tmp3, @tmp);
    }
}
close INPUT3;
### genes=[geneID, length, taxonID, specie, phylum, #reads]
#print "number of genes/proteins:  " . keys( %genes ) . ".\n";


print "total number of mapped reads:  @N_total.\n\n\n";

open(INPUT4, $infile4) or die "Error opening $infile4 : $!\n";
while  ( my $line = <INPUT4>) {
    chomp($line);
    my @tmp = split(/\t/, $line);
    if (exists $genes{$tmp[0]}) {
        my $tmp4 = $genes{$tmp[0]};
        $genes{$tmp[0]} = join("\t", $tmp4, $tmp[1]);
    }
}
close INPUT4;
### genes=[geneID, length, taxonID, specie, phylum, #reads, PPI]
#print "number of genes/proteins:  " . keys( %genes ) . ".\n";



unlink($outfile);
open(OUTPUT, '>>', $outfile) or die "Error opening $outfile : $!\n";
foreach my $key (sort keys %genes) {
    my @tmp = split(/\t/, $genes{$key});
    if (scalar(@tmp) == (5 + $N)) {
        push @tmp, "-";
    }
    print OUTPUT join("\t", @tmp)."\n";
}
close OUTPUT;

exit 0;