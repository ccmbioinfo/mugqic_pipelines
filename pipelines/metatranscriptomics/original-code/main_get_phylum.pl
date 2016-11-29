#!/usr/local/bin/perl -w
  

my ($read_files, $strain_type)=@ARGV;


my $Datapath="~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath="";


my %db_strain = ('micro_cds'=>"microbial_cds_sub", 'nr'=>"nr_all_sub");	



my $IDfile=$ENV{'BLASTDB'}."/taxid_all_categories.txt";	#[taxonID phylum taxonIDs]
my $infile=$Workpath.$db_strain{$strain_type}."_IDs_map_taxid.txt";	##[geneID/proteinID refID/giID specie taxonID]
my $outfile=$Workpath.$db_strain{$strain_type}."_IDs_map_taxid_phylum.txt";	##[geneID/proteinID refID/giID specie taxonID phylum]
print "Outputfile: $outfile\n\n";


### for mapped genes/proteins
my %genes;
open(INPUT, $infile) or die "Error opening $infile : $!\n";
while  ( my $line = <INPUT>) {
	chomp($line);
	my @tmp=split(/\t/,$line);

	my $mycmd="grep -E \',".$tmp[3].",\|,".$tmp[3]."\$\' ".$IDfile; 
	my @out1=`$mycmd`; 
	chomp(@out1);
	if (scalar(@out1)>0) {
		my @tmp1=split(/\t/,$out1[-1]);
		$genes{$tmp[0]}=join("\t",@tmp,$tmp1[1]);
	}else{
		$mycmd="grep -P \'\t".$tmp[3].",\' ".$IDfile; 
		@out1=`$mycmd`; 
		chomp(@out1);
		if (scalar(@out1)>0) {
			my @tmp2=split(/\t/,$out1[-1]);
			$genes{$tmp[0]}=join("\t",@tmp,$tmp2[1]);
		}else{
			$genes{$tmp[0]}=join("\t",@tmp,'others');
		}
	}
}
close INPUT; 
### ##[geneID/proteinID refID/giID specie taxonID phylum]
print "number of genes/proteins:  " . keys( %genes ) . ".\n\n\n";


unlink($outfile); 
open(OUTPUT,'>>', $outfile) or die "Error opening $outfile : $!\n";
foreach $key (sort keys %genes){
	print OUTPUT $genes{$key}."\n";
}
close OUTPUT;


exit 0;