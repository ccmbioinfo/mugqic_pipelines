#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;
  

my ($read_files, $db_type, $map_type, $read_type, $strain_type)=@ARGV;


my $Datapath="~/CourseData/metagenomics/metatranscriptomics/";
my $Workpath="";



my $II=1;
if (($read_type eq 'pairs') or ($read_type eq 'singletons')){
	$II=2;
}


for (my $i=1; $i<=$II; $i++) {
	my $IDfile;
	my $infile;
	my $outfile1;
	my $outfile2;
	if ($db_type eq 'rRNA') {
		$IDfile=$Workpath.$read_files."_rRNA_infernal_IDs.txt";

		$infile=$Workpath.$read_files.$i."_qual_all_unique.fastq";

		$outfile1=$Workpath.$read_files.$i."_qual_unique_rRNA.fastq";
		$outfile2=$Workpath.$read_files.$i."_qual_unique_n_rRNA.fastq";
	}
	elsif ($db_type eq 'host') {
		$IDfile = $Workpath.$read_files."_host_bwa_IDs.txt";
		$infile=$Workpath.$read_files.$i."_qual_unique_n_rRNA.fastq";

		$outfile1=$Workpath.$read_files.$i."_qual_unique_n_rRNA_host.fastq";
		$outfile2=$Workpath.$read_files.$i."_qual_unique_n_rRNA_n_host.fastq";
	}
	elsif ($db_type eq 'assembly') {
		$IDfile=$Workpath.$read_files."_trinity_bwa_IDs.txt";
		$infile=$Workpath.$read_files.$i."_mRNA.fastq";

		$outfile1=$Workpath.$read_files.$i."_mRNA_mappedreads.fastq";
		$outfile2=$Workpath.$read_files.$i."_singletons.fastq";
	}
	elsif ($db_type eq 'microgenes') {
		$IDfile=$Workpath.$read_files."_".$read_type."_".$strain_type."_".$map_type."_IDs.txt";
		if ($read_type eq 'contigs') {
			$infile = $Workpath.$read_files."_".$read_type.".fasta";
			$outfile1 = $Workpath.$read_files."_".$read_type."_".$strain_type.".fasta";
			$outfile2 = $Workpath.$read_files."_".$read_type."_n_".$strain_type.".fasta";
		}else{
			$infile = $Workpath.$read_files.$i."_".$read_type.".fastq";
			$outfile1 = $Workpath.$read_files.$i."_".$read_type."_".$strain_type.".fasta";
			$outfile2 = $Workpath.$read_files.$i."_".$read_type."_n_".$strain_type.".fasta";
		}
	}
	elsif ($db_type eq 'microgenes_blat') {
		$IDfile=$Workpath.$read_files."_".$read_type."_n_".$strain_type."_".$map_type."_IDs.txt";
		if ($read_type eq 'contigs') {
			$infile = $Workpath.$read_files."_".$read_type."_n_".$strain_type.".fasta";
			$outfile1 = $Workpath.$read_files."_".$read_type."_n_".$strain_type."_blat.fasta";
			$outfile2 = $Workpath.$read_files."_".$read_type."_n_".$strain_type."_rest.fasta";
		}else{
			my $tmp=$Workpath.$read_files.$i."_singletons";
			$infile = $tmp."_n_".$strain_type.".fasta";
			$outfile1 = $tmp."_n_".$strain_type."_blat.fasta";
			$outfile2 = $tmp."_n_".$strain_type."_rest.fasta";
		}
	}
	elsif ($db_type eq 'nr') {
		$IDfile=$Workpath.$read_files."_nr_IDs.txt";
		if ($read_type eq 'contigs') {
			$infile = $Workpath.$read_files."_".$read_type."_n_".$strain_type."_rest.fasta";
			$outfile1 = $Workpath.$read_files."_".$read_type."_nr.fasta";
			$outfile2 = $Workpath.$read_files."_".$read_type."_n_nr.fasta";
		}else{
			my $tmp=$Workpath.$read_files.$i."_singletons";
			$infile = $tmp."_n_".$strain_type."_rest.fasta";
			$outfile1 = $tmp."_nr.fasta";
			$outfile2 = $tmp."_n_nr.fasta";
		}
	}
	print "Outputfiles: $outfile1\n$outfile2\n\n";


	my %IDs;
	open(INPUT0, $IDfile) or die "Error opening $IDfile : $!\n";
	while  (my $line = <INPUT0>) {
		chomp($line);
		if ($db_type eq 'rRNA'){
			my @line1 = split(/\t/, $line);
			$IDs{$line1[0]} = $line1[1];
		}else{
			$IDs{$line} = 1;
		}
	}
	close INPUT0; 
	print "number of reads for ".$db_type.":  " . keys( %IDs ) . ".\n\n\n";



	my ($root,$suff_in)=$infile=~m/(.*)\.(.*)/;
	my ($root,$suff_out)=$outfile1=~m/(.*)\.(.*)/;

	my $seq_out1 = Bio::SeqIO->new(-file=>">$outfile1",-format=>$suff_out);
	my $seq_out2 = Bio::SeqIO->new(-file=>">$outfile2",-format=>$suff_out);
	my $seq_in=Bio::SeqIO->new(-file=>$infile, -format=>$suff_in);
	while (my $seqobj=$seq_in->next_seq) {
		my $head=$seqobj->display_id;
		my $seq=$seqobj->seq;

		my @line1 = split('\/', $head);

		if (exists $IDs{$line1[0]}) {
			$seq_out1->write_seq($seqobj);
		}else{
			$seq_out2->write_seq($seqobj);
		}
	}
}

exit 0;
