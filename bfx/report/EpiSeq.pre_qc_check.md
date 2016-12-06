### FastQC for Bisulfite Sequencing

[FastQC] runs on sequencing data and measure various statistics to inform the user about the quality of the raw data and what precautions to consider going forward. Information about library contamination and poor base quality can be used to guide the trimming and filtering processes [@fastqc].

#### HTML Reports

| Sample | Readset | File 1 | File 2 | Date Submitted | Date Completed |
|:------:|:-------:|:-------|:-------|:---------------|:---------------|
$data_table$

#### Note:
Since reads have been altered by the bisulfite reagent, the QC metrics needs to be interpreted differently. In particular, the percent composition of each nucleotide will be altered because all unmethylated cytosines will be converted to uracil. This will be interpreted as thymine by the sequencer. Other warnings like repetitive sequences can sometimes be ignored as well. For example, reads from a RRBS dataset will always start with a portion of the recognition site for _Mspl_.

#### For more information:

[Guide for WGBS datasets](http://www.epigenesys.eu/images/stories/protocols/pdf/20120720103700_p57.pdf)

[Guide for RRBS datasets](http://www.bioinformatics.babraham.ac.uk/projects/bismark/RRBS_Guide.pdf)
