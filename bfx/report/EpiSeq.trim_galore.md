### Read Trimming, Filtering and PostQC Metrics

Trimming and filtering is handled by [Trim Galore!]([@trim_galore]). At this stage of processing, raw sequencing reads are analysed and shortened to ensure high reliability. Trim Galore trims and filters reads to improve data quality.

Trimming reads involves the process of removing low-quality information and sequencing artifacts from the raw data. During the PCR process, a polymerase extends a primer to create the complement of the template sequence. As a sequence grows longer, the polymerase may make more errors in the sequence. This reduces the reliability of a sequence. Therefore, limiting the size of the read will ensure that the low-quality bases are removed from the sequence. Trimming also removes library barcodes from the start of a sequence, which helps determine which sequences belong together. Since these barcodes are artificially place, they must be removed beore the reads are aligned to the genome.

After processing the data, another [FastQC]([@fastqc]) run is done to provide the before and after quality metrics.

#### Reports

| Sample | Readset | Trimming Report | Post-QC Report| Date Submitted | Date Completed |
|--------|---------|-----------------|---------------|----------------|----------------|
$data_table$

#### Download Reports
[Trimming report logs](trim_galore.zip)

[Post-trim QC Metrics](trim_galore_qc.zip) -- only available if `--fastqc` or `--fastqc_args` is specified
