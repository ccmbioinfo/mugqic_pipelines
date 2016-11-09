"""
bismark_merge_reports.py by Michael Li

This script helps merge Bismark\'s alignment reports. This is often required if a sample was sequenced
with multiplexing. In particular, being sequenced by multiple libraries will result in multiple FASTQ
files. By GATK's best practices, it is best process each lane/run separately and merge after alignment is
done. This should be followed by deduplication in non-RRBS datasets.

This method allows maximum parallelization and respects the need to keep read groups seperated.

"""
import argparse
import os.path as path
import re
from shutil import copyfile

# Constant strings containing needed regular expressions to capture data
SEARCH_TOTAL_SEQS = '^Sequence.+? analysed in total.*?\s+([0-9]+)'
SEARCH_UNIQ_HIT = '^Number of .+? alignments with a unique best hit.*?\s+([0-9]+)'
SEARCH_NO_ALIGN = '^Sequence.+? with no alignments under any condition.*?\s+([0-9]+)'
SEARCH_NOT_UNIQ = '^Sequence.+? did not map uniquely.*?\s+([0-9]+)'
SEARCH_DISCARDED = 'Sequence.+? which were discarded because genomic sequence.*?\s+([0-9]+)'
SEARCH_CT_CT = 'CT/GA/CT.*?\s+?([0-9]+)\s+?.*'
SEARCH_GA_CT = 'GA/CT/CT.*?\s+?([0-9]+)\s+?.*'
SEARCH_GA_GA = 'GA/CT/GA.*?\s+?([0-9]+)\s+?.*'
SEARCH_CT_GA = 'CT/GA/GA.*?\s+?([0-9]+)\s+?.*'
SEARCH_DIRECT = 'complementary strands being rejected in total.*?\s+([0-9]+)'
SEARCH_TOTALC = 'Total number of .+? analysed.*?\s+([0-9]+)'
SEARCH_MC_CPG = 'methylated .+? in CpG context.*?\s+([0-9]+)'
SEARCH_MC_CHG = 'methylated .+? in CHG context.*?\s+([0-9]+)'
SEARCH_MC_CHH = 'methylated .+? in CHH context.*?\s+([0-9]+)'
SEARCH_MC_UNK = 'methylated .+? in Unknown context.*?\s+([0-9]+)'
SEARCH_C_CPG = 'unmethylated C\'s in CpG context.*?\s+([0-9]+)'
SEARCH_C_CHG = 'unmethylated C\'s in CHG context.*?\s+([0-9]+)'
SEARCH_C_CHH = 'unmethylated C\'s in CHH context.*?\s+([0-9]+)'
SEARCH_C_UNK = 'unmethylated C\'s in Unknown context.*?\s+([0-9]+)'

# Global counters for our desired criteria
total_seqs, total_c, direction_rejected = (0, 0, 0)
uniq_hit, no_align, not_uniq, discarded = (0, 0, 0, 0)
ct_ct, ga_ct, ga_ga, ct_ga = (0, 0, 0, 0)
mc_cpg, mc_chg, mc_chh, mc_unk = (0, 0, 0, 0)
c_cpg, c_chg, c_chh, c_unk = (0, 0, 0, 0)


def merge_logs(output_report, name, log_reports):
    """
    The main (and only) method that reads all given log reports and adds up various values to produce
    a merged output report. This function has a side-effect of writing an output file at a given
    path.

    :param output_report: The file path where the log report should be written to. The filename need
    not be related to the sample name.
    :type output_report: str
    :param name: The name of the sample to put at the top of the merged report file.
    :type name: str
    :param log_reports: A list of files containing the align reports that will be merged.
    :type log_reports: list(str)
    :return: None
    :rtype: None
    """
    # Check arg values
    if not output_report and name:
        output_report = path.join('.', name + '_aligned_report.txt')
    elif output_report and not name:
        name = path.splitext(path.basename(output_report))
    elif not (output_report and name):
        name = ' '.join([path.splitext(path.basename(log))[0] for log in log_reports])
        output_report = path.join('.', name + '_aligned_report.txt')

    if len(log_reports) == 1:
        copyfile(log_reports[0], output_report)
        return 0

    # Enumerate desired counters
    value_all = [total_seqs, total_c, direction_rejected,
                 uniq_hit, no_align, not_uniq, discarded,
                 ct_ct, ga_ct, ga_ga, ct_ga,
                 mc_cpg, mc_chg, mc_chh, mc_unk,
                 c_cpg, c_chg, c_chh, c_unk]
    search_all = [SEARCH_TOTAL_SEQS, SEARCH_TOTALC, SEARCH_DIRECT,
                  SEARCH_UNIQ_HIT, SEARCH_NO_ALIGN, SEARCH_NOT_UNIQ, SEARCH_DISCARDED,
                  SEARCH_CT_CT, SEARCH_GA_CT, SEARCH_GA_GA, SEARCH_CT_GA,
                  SEARCH_MC_CPG, SEARCH_MC_CHG, SEARCH_MC_CHH, SEARCH_MC_UNK,
                  SEARCH_C_CPG, SEARCH_C_CHG, SEARCH_C_CHH, SEARCH_C_UNK]

    # Read all files and store in memory
    for log in log_reports:
        with open(log) as log_handle:
            file_data = log_handle.readlines()  # For each metric, parse strings to get desired values.
            item = 0
            for val in search_all:
                for each_line in file_data:
                    result = re.search(val, each_line)
                    if result:
                        value_all[item] += int(result.group(1))
                        item += 1
                        break
    # Write out results
    with open(output_report, 'w') as writer:
        writer.write("""\
Merged Bismark report for {sample}: {readsets}

Final Alignment report
======================
Sequences analysed in total:\t{seqs}
Number of alignments with a unique best hit from the different alignments:\t{uniq}
Mapping efficiency:\t{efficient:.1%}
Sequences with no alignments under any condition:\t{no_align}
Sequences did not map uniquely:\t{no_uniq}
Sequences which were discarded because genomic sequence could not be extracted:\t{discarded}

Number of sequences with unique best (first) alignment came from the bowtie output:
CT/GA/CT:\t{ct_ct}\t((converted) top strand)
GA/CT/CT:\t{ga_ct}\t(complementary to (converted) top strand)
GA/CT/GA:\t{ga_ga}\t(complementary to (converted) bottom strand)
CT/GA/GA:\t{ct_ga}\t((converted) bottom strand)

Number of alignments to (merely theoretical) complementary strands being rejected in total:\t{reject}

Final Cytosine Methylation Report
=================================
Total number of C's analysed:\t{total_c}

Total methylated C's in CpG context:\t{mc_cpg}
Total methylated C's in CHG context:\t{mc_chg}
Total methylated C's in CHH context:\t{mc_chh}
Total methylated C's in Unknown context:\t{mc_unk}


Total unmethylated C's in CpG context:\t{c_cpg}
Total unmethylated C's in CHG context:\t{c_chg}
Total unmethylated C's in CHH context:\t{c_chh}
Total unmethylated C's in Unknown context:\t{c_unk}

C methylated in CpG context:\t{rate_cpg:.1%}
C methylated in CHG context:\t{rate_chg:.1%}
C methylated in CHH context:\t{rate_chh:.1%}
C methylated in Unknown context (CN or CHN):\t{rate_unk:.1%}
        """.format(sample=name, readsets=' '.join(log_reports), seqs=value_all[0],
                   uniq=value_all[3], efficient=float(value_all[3]) / float(value_all[0]),
                   no_align=value_all[4], no_uniq=value_all[5], discarded=value_all[6],
                   ct_ct=value_all[7], ga_ct=value_all[8], ga_ga=value_all[9], ct_ga=value_all[10],
                   reject=value_all[2], total_c=value_all[1], mc_cpg=value_all[11],
                   mc_chg=value_all[12], mc_chh=value_all[13],  mc_unk=value_all[14],
                   c_cpg=value_all[15], c_chg=value_all[16], c_chh=value_all[17], c_unk=value_all[18],
                   rate_cpg=float(value_all[11]) / float(value_all[11] + value_all[15]),
                   rate_chg=float(value_all[12]) / float(value_all[12] + value_all[16]),
                   rate_chh=float(value_all[13]) / float(value_all[13] + value_all[17]),
                   rate_unk=float(value_all[14]) / float(value_all[14] + value_all[18])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""This script helps merge Bismark\'s alignment
    reports. This is often required if a sample was sequenced with multiplexing. In particular,
    being sequenced by multiple libraries will result in multiple FASTQ files. By GATK's best practices,
    it is best process each lane/run separately and merge after alignment is done. This should be
    followed by deduplication in non-RRBS datasets.

    This method allows maximum parallelization and respects the need to keep read groups seperated.""",
                                     epilog='By Michael Li, Oct. 2016')
    parser.add_argument('-o', '--output', action='store', default='', type=str,
                        required=False, metavar='file', dest='out_file', help="""The output file
                        path where the merged report should be written to. If not provided, the
                        program will use the current working directory and use the name option
                        to determine the filename. If name is also not provided, the file name
                        is determined by the ID of each log file given.""")
    parser.add_argument('-n', '--name', action='store', default='', type=str,
                        required=False, metavar='sample_id', dest='sample_name', help="""The name of the sample or
                        group to identify the new merged log report. If name is not provided, the file name is
                        determined  by the ID of each log file given. (ex.
                        SRRXXXX_SRRXXXX_...SRRXXXX_aligned_report.txt """)
    parser.add_argument('log_reports', nargs='+', type=str,
                        metavar='Log', help="""A space separated list of Bismark's
                        alignment report files that you want to merge.""")
    args = parser.parse_args()
    merge_logs(args.out_file, args.sample_name, args.log_reports)
