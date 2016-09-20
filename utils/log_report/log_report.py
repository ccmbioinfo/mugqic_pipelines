#!/usr/bin/env python

"""
Plays the role of both "view" (report) and "controller" (main)
"""

import argparse

from job_logs import create_job_logs
from table_report import get_log_text_report


###################################################################################################
# MAIN
###################################################################################################

def parse_args():
    '''
    Build the argparse.ArgumentParser and parse arguments
    :return: arguments
    '''
    parser = argparse.ArgumentParser(description='Display information about a set of jobs')
    parser.add_argument('job_file', help='Path to file containing jobs and associated info')
    parser.add_argument('-r', '--report', help='Display a text report summarizing the jobs and their status', action='store_true')

    # We can't allow 'success', 'nosucess', and 'minimal_detail' together, since minimal_detail is required for status to be defined
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--success', help='Show successful jobs only', action='store_true')
    group.add_argument('-nos', '--nosuccess', help='Show unsuccessful jobs only i.e. failed or uncompleted jobs', action='store_true')
    group.add_argument('-m', '--minimal_detail', help="Only aggregate minimal information about each job", action='store_true')

    parser.add_argument('-t', '--top_to_bottom', help="Display an indented text showing job dependencies from first completed to last completed", action='store_true')
    parser.add_argument('-b', '--bottom_to_top', help="Display an indented text showing job dependencies from last completed to first completed", action='store_true')

    args = parser.parse_args()

    return args.job_file, args.report, args.success, args.nosuccess, args.minimal_detail


def main():
    job_list_file, report, success_option, no_success_option, minimal_detail = parse_args()

    job_logs = create_job_logs(job_list_file, success_option, no_success_option, minimal_detail=minimal_detail)

    if report:
        print(get_log_text_report(job_logs, summary=minimal_detail))

if __name__ == '__main__':
    main()
