#!/usr/bin/env python

"""
Call this script to check up on job progress
Give it a "job list file" (usually located in the "job_output" folder)
Plays "controller" role of MVC
"""

import argparse

# Scripts that should be located in the script directory
from job_logs import create_job_logs
from table_report import get_log_text_report
from dependency_graph import dependency_graph


###################################################################################################
# MAIN
###################################################################################################

def parse_args():
    '''
    Build the argparse.ArgumentParser and parse arguments
    :return: args.{job_file, report, dependency_first, depender_first, success, nosuccess, minimal_detail}
    '''
    parser = argparse.ArgumentParser(description='Display information about a set of jobs')
    parser.add_argument('job_file', help='Path to file containing jobs and associated info')

    # The user must specify at least one of '--report', '--dependendency-first', or '--depender-first'
    parser.add_argument('-r', '--report', help='Display a text report summarizing the jobs and their status', action='store_true')
    parser.add_argument('-dy', '--dependency-first', help="Display an indented text showing job dependencies from first completed to last completed", action='store_true')
    parser.add_argument('-dr', '--depender-first', help="Display an indented text showing job dependencies from last completed to first completed", action='store_true')

    # We can't allow 'success', 'nosucess', and 'minimal_detail' together, since minimal_detail is required for status to be defined
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--success', help='Show successful jobs only', action='store_true')
    group.add_argument('-nos', '--nosuccess', help='Show unsuccessful jobs only i.e. failed or uncompleted jobs', action='store_true')
    group.add_argument('-m', '--minimal-detail', help="Only aggregate minimal information about each job", action='store_true')

    args = parser.parse_args()

    if not (args.report or args.dependency_first or args.depender_first):
        parser.error('One of --report, --dependency-first, or --depender-first is required')

    return args


def main():
    args = parse_args()

    job_logs = create_job_logs(args.job_file, args.success, args.nosuccess, args.minimal_detail)

    if args.report:
        print(get_log_text_report(job_logs, summary=not args.minimal_detail))

    if args.dependency_first:
        print(dependency_graph(job_logs, dependencies_first=True))
    if args.depender_first:
        print(dependency_graph(job_logs, dependencies_first=False))


if __name__ == '__main__':
    main()
