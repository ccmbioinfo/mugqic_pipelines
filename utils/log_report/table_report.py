from job_logs import JobLog


# All dates are printed in this format, eg. '2016-01-01T09:02:03'
OUTPUT_DATE_FORMAT = '%Y-%m-%dT%H-%M-%S'

# If a field is undefined, represent it as follows:
UNDEFINED = 'N/A'


###################################################################################################
# SUMMARIZE JOBS
###################################################################################################

def min_max_field(job_logs, field):
    """
    Return name and field of the min and max JobLogs with respect to field
    Return UNDEFINED * 4 if we can't find values
    :param field: JobLog field - must be able to compare values (<)
    :return: min_name, min_field, max_name, max_field strings
    """
    # Sort the jobs by 'field'
    jobs = sorted(filter(lambda log: hasattr(log, field), job_logs), key=lambda log: getattr(log, field))
    if len(jobs) > 0:
        min_job = jobs[0]
        max_job = jobs[-1]
        return min_job.job_name, str(getattr(min_job, field)), max_job.job_name, str(getattr(max_job, field))
    else:
        return UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED


def summarize_status(job_logs):
    """
    return the number of jobs that have status:
             SUCCESS, ACTIVE, INACTIVE, FAILED
    If 'status' is not defined on all jobs, return UNDEFINED's instead
    """
    # If every job has a status, return a breakdown of the jobs
    if len(filter(lambda log: hasattr(log, 'status'), job_logs)) > 0:
        return sum(log.status == 'SUCCESS' for log in job_logs), \
               sum(log.status == 'ACTIVE' for log in job_logs), \
               sum(log.status == 'INACTIVE' for log in job_logs), \
               sum(log.status == 'FAILED' for log in job_logs)
    else:
        return UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED


def get_start_date(job_logs):
    """
    return None if we can't find a start date
    """
    start_dates = [log.start_date for log in job_logs if hasattr(log, 'start_date')]

    if len(start_dates) > 0:
        return min(start_dates)
    return None


def get_end_date(job_logs):
    """
    return None if we can't find a start date
    """
    end_dates = [log.end_date for log in job_logs if hasattr(log, 'end_date')]
    if len(end_dates) > 0:
        return min(end_dates)
    return None


###################################################################################################
# FORMAT LOG
###################################################################################################

def format_field(log, field):
    """
    Lookup the field in log, and format the result for displaying in a table

    :param field: one of the fields of a JobLog
    :return: string representation of log's field
    """
    if hasattr(log, field):
        if field == 'start_date' or field == 'end_date':
            # Always format dates as OUTPUT_DATE_FORMAT
            return str(getattr(log, field).strftime(OUTPUT_DATE_FORMAT))
        else:
            return str(getattr(log, field))
    else:
        # If the field is not defined, display the field as 'N/A'
        return UNDEFINED


def format_log(log):
    """
    :param log: JobLog
    :return: tab-delimited string in order of JobLog.__slots__
    """
    # Eg:
    # 123   456 trimmomatic 789 SUCCESS 0   0   ...
    return '\t'.join(format_field(log, field) for field in JobLog.__slots__)


###################################################################################################
# REPORT
###################################################################################################

def get_log_text_report(job_logs, summary=False):
    """
    Unassigned fields are given the value 'N/A'

    :param job_logs: list of JobLogs
    :param summary: whether to include a summary
    :return: report string
    """
    report = ''

    if not summary:
        start_date = get_start_date(job_logs)
        end_date = get_end_date(job_logs)

        min_mem_name, min_mem, max_mem_name, max_mem = min_max_field(job_logs, 'mem')
        min_walltime_name, min_walltime, max_walltime_name, max_walltime = min_max_field(job_logs, 'walltime')

        num_successful, num_active, num_inactive, num_failed = summarize_status(job_logs)

        if start_date and end_date:
            exec_time = start_date.strftime(OUTPUT_DATE_FORMAT) + ' - ' + end_date.strftime(OUTPUT_DATE_FORMAT) + \
                        ' (' + str(end_date - start_date) + ')'
        else:
            exec_time = UNDEFINED


        report += '''\
# Number of jobs: {num_jobs}
# Number of successful jobs: {successful}
# Number of active jobs: {active}
# Number of inactive jobs: {inactive}
# Number of failed jobs: {failed}

# Execution time: {exec_time}

# Shortest job: {shortest_name} ({shortest_duration})
# Longest job: {longest_name} ({longest_duration})

# Lowest memory job: {min_mem_name} ({min_mem})
# Highest memory job: {max_mem_name} ({max_mem})

'''.format(num_jobs=len(job_logs), successful=num_successful, active=num_active, inactive=num_inactive, failed=num_failed,
           exec_time=exec_time,
           shortest_name=min_walltime_name, shortest_duration=min_walltime,
           longest_name=max_walltime_name, longest_duration=max_walltime,
           min_mem_name=min_mem_name, min_mem=min_mem,
           max_mem_name=max_mem_name, max_mem=max_mem)

    report += '''\
#JOB_ID\tJOB_FULL_ID\tJOB_NAME\tJOB_DEPENDENCIES\tSTATUS\tJOB_EXIT_CODE\tCMD_EXIT_CODE\t\
REAL_TIME\tSTART_DATE\tEND_DATE\tCPU_TIME\tCPU_REAL_TIME_RATIO\tPHYSICAL_MEM\tVIRTUAL_MEM\t\
EXTRA_VIRTUAL_MEM_PCT\tLIMITS\tQUEUE\tUSERNAME\tGROUP\tSESSION\tACCOUNT\tNODES\tPATH
'''

    report += '\n'.join([format_log(log) for log in job_logs])

    return report


