#!/usr/bin/env python
"""
Pipeline Performance
"""
from __future__ import division

import argparse
import re
import sys
from collections import OrderedDict
from datetime import timedelta
from numpy import percentile


class Step(object):
    """
    The step class is a container of tasks that are group together by the step name. This is ideally used in pipelines
    as multiple jobs will share the same root name if they are of the same step.

    This class holds some aggregate data about the group of tasks assigned to the object.
    """

    def __init__(self, name=''):
        """

        :param name: Assigns a name to a new Step object. Empty string otherwise.
        :type name: str
        """
        # Variable Members
        self.tasks = 0  # type: int
        self.success = []  # type: [Task]
        self.failure = []  # type: [Task]

        # Quasi-Constant Members
        self.max_time = None  # type: timedelta
        self.max_mem = 0  # type: int
        self.max_vmem = 0  # type: int
        self.name = name  # type: str

    def __str__(self):
        val = '*' * 100 + '\n'
        val += "Step Name: " + self.name + '\n'
        val += "Total jobs found: " + str(len(self.success) + len(self.failure)) + "\n"
        val += "Success: " + str(len(self.success)) + "\tFailure: " + str(len(self.failure))
        val += "Mem Limit: " + str(self.max_mem) + "\tVmem Limit: " + str(self.max_vmem) + "\n"
        val += "Walltime Limit: " + str(self.max_time) + '\n'
        for task in (self.success + self.failure):
            val += str(task)
        return val

    def get_stat(self, attr):
        """
        Given a function and an attribute, generates a list with the corresponding attribute and
        returns the result from the function against the list.

        :param attr: The attribute to obtain from all member tasks.
        :type attr: str
        :return: The result of the function on the created list.
        :rtype: [any]
        """
        attr_vals = filter(None, [val.get(attr) for val in self.success + self.failure])
        if re.search('^\d\d:\d\d:\d\d$', attr_vals[0]):
            attr_vals = [get_delta(val).total_seconds() for val in attr_vals]
        elif re.search('^\d+\.\d+\s?[A-Z]|[a-z]i?b$', attr_vals[0]):
            attr_vals = [get_bytes(val) for val in attr_vals]
        return attr_vals

    def row_print(self):
        """

        :return: A row of values containing resource data for the given Step object.
        :rtype: str
        """
        row = "{:<35}" + "{:>10}" * 3 + "{:>17}" * 3 + "{:>15}" * 6
        time_list = self.get_stat('Wallclock Duration')
        vmem_list = self.get_stat('vmem Used')
        pmem_list = self.get_stat('Memory Used')
        stat_general = [self.name, self.tasks, len(self.success), len(self.failure)]
        stat_time = [str(self.max_time),
                     str(timedelta(seconds=max(time_list))),
                     str(timedelta(seconds=int(percentile(time_list, 60))) // 1)]
        stat_vmem = [humansize(self.max_vmem),
                     humansize(max(vmem_list)),
                     humansize(int(percentile(vmem_list, 60)))]
        stat_pmem = [humansize(self.max_mem),
                     humansize(max(pmem_list)),
                     humansize(int(percentile(pmem_list, 60)))]
        stat_all = stat_general + stat_time + stat_vmem + stat_pmem
        return row.format(*stat_all)

    def add_task(self, task):
        """
        Add a new task instance or instances to the step object.

        :param task: The task object(s) to add to the list of tasks in the step object.
        :type task: Union(list, Task)
        :rtype: None
        """
        if type(task) is list:
            inst = task.pop()
            self.add_task(task)  # Recursion!
        else:
            inst = task

        # Add metrics if it is a successful run, converting strings to objects as needed.
        if inst.get('Exit Code') == '0':
            self.success.append(inst)
        else:
            self.failure.append(inst)

        # Generate new limit information
        self.set_limits()
        self.tasks += 1

    def remove_task(self, task_name=None, task_id=None):
        """
        Deletes any tasks with the given name and/or id. If only one is given, then only that attribute is used to
        check.

        :param task_name: The name of the task to delete
        :type task_name: str
        :param task_id: The task id number to delete.
        :type task_id: int
        :rtype: None
        """
        if task_name or task_id:
            for i, task in enumerate(self.success + self.failure):
                check_name = task_name == task.get('Job Name')
                check_id = str(task_id) == task.get('Job Id')

                if (task_name and task_id and check_name and check_id) or \
                        (task_name and not task_id and check_name) or \
                        (not task_name and task_id and check_id):
                    if check_name and check_id:
                        if task.get('Exit Code') == 0:
                            del self.success[i]
                        else:
                            del self.failure[i-len(self.success)]
                        self.set_limits()
                        self.tasks -= 1

    def set_limits(self):
        """
        Determines the highest PBS resource limit requested from the list of tasks. Mutates the object.

        :rtype: None
        """
        if not self.success:
            return
        self.max_time = max([get_delta(task.get('Wallclock Limit')) for task in self.success])
        self.max_mem = max([get_bytes(task.get('Memory Limit')) for task in self.success])
        self.max_vmem = max([get_bytes(task.get('vmem Limit')) for task in self.success])


class Task(object):
    """
    The Task class holds runtime statistics of a particular job from PBS/Torque. The various metrics are stored
    as keys to a dictionary and the values from the job are matched to its related key.
    """

    def __init__(self, section=None):
        """

        :param section: A list of strings containing one job entry from showjobs.
        :type section: List[str]
        """
        self.data = OrderedDict()
        for line in section:
            trim = line.strip()
            if trim:
                (attr, value) = line.strip().split(': ')
                self.data[attr.strip()] = value.strip()
        self.name = self.data.setdefault('Job Name', '')
        self.name += '.' + self.data.setdefault('Job Id', '')
        self.name += '.' + self.data.setdefault('Start Time', '')

    def __str__(self):
        val = '-' * 64 + '\n'
        val += self.name + '\n'
        for item in self.data.items():
            val += item[0] + ':\t' + item[1] + '\n'
        val += '-' * 64 + '\n'
        return val

    def get(self, attr):
        """
        A getter function for Task objects. In theory, task.data[key] is also sufficient.

        :param attr: The key/attribute to query with.
        :type attr: str
        :return: The value that is associated with attr in the task's data dictionary.
        :rtype: any
        """
        return self.data.setdefault(attr)

    def set_value(self, attr, val):
        """
        A setter function to go with a getter function. In theory, task.data[key] = val will also work.

        :param attr: The attribute which will be set/added to the dictionary
        :type attr: str
        :param val: The value to associate with the atrribute
        :type val: any
        :rtype: None
        """
        self.data[attr] = val


class StatsManager(object):
    """
    v
    """
    template = """
Performance and Resource Statistics for Pipeline {pipeline_name}

{table}


    """
    metrics = ["Step Name", "Total Runs", "Success", "Failed",
               "Walltime Limit", "Max Walltime", "Avg Walltime",
               "V. Mem Limit", "Max V. Mem", "Avg V. Mem",
               "Memory Limit", "Max Memory", "Avg Memory"]
    header = ("{:<35}" + "{:>10}" * 3 + "{:>17}" * 3 + "{:>15}" * 6).format(*metrics)

    def __init__(self, step_list, name=None, pipeline=None):
        """

        :param step_list: A mapping between the step name and the step objects.
        :type step_list: Dict[str, Step]
        :param name: A string that will be used to name the StatsManager object.
        :type name: str
        :param pipeline: A list of strings containing the proper name of each step in a pipeline.
        :type pipeline: List[str]
        """
        self.name = name
        self.steps = step_list
        if pipeline:
            with open(pipeline) as f:
                self.pipeline = f.readlines()
        else:
            self.pipeline = [step.name for step in step_list]
        if pipeline and step_list:
            self.filter_steps()

    def filter_steps(self):
        """
        Helps discard all irrelevant jobs that is not listed in self.pipeline.

        :rtype: None
        """
        if not self.pipeline:
            raise AttributeError('Missing object member: pipeline')
        new_step = dict()
        for step in self.pipeline:
            step = step.strip()
            new_step[step] = self.steps[step]
        if new_step:
            self.steps = new_step
        else:
            raise ValueError('No steps found that matches the steps in your pipeline')

    def print_table(self):
        """
        Prints aggregate data into a table for easy reading.

        :rtype: None
        """
        out_str = self.header + '\n'  # type: str
        for step in self.steps:
            out_str += self.steps[step].row_print() + '\n'
        print (self.template.format(pipeline_name=self.name, table=out_str))


def parse_pbs_job_data(in_file):
    """
    Either from stdin or a file, parse data from showjobs and produce a set of Step objects or analysis.

    :param in_file: A log file containing the data from showjobs on the PBS system.
    :type in_file: str
    :return: The dictionary object containing the list of steps that are identified from the input log.
    :rtype: Union[None, Dict[str, Step]
    """
    delim = '--------------------------------------------------------------------------------'
    entries = []  # type: [Task]
    entry = []  # type: list(str)

    if in_file:  # File log is given
        with open(in_file) as f:
            for line in f:
                if line.startswith(delim):
                    if entry:  # If we hit a delimiter, create object
                        entries.append(Task(entry))
                        entry = []  # reset buffer
                        continue  # Skip delimiter line
                elif line:
                    entry.append(line)  # Don't add blank lines
            if entry:
                entries.append(Task(entry))
    else:  # Read from stdin
        f = sys.stdin
        for line in f:
            if line.startswith(delim):
                if entry:  # If we hit a delimiter, create object
                    entries.append(Task(entry))
                    entry = []  # reset buffer
                    continue  # Skip delimiter line
            elif line:
                entry.append(line)  # Don't add blank lines
        if entry:
            entries.append(Task(entry))
    steps_handle = generate_steps(entries)
    for item in steps_handle.values():  # type: Step
        item.set_limits()
    return steps_handle


def generate_steps(jobs):
    """
    Given a list of jobs, creates several steps by grouping jobs by name. The steps are stored as a dictionary.

    :param jobs: A list of jobs parsed from the input log data.
    :type jobs: List[Task]
    :return: A mapping from a step name to the associated Step object.
    :rtype: Dict[str, Step]
    """
    step_map = dict()

    for job in jobs:  # type: Task
        if len(job.data) <= 3:
            continue  # This job doesn't have enough data for anything
        name = job.get('Job Name').split('.')[0]
        if name not in step_map:
            step_map[name] = Step(name=name)  # Set limits here.
        step_map[name].add_task(job)

    return step_map


def get_bytes(size_string):
    """
    Converts any data size from string to integer as bytes.

    :param size_string: The input string containing the data size.
    :type size_string: str
    :return: The number of bytes equivalent to size_string
    :rtype: int
    """
    try:
        size_string = size_string.lower().replace(',', '')
        size = float(re.search('^(\d+(\.\d+)?)\s?[a-z]?b$', size_string).groups()[0])
        suffix = re.search('^\d+(\.\d+)?\s?([kmgtp])?b$', size_string).groups()[1]
    except AttributeError:
        raise
    symbols = ('b', 'k', 'm', 'g', 't', 'p', 'e', 'z', 'y')
    prefix = {symbols[0]: 1}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i + 1) * 10
    return int(size * prefix[suffix])


# Attributed to: nneonneo @ StackOverflow
# (http://stackoverflow.com/questions/14996453/python-libraries-to-calculate-human-readable-filesize-from-bytes)
def humansize(nbytes):
    """
    Converts a number of bytes from an integer to a human readable value.

    :param nbytes: The input string containing the data size.
    :type nbytes: int
    :return: The number of bytes equivalent to nbytes
    :rtype: str
    """
    if nbytes == 0:
        return '0 B'
    i = 0
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    while nbytes >= 1024 and i < len(suffixes) - 1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('.')
    return '%s %s' % (f, suffixes[i])


def get_delta(time_str):
    """
    Converts a given amount of time as a string and produces the equivalent datetime.timedelta object.

    :param time_str: The amount of time in the format dd:hh:mm:ss
    :type time_str: str
    :return: A datetime.timedelta object
    :rtype: datetime.timedelta
    """
    delta_obj = timedelta()
    if time_str:
        try:
            components = filter(None, time_str.split(':'))
            components = [int(item) for item in components]
            if len(components) == 4:
                delta_obj = timedelta(days=components[0], hours=components[1],
                                      minutes=components[2], seconds=components[3])
            elif len(components) == 3:
                delta_obj = timedelta(hours=components[0], minutes=components[1], seconds=components[2])
            elif len(components) == 2:
                delta_obj = timedelta(minutes=components[0], seconds=components[1])
            elif len(components) == 1:
                delta_obj = timedelta(seconds=components[0])
            else:
                raise ValueError
        except ValueError:
            raise ValueError('Unknown time string: ' + time_str)
    return delta_obj


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""This script reads in data from showjobs and create aggregate
    performance statistics. It accepts from stdin or files, as needed.""",
                                     epilog='By Michael Li, Nov. 2016')
    parser.add_argument('-i', '--input', type=str, metavar='log_file', action='store', dest='in_file',
                        help="""Output from showjobs, stored in a file.""")
    parser.add_argument('-s', '--steps', type=str, metavar='list_steps', action='store', dest='pipe_steps',
                        help="""Names of steps of a pipeline in a file, with one step per line.""")
    args = parser.parse_args()

    steps = parse_pbs_job_data(args.in_file)
    StatsManager(steps, '', args.pipe_steps).print_table()
    exit(0)
