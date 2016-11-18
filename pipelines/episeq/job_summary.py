#!/bin/usr/env python
"""
Pipeline Performance
"""
import argparse
import datetime
import string
import re


class Step(object):
    """
    a
    """
    name = ''
    tasks = []
    success = 0
    failure = 0
    walltimes = []
    cputimes = []
    vmem = []
    pmem = []

    def __init__(self, name, max_time=datetime.timedelta(), max_vmem=0, max_pmem=0):
        self.name = name
        self.max_time = max_time
        self.max_vmem = max_vmem
        self.max_pmem = max_pmem

    def add_task(self, task):
        """

        :param task:
        :type task: Union(list, Task)
        :return:
        :rtype:
        """
        if type(task) is list:
            inst = task.pop()
            self.add_task(task)  # Recursion!
        else:
            inst = task

        if inst.get('Exit Code') == '0':
            self.success += 1
            self.walltimes.append(inst.get('Wallclock Duration'))
            self.cputimes.append(inst.get('CPUTime'))
            self.vmem.append(inst.get('vmem Used'))  # TODO normalize memory values
            self.pmem.append(inst.get('Memory Used'))
        else:
            self.failure += 1
        pass


class Task(object):
    """
    b
    """
    def __init__(self, section=None):
        self.data = dict()
        for attr, value in [line.split(':') for line in section]:
            self.data[attr.strip()] = value.strip()

    def get(self, attr):
        """

        :param attr:
        :type attr: str
        :return:
        :rtype: None
        """
        if attr not in self.data:
            raise KeyError('Attribute ' + attr + ' is not in this task. Use one of ' + ', '.join(self.data.keys()))
        return self.data.setdefault(attr)

    def set_value(self, attr, val):
        """

        :param attr:
        :type attr: str
        :param val:
        :type val: any
        :return:
        :rtype: None
        """
        self.data[attr] = val


class StatsManager(object):
    """
    v
    """
    def __init__(self):
        pass


def parse_pbs_job_data(in_file):
    """

    :param in_file:
    :type in_file: str
    :return:
    :rtype:
    """
    delim = '--------------------------------------------------------------------------------'
    entries = []  # type: [Task]
    entry = []  # type: list(str)

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
    steps_handle = generate_steps(entries)


def generate_steps(jobs):
    """

    :param jobs:
    :type jobs: [Task]
    :return:
    :rtype: dict
    """
    step_map = dict()

    for job in jobs:  # type: Task
        name = job.get('Job Name').split[0]
        if name not in step_map:
            step_map[name] = Step(name=name)  # Set limits here.
        step_map[name].add_job(job)

    return step_map


def get_bytes(size_string):
    try:
        size_string = size_string.lower().replace(',', '')
        size = re.search('^(\d+)[a-z]i?b$', size_string).groups()[0]
        suffix = re.search('^\d+([kmgtp])i?b$', size_string).groups()[0]
    except AttributeError:
        raise ValueError("Invalid Input")
    shft = suffix.translate(string.maketrans('kmgtp', '12345')) + '0'
    return int(size) << int(shft)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""This script reads in data from showjobs and create aggregate
    performance statistics.""",
                                     epilog='By Michael Li, Nov. 2016')
    parser.add_argument('in_file', type=str, metavar='Log',
                        help="""Output from showjobs.""")
    args = parser.parse_args()

    parse_pbs_job_data(args.in_file)
