#!/usr/bin/env python
"""
Pipeline Performance
"""
import argparse
from datetime import timedelta
import string
import re
from collections import OrderedDict


class Step(object):
    """
    The step class is a container of tasks that are group together by the step name. This is ideally used in pipelines
    as multiple jobs will share the same root name if they are of the same step.

    This class holds some aggregate data about the group of tasks assigned to the object.
    """

    def __init__(self, name=''):
        # Variable Members
        self.tasks = []  # type: [Task]
        self.success = 0  # type: int
        self.failure = 0  # type: int

        # Quasi-Constant Members
        self.max_time = None  # type: [timedelta]
        self.max_mem = 0  # type: int
        self.max_vmem = 0  # type: int
        self.name = name  # type: str

    def __str__(self):
        val = '*' * 100 + '\n'
        val += "Step Name: " + self.name + '\n'
        val += "Total jobs found: " + str(len(self.tasks)) + "\n"
        val += "Success: " + str(self.success) + "\tFailure: " + str(self.failure)
        val += "Mem Limit: " + str(self.max_mem) + "\tVmem Limit: " + str(self.max_vmem) + "\n"
        val += "Walltime Limit: " + str(self.max_time) + '\n'
        for task in self.tasks:
            val += str(task)
        return val

    def get_stat(self, func, attr):
        """

        :param func:
        :type func:
        :param attr:
        :type attr:
        :return:
        :rtype:
        """
        attr_vals = filter(None, [val.get(attr) for val in self.tasks])
        if type(attr_vals[0]) == timedelta:
            attr_vals = [val.total_seconds() for val in attr_vals]
        else:
            attr_vals = [get_bytes(val) for val in attr_vals]
        return func(attr_vals)

    def row_print(self):
        """

        :return:
        :rtype:
        """
        row = "{:^15}" + "{:>15}" * 15
        stat_general = [self.name, len(self.tasks), self.success, self.failure]
        stat_time = [self.max_time, timedelta(self.get_stat(max, 'Wallclock Duration')),
                     timedelta(self.get_stat(min, 'Wallclock Duration'))]
        stat_vmem = [self.max_vmem, self.get_stat(max, 'vmem Used'), self.get_stat(min, 'vmem Used')]
        stat_pmem = [self.max_mem, self.get_stat(max, 'Memory Used'), self.get_stat(min, 'Memory Used')]
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
            self.success += 1
        else:
            self.failure += 1

        # Add task to list and generate new limit information
        self.tasks.append(inst)
        self.set_limits()

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
            for i, task in enumerate(self.tasks):
                check_name = task_name == task.get('Job Name')
                check_id = str(task_id) == task.get('Job Id')

                if (task_name and task_id and check_name and check_id) or \
                        (task_name and not task_id and check_name) or \
                        (not task_name and task_id and check_id):
                    if check_name and check_id:
                        if task.get('Exit Code') == 0:
                            self.success -= 1
                        else:
                            self.failure -= 1
                        del self.tasks[i]
                        self.set_limits()

    def set_limits(self):
        """
        Determines the highest PBS resource limit requested from the list of tasks.

        :rtype: None
        """
        self.max_time = max([get_delta(task.get('Wallclock Limit')) for task in self.tasks])
        self.max_mem = max([get_bytes(task.get('Memory Limit')) for task in self.tasks])
        self.max_vmem = max([get_bytes(task.get('vmem Limit')) for task in self.tasks])


class Task(object):
    """
    The Task class holds runtime statistics of a particular job from PBS/Torque. The various metrics are stored
    as keys to a dictionary and the values from the job are matched to its related key.
    """

    def __init__(self, section=None):
        self.data = OrderedDict()
        for attr, value in [line.split(':') for line in section]:
            self.data[attr.strip()] = value.strip()
        self.name = self.data.setdefault('Job Name')
        self.name += self.data.setdefault('Job Id')
        self.name += self.data.setdefault('Start Time')

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
        if attr not in self.data:
            raise KeyError('Attribute ' + attr + ' is not in this task. Use one of ' + ', '.join(self.data.keys()))
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
    metrics = ["Step Name", "Total Runs", "Successful Runs", "Failed Runs",
               "Walltime Limit", "Max Walltime", "Average Walltime", "Min Walltime",
               "V. Mem Limit", "Max V. Mem", "Average V. Mem", "Min V. Mem",
               "Memory Limit", "Max Memory", "Average Memory", "Min Memory"]
    header = ("{:^15}" * 16).format(*metrics)

    def __init__(self, step_list, name=None, pipeline=None):
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
        new_step = list()
        for step in self.pipeline:
            new_step += [item for item in self.steps if item.name == step]
        if new_step:
            self.steps = new_step
        else:
            raise ValueError('No steps found that matches the steps in your pipeline')

    def print_table(self):
        """

        :return:
        :rtype:
        """
        out_str = self.header  # type: str
        for step in self.steps:
            out_str += step.row_print()
        print (self.template.format(pipeline_name=self.name, table=out_str))


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
    for item in steps_handle.values():  # type: Step
        item.set_limits()
    return steps_handle


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
    """

    :param size_string:
    :type size_string:
    :return:
    :rtype:
    """
    try:
        size_string = size_string.lower().replace(',', '')
        size = re.search('^(\d+)[a-z]i?b$', size_string).groups()[0]
        suffix = re.search('^\d+([kmgtp])i?b$', size_string).groups()[0]
    except AttributeError:
        raise ValueError("Invalid Input")
    shft = suffix.translate(string.maketrans('kmgtp', '12345')) + '0'
    return int(size) << int(shft)


def get_delta(time_str):
    """

    :param time_str:
    :type time_str: str
    :return:
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
    performance statistics.""",
                                     epilog='By Michael Li, Nov. 2016')
    parser.add_argument('in_file', type=str, metavar='Log', help="""Output from showjobs.""")
    parser.add_argument('pipe_steps', type=str, metavar='File', help="""Names of steps of a pipeline in a file""")
    args = parser.parse_args()

    steps = parse_pbs_job_data(args.in_file)
    
    with open(args.pipe_steps) as f:
        pipeline = f.readlines()
    StatsManager(steps, '', pipeline)
