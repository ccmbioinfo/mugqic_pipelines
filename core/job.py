#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import collections
import datetime
import logging
import os
import ntpath

# MUGQIC Modules
from config import *

log = logging.getLogger(__name__)

class Job:

    def __init__(self, input_files=[], output_files=[], module_entries = [], name="", command="", report_files=[], removable_files=[], local=""):
        # Remove undefined input/output/removable files if any
        self._input_files = filter(None, input_files)
        self._output_files = filter(None, output_files)
        self._report_files = filter(None, report_files)
        self._removable_files = filter(None, removable_files)

        # Retrieve modules from config, removing duplicates but keeping the order
        self._modules = list(collections.OrderedDict.fromkeys([config.param(section, option) for section, option in module_entries]))
        self._name = name
        self._local = True if (local.lower() in ["yes", "true", "y", "t", "1"]) else False

        # Get the references for add_localhd before command is called
        if self._local:
            sections = list(set([section for section, option in module_entries]))
            self._refs = [config.param(section, 'ref_dirs', required=False) for section in sections]

        self._command = command
        self._orig_command = command # Original command kept for piping jobs

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._name

    @property
    def output_dir(self):
        return self._output_dir

    @property
    def input_files(self):
        return self._input_files

    @property
    def output_files(self):
        return self._output_files

    @property
    def report_files(self):
        return self._report_files

    @property
    def removable_files(self):
        return self._removable_files

    @property
    def done(self):
        return self._done

    @property
    def dependency_jobs(self):
        return self._dependency_jobs

    @property
    def modules(self):
        return self._modules

    @property
    def local(self):
        return self._local

    @property
    def command(self):
        return self._command

    @property
    def orig_command(self):
        return self._orig_command

    @property
    def command_with_modules(self):
        command = self.command
        if self.modules:
            command = "module load " + " ".join(self.modules) + " && \\\n" + command
        return command

    # Modify the command when it's been set
    def __setattr__(self, name, value):
        self.__dict__[name] = value

        # If the command is being changed...
        if name == "command" or name == "_command":
            if self.local:
                self._orig_command = value
                self.__dict__[name] = self.add_localhd(value)
        # If the input files associated with the job are being changed... Redo the command setting
        if name == "input_files" or name == "_input_files":
            if ("_local" in self.__dict__) and (self.local):
                    self.command = self.add_localhd(self.orig_command)

    def abspath(self, file):
        tmp_file = os.path.expandvars(file)
        if not os.path.isabs(tmp_file):
            # File path is relative to the job output directory
            tmp_file = os.path.normpath(os.path.join(self.output_dir, tmp_file))
        return tmp_file

    def is_up2date(self):
        # If job has dependencies, job is not up to date
        if self.dependency_jobs:
            log.debug("Job " + self.name + " NOT up to date")
            log.debug("Dependency jobs:\n  " + "\n  ".join([job.name for job in self.dependency_jobs]) + "\n")
            return False

        # Retrieve absolute paths for .done, input and output files to avoid redundant OS function calls
        abspath_done = self.abspath(self.done)
        abspath_input_files = [self.abspath(input_file) for input_file in self.input_files]
        abspath_output_files = [self.abspath(output_file) for output_file in self.output_files]

        # If any .done, input or output file is missing, job is not up to date
        for file in [abspath_done] + abspath_input_files + abspath_output_files:
            # Use 'exists' instead of 'isfile' since input/output files can be directories
            if not os.path.exists(file):
                log.debug("Job " + self.name + " NOT up to date")
                log.debug("Input, output or .done file missing: " + file)
                return False

        # Remove any files in the list of input files if it also appears in the output file list
        # This is because the time of a file that is in the input and output will cause issues with the job checking
        if config.param('DEFAULT', 'rm_infile_if_out', required=False).lower() in ['yes', 'true', 't', '1']:
            log.debug("Not checking for time integrity if an output file is the same as the input file")
            abspath_input_files = [file for file in abspath_input_files if file not in abspath_output_files]

        # Retrieve latest input file by modification time i.e. maximum stat mtime
        # Use lstat to avoid following symbolic links
        latest_input_file = max(abspath_input_files, key=lambda input_file: os.lstat(input_file).st_mtime)
        latest_input_time = os.lstat(latest_input_file).st_mtime

        # Same with earliest output file by modification time
        earliest_output_file = min(abspath_output_files, key=lambda output_file:os.lstat(output_file).st_mtime)
        earliest_output_time = os.lstat(earliest_output_file).st_mtime

        # If any input file is strictly more recent than all output files, job is not up to date
        if latest_input_time > earliest_output_time:
            log.debug("Job " + self.name + " NOT up to date")
            log.debug("Latest input file modification time: " + latest_input_file + " " + datetime.datetime.fromtimestamp(latest_input_time).isoformat() + " > earliest output file modification time: " + earliest_output_file + " " + datetime.datetime.fromtimestamp(earliest_output_time).isoformat() + "\n")
            return False

        # If all previous tests passed, job is up to date
        return True

    # Helper function for add_localhd
    # in_files = input files with relative paths (ignore absolute paths)
    def setup_localhd(self, in_files):
        # Get the name and directory for bam files so can get bai files
        bam_files = [x.rsplit('.',1)[0] for x in in_files if x.endswith('.bam')]
        new_cmd = ""

        # Copy files to the local space (if any relative paths)
        if in_files:
            new_cmd = new_cmd+"""\
cp -r -L --parents {input_files} $TMPDIR/ && \\\n""".format(input_files = " ".join(in_files))

            # Copy the bai files for any bam files that are copied over - don't want to && b/c egrep returning nothing is okay
            for bam_file in bam_files:
                new_cmd = new_cmd+"""\
to_copy=$(ls {bam_file}.* | egrep '{bam_file}.(bam.)?bai')
rm_copy="$rm_copy $to_copy"
if [[ -n $to_copy ]]; then cp --parents $to_copy $TMPDIR/; fi && \\\n""".format(bam_file = bam_file)

        # Move to the TMPDIR directory
        new_cmd = new_cmd+ "cd $TMPDIR && \\\n"

        # Create the subdirectories in $TMPDIR for output
        output_directories = []
        for out_file in self.output_files:
            out_dir = ntpath.split(out_file)[0]
            if out_dir not in output_directories:
                new_cmd = new_cmd+"""mkdir -p {output_dir} && \\\n""".format(output_dir = out_dir)
                output_directories.append(out_dir)

        # Return the now set up new command
        return new_cmd


    # Add local-hd step
    def add_localhd(self, command):
        # Before modify - make sure input_files, output_files, and command not null
        if not(self.input_files) or not(self.output_files) or not(command):
            return command

        # Make sure localhd hasn't been implemented yet
        if (command.find("$TMPDIR") >= 0) or \
                (command.find(config.param('DEFAULT', 'tmp_dir')) >= 0):
            return command

        # Initialization
        # Only take note of relative input paths (ignore absolute)
        input_files = [x for x in self.input_files if not x.startswith('/')]
        mod_cmd = command
        new_cmd = self.setup_localhd(input_files)

        # If the command makes a directory at the beginning, put before moving to $TMPDIR
        mk_idx = mod_cmd.find("mkdir") if mod_cmd.find("mkdir") == 0 else -1
        ret_idx = mod_cmd.find("\n", mk_idx) if mk_idx >= 0 else -1
        if ret_idx >= 0:
            new_cmd = ''.join([mod_cmd[mk_idx:ret_idx+1], new_cmd])
            mod_cmd = ''.join([mod_cmd[:mk_idx], mod_cmd[ret_idx+1:]])

        # Add beginning path to reference directories (defined in config file)
        for ref in self._refs:
            for rdir in ref.split(","):
                if rdir:
                    rname = "$PBS_O_INITDIR/"+rdir.lstrip('/')
                    mod_cmd = re.sub(r"(\s|=)"+rdir+r"(/|\s)", r"\1"+rname+r"\2", mod_cmd)

        # Copy files back to the original space
        new_cmd = ''.join([new_cmd, mod_cmd, " && \\\n"])
        new_cmd = new_cmd + """\
rm -rf {remove_files} $rm_copy && \\
cp -r $TMPDIR/* $PBS_O_INITDIR/ && \\
cd $PBS_O_INITDIR""".format(
            remove_files=" ".join(input_files)
        )

        # Return the new command
        return new_cmd


# Create a new job by concatenating a list of jobs together
def concat_jobs(jobs, name=""):

    # Merge all input/output/report/removable files and modules
    input_files = []
    output_files = []
    report_files = []
    removable_files = []
    modules = []
    full_input_files = []
    for job_item in jobs:
        input_files.extend([input_file for input_file in job_item.input_files if input_file not in input_files and input_file not in output_files])
        output_files.extend([output_file for output_file in job_item.output_files if output_file not in output_files])
        report_files.extend([report_file for report_file in job_item.report_files if report_file not in report_files])
        removable_files.extend([removable_file for removable_file in job_item.removable_files if removable_file not in removable_files])
        modules.extend([module for module in job_item.modules if module not in modules])

        # If using localhd - update the input files
        full_input_files.extend([input_file for input_file in job_item.input_files if input_file not in full_input_files])
        if job_item.local:
            job_item.input_files = full_input_files

    job = Job(input_files, output_files, name=name, report_files=report_files, removable_files=removable_files)
    job.modules = modules

    # Merge commands
    job.command = " && \\\n".join([job_item.command for job_item in jobs])

    return job

# Create a new job by piping a list of jobs together
def pipe_jobs(jobs, name=""):

    job = Job(jobs[0].input_files, jobs[-1].output_files, name=name)

    # Merge all report/removable files and modules
    report_files = []
    removable_files = []
    modules = []
    for job_item in jobs:
        report_files.extend(job_item.report_files)
        removable_files.extend(job_item.removable_files)
        modules.extend(job_item.modules)

    # Remove duplicates if any, keeping the order
    report_files = list(collections.OrderedDict.fromkeys([report_file for report_file in report_files]))
    job.report_files = report_files
    removable_files = list(collections.OrderedDict.fromkeys([removable_file for removable_file in removable_files]))
    job.removable_files = removable_files
    modules = list(collections.OrderedDict.fromkeys([module for module in modules]))
    job.modules = modules

    # Merge commands
    job.command = " | \\\n".join([job_item.orig_command for job_item in jobs])

    return job


# Create a job that creates the directory needed for output_file
def mkdir(output_file):
    return Job(name='mkdir.{}'.format(output_file),
               command='mkdir -p {dir}'.format(dir=os.path.dirname(output_file)))
