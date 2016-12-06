#!/usr/bin/env python


from core.config import config
from core.job import Job


def cmscan(rfam_path, query, tblout, log_path, name=None):
    return Job(name=name,
               input_files=[rfam_path, query],
               output_files=[tblout, log_path],
               module_entries=[['cmscan', 'infernal']],
               command='cmscan -o {log} '
                       '--tblout {tblout} '
                       '{options} '
                       '{rfam_path} '
                       '{query}'.format(log=log_path,
                                        tblout=tblout,
                                        options=config.param('cmscan', 'options', required=False),
                                        rfam_path=rfam_path,
                                        query=query))