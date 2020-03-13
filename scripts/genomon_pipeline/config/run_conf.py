#! /usr/bin/env python

class Run_conf(object):
    """
    class for job related parameters
    """

    def __init__(self, sample_conf_file = None, 
                        project_root = None, 
                        analysis_type = None,
                        genomon_conf_file = None):

        self.sample_conf_file = sample_conf_file
        self.project_root = project_root
        self.genomon_conf_file = genomon_conf_file 
