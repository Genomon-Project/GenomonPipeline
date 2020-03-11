#! /usr/bin/env python

import os, pwd
import configparser
from genomon_pipeline.__init__ import __version__

class Genomon_conf(object):
    """
    class for job related parameters
    """

    def __init__(self, conf, analysis_date):
            
        self.software_version ={'genomon_pipeline':'genomon_pipeline-'+__version__} 
        self.analysis_date = analysis_date
        self.genomon_conf = configparser.SafeConfigParser()
        self.genomon_conf.read(conf)
        
    def _conf_check(self, target_section = None, target_option = None):
        err_msg = 'No target File : \'%s\' for the %s key in the section of %s' 
        
        def __path_check(section, option):
            value = self.genomon_conf.get(section, option)
            if value == "":
                return True
            if os.path.exists(value):
                return True
            raise ValueError(err_msg % (value, option, section))
    
        if target_section == None:
            for section in self.genomon_conf.sections():
                if target_option == None:
                    for option in self.genomon_conf[section]:
                        __path_check(section, option)
                else:
                    if target_option in self.genomon_conf[section]:
                        __path_check(section, target_option)                  
        else:
            if target_option == None:
                for option in self.genomon_conf[target_section]:
                    __path_check(target_section, option)
            else:
                __path_check(target_section, target_option) 
                        
    def genomon_conf_check(self):
    
        self._conf_check(target_section = "REFERENCE")
        self._conf_check(target_option = "image")
    
    def software_version_set(self):
        for section in self.genomon_conf.sections():
            if "image" in self.genomon_conf[section]:
                value = self.genomon_conf.get(section, "image")
                if value == "":
                    continue
                image = value.replace(".simg", "").split("/")[-1]
                self.software_version[section] = image
        
    def get_version(self, key):
        return self.software_version[key]
    
    def get_meta_info(self, softwares):
    
        print_meta_info = "# Version: " + ' '.join([self.software_version[x] for x in softwares])
        print_meta_info = print_meta_info + '\n' + "# Analysis Date: " + self.analysis_date
        print_meta_info = print_meta_info + '\n' + "# User: " + pwd.getpwuid(os.getuid()).pw_name
       
        return print_meta_info

