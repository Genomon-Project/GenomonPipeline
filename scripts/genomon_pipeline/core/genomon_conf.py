#! /usr/bin/env python

import os
import pwd
import datetime
import configparser
from genomon_pipeline.__init__ import __version__

date_format = "{year:0>4d}-{month:0>2d}-{day:0>2d} {hour:0>2d}:{min:0>2d}:{second:0>2d}"
timestamp_format = "{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}{second:0>2d}_{msecond:0>6d}"

class Genomon_conf(object):
    """
    class for job related parameters
    """

    def __init__(self, conf, exist_check = True):
        
        self.exist_check = exist_check
        self.software_version ={'genomon_pipeline':'genomon_pipeline-'+__version__} 
        self.genomon_conf = configparser.SafeConfigParser()
        self.genomon_conf.read(conf)
        
        now = datetime.datetime.now()
        self.analysis_date = date_format.format(
             year = now.year,
             month = now.month,
             day = now.day,
             hour = now.hour,
             min = now.minute,
             second = now.second
        )
        self.analysis_timestamp = timestamp_format.format(
             year = now.year,
             month = now.month,
             day = now.day,
             hour = now.hour,
             min = now.minute,
             second = now.second,
             msecond = now.microsecond
        )
        
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

    def get(self, section, option):
        return self.genomon_conf.get(section, option)

    def getboolean(self, section, option):
        return self.genomon_conf.getboolean(section, option)

    def path_get(self, section, option):
        path = self.genomon_conf.get(section, option)
        if self.exist_check and not os.path.exists(path):
            err_msg = "[%s] %s: path %s is not exists" % (section, option, path)
            raise ValueError(err_msg)

        return self.genomon_conf.get(section, option)
    
    def safe_get(self, section, option):
        try:
            return self.genomon_conf.get(section, option)
        except Exception:
            pass
        return None
