#! /usr/bin/env python

import os, pwd, sys
if sys.version_info.major == 2:
    import ConfigParser as cp
else:
    import configparser as cp
from genomon_pipeline.__init__ import __version__
import genomon_pipeline.config.run_conf as rc

global genomon_conf

genomon_conf = cp.SafeConfigParser()

software_version ={'genomon_pipeline':'genomon_pipeline-'+__version__}

dna_software_version = {} 
rna_software_version = {} 

err_msg = 'No target File : \'%s\' for the %s key in the section of %s' 


def _conf_check(target_section = None, target_option = None):
    def __path_check(section, option):
        value = genomon_conf.get(section, option)
        if value == "":
            return True
        if os.path.exists(value):
            return True
        raise ValueError(err_msg % (value, option, section))

    if target_section == None:
        for section in genomon_conf.sections():
            if target_option == None:
                for option in genomon_conf[section]:
                    __path_check(section, option)
            else:
                if target_option in genomon_conf[section]:
                    __path_check(section, target_option)                  
    else:
        if target_option == None:
            for option in genomon_conf[target_section]:
                __path_check(target_section, option)
        else:
            __path_check(target_option, target_option) 
                    
def dna_genomon_conf_check(section):
    """
    function for checking the validity of genomon_conf for DNA analysis
    """

    _conf_check(target_section = "REFERENCE")
    _conf_check(target_option = "image")
    
    if genomon_conf.has_option("annotation", "active_annovar_flag") :
        if genomon_conf.get("annotation", "active_annovar_flag"):
            _conf_check(target_section = "SOFTWARE", target_option = "annovar")
                

def rna_genomon_conf_check():
    """
    function for checking the validity of genomon_conf for RNA analysis
    """

    _conf_check(target_section = "REFERENCE")
    _conf_check(target_option = "image")

def _image_version():
    for section in genomon_conf.sections():
        if "image" in genomon_conf[section]:
            value = genomon_conf.get(section, "image")
            if value == "":
                continue
            image = value.replace(".simg", "").split("/")[-1]
            software_version[section] = image

def dna_software_version_set():
    _image_version()
    
def rna_software_version_set():
    _image_version()
    
def get_version(key):
    return software_version[key]

def get_meta_info(softwares):

    print_meta_info = "# Version: " + ' '.join([software_version[x] for x in softwares])
    print_meta_info = print_meta_info + '\n' + "# Analysis Date: " + rc.run_conf.analysis_date
    print_meta_info = print_meta_info + '\n' + "# User: " + pwd.getpwuid(os.getuid()).pw_name
   
    return print_meta_info

