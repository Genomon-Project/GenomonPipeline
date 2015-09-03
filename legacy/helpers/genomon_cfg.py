#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015
import sys
import os
import ConfigParser
import yaml
import job_check

"""
    Genomon system configuration file parser
    How to use it:
        ge_cfg = genomon_config( 'path to file' )
        ge_cfg.get( 'REFERENCE', 'hg19_fasta' )

"""
class genomon_config( object ):
    def __init__( self, config_file = None, log = None):

        self.__log = log
        if config_file != None:
            self.__config_file = config_file
            self.open_cfg( config_file )

    def open_cfg( self, config_file ):
        self.__config_file = config_file

        try:
            if self.__config_file != None:
                self.__conf = ConfigParser.SafeConfigParser()
                self.__conf.read( self.__config_file )
            else:
                self.__log.error( "genomon_config.open_cfg: configuration file is not loaded properly." )
                raise

        except IOError as (errno, stderror ):
            self.__log.error( "genomon_config.open_cfg: IOError: error number: {num}, std_error: {stderr}".format(
                        num = errno, stderr = stderror ) )
            raise
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.__log.error( "genomon_config.open_cfg: Unexpected error" )
            self.__log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            raise


    def get( self, section, item ):
        try:
            return_string = ''
            if self.__conf != None:
                return_string =  self.__conf.get( section, item )
                if None == return_string:
                    self.__log.error( "genomon_config.get: {sec}: {item} is not defined in system config file.".format(
                                    sec = section,
                                    item = item
                                ))
            return return_string

        except:
            self.__log.error( "genomon_config.get: configuration file is not loaded properly." )
            return None

    def check_file( self, keyword_file ):
        try:
            f = open( keyword_file )
            f_yaml = yaml.load( f )
            if job_check.System_config_file_check( self.__conf, f_yaml ):
                return_value = True
            else:
                return_value = False

            f.close()

            return return_value

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.__log.error( "genomon_job.open_job: unexpected error:", sys.exc_info()[0] )
            self.__log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            raise

    def items( self, section ):
        try:
            return_list = []
            if self.__conf != None:
                return_list = self.__conf.items( section )
                if [] == return_list:
                    self.__log.error( "genomon_config.get: {sec}: {item} is not defined in system config file.".format(
                                    sec = section,
                                    item = item
                                ))
            return return_string

        except:
            self.__log.error( "genomon_config.get: configuration file is not loaded properly." )
            return None