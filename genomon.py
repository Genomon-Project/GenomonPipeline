#!/usr/bin/python
"""

Genome Analysis Pipeline


"""

import sys
import os
import shutil
from datetime import datetime
import argparse
from ruffus import *
import yaml

################################################################################
#
# Globals
#
global log
log = None
global Geno
class Geno(object):
    now     = None
    options = None
    conf    = None
    job     = None
    dir     = {}
    RT      = None
    cwd     = None
    dir_mode = 0755
    status = None


################################################################################
#
# Private modules
#
from resource import genomon_rc as res
from helpers.genomon_cfg import genomon_config as ge_cfg
from helpers.genomon_job import genomon_job as ge_job
from helpers.runtask import RunTask
from helpers.utils import *
from helpers.status import genomon_status as ge_status

################################################################################
#
# Subroutines
#
def construct_arguments( ):
    """
    Call argparse and create argument object
    
    """

    parser = cmdline.get_argparse( description='Genome Analysis Pipeline' )

    ge_arg = parser.add_argument_group( 'genomon', 'Genomon options' );

    ge_arg.add_argument( '-s', "--config_file",  help = "Genomon pipeline configuration file",      type = str )
    ge_arg.add_argument( '-f', "--job_file",     help = "Genomon pipeline job configuration file",  type = str )
    ge_arg.add_argument( '-p', "--param_file",   help = "Genomon pipeline analysis parameter file", type = str )

    ge_arg.add_argument( '-m', "--mpi",          help = "Use MPI job submission",       action ='store_true',
                                                                                        default = False )
    ge_arg.add_argument( '-d', "--drmaa",        help = "Use DRMAA job submission",     action ='store_true',
                                                                                        default = False )
    ge_arg.add_argument( '-l', "--abpath",       help = "Use absolute path in scripts", action ='store_true',
                                                                                        default = False )

    return parser

########################################
def printheader( myself, options ):
    """
    Print infomration about this run

    """
    Geno.now = datetime.now()

    log.info( "\nGenerated by {my}".format(my = myself ) )
    log.info( "Input config file = {input}".format( input = options.config_file  ) )
    log.info( "Input job file    = {input}".format( input = options.job_file  ) )


def make_directories():
    """
    Make Directory Tree Structure
       Read directory structure from resource file
       Make sure that input data exists.
       Create results directory if necesesary

    """

    #
    # Make sure the input data is availalbe.
    #
    

    error_message = ''
    try:
        project_root = os.path.expanduser( Geno.job.get_job( 'project_root' ))
        Geno.dir[ 'project_root' ] = project_root
        if not os.path.exists( project_root ):
            log.error( "Dir: {dir} not found.".format( dir = project_root ) ) 
            raise


        #
        # get directory tree, directory permission
        #
        dir_tree  = Geno.job.get_job( 'project_dir_tree' )
        if not dir_tree:
            dir_tree  = yaml.load( res.dir_tree_resource  )

        dir_permit  = Geno.job.get_job( 'directory_permission' )
        if dir_permit and dir_permit == 'group':
            Geno.dir_mode = 0775
        elif dir_permit and dir_permit == 'all':
            Geno.dir_mode = 0777

        #
        # get directory locations
        #
        Geno.dir[ 'cwd' ] = os.getcwd()
        cwd = project_root
        if not Geno.options.abpath:
            os.chdir( cwd )
            cwd = '.'

        sample_subdir = Geno.job.get_job( 'sample_subdir' )
        if sample_subdir:
            make_input_target( sample_subdir, dir_tree, cwd, Geno )
        else:
            make_input_target( '', dir_tree, cwd, Geno)


        #
        # data diretory
        # make symbolic link from the original input_file_dir
        #   if input_file_dir is not the same as data dir
        #
        make_dir( "{data}/{sample_date}".format(
                                data = get_dir( dir_tree, cwd, 'data', Geno ),
                                sample_date = Geno.job.get_job( 'sample_date' )
                                ),
                  Geno )
        Geno.dir[ 'data' ] = "{data}/{sample_date}/{sample_name}". format(
                                data = get_dir( dir_tree, cwd, 'data', Geno ),
                                sample_date = Geno.job.get_job( 'sample_date' ),
                                sample_name = Geno.job.get_job( 'sample_name' ) )
        if ( not os.path.exists( Geno.dir[ 'data' ] ) and
             Geno.job.get_job( 'input_file_dir' ) != Geno.dir[ 'data' ] ):
                os.symlink( Geno.job.get_job( 'input_file_dir' ), Geno.dir[ 'data' ] )

    except IOError, (errno, strerror):
        log.error( "make_directories failed." )
        log.error( "IOError {0}]{1}".format( errno, strerror ) )
        raise


    except Exception, e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "make_directories failed." )
        log.error( "Unexpected error: {1}".format( error_message ) )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        raise

########################################
def copy_config_files():
    """
    Copy genomon system configuration file and job configuration file
    to results directory

    """
    global Geno

    config_dir = Geno.dir[ 'config' ]

    #
    # Backup system configuration file
    #
    header = 'genomon'
    timestamp = res.timestamp_format.format(
                                        year=Geno.now.year,
                                        month=Geno.now.month,
                                        day=Geno.now.day,
                                        hour=Geno.now.hour,
                                        min=Geno.now.minute,
                                        msecond=Geno.now.microsecond )
    for src, ext in ( ( Geno.options.config_file, '.cfg' ),
                      ( Geno.options.job_file, '_job.yaml' ),
                      ( Geno.options.param_file, '_param.yaml' ) ):

        dest = "{dir}/{basename}_{timestamp}{ext}".format(
                                dir = config_dir,
                                basename = header,
                                timestamp = timestamp,
                                ext = ext )
        if not os.path.isabs( src ):
            src = Geno.dir[ 'cwd' ] + '/' + src
        shutil.copyfile( src, dest )

    #
    # Initialize status object
    #
    Geno.status = ge_status( config_dir, header, timestamp)

########################################
def copy_script_files():
    """
    Copy genomon script files to script directory

    """
    global Geno

    for script_file in res.script_files:
        src = "{dir}/{file}".format(
                    dir = Geno.dir[ 'genomon' ],
                    file = script_file )
        dest = "{dir}/{file}".format(
                    dir = Geno.dir[ 'script' ],
                    file = os.path.basename( script_file )
                )
        shutil.copy( src, dest )

    
########################################
def set_env_variables():
    """
    Set environment variables for some tools

    """
    global Geno

    for tool_env in res.env_list.keys():
        env_value = Geno.conf.get( 'ENV', tool_env )
        if None != env_value:
            for env_name in res.env_list[ tool_env ]:
                if env_name in os.environ:
                    tmp = os.environ[ env_name ]
                    os.environ[ env_name ] = tmp + ':' + env_value
                else:
                    os.environ[ env_name ] = env_value
            return_value = True
        else:
            return_value = False

    if 'PYTHONPATH' in os.environ:
        tmp = os.environ[ 'PYTHONPATH' ] +  ':'
    else:
        tmp = ''
    os.environ[ 'PYTHONPATH' ] = tmp + Geno.dir[ 'genomon' ]

    return return_value

########################################
def cleanup_intermediates():
    """
    Clean up intermediate files

    * Split fastq files
    * Unnecessary bam files
    * Split mutation call results

    """
    pass

###############################################################################
#
# main
#
def main():
    global log
    global log_mutex
    global Geno

    try:
        #
        # Argument parse
        #
        argvs = sys.argv
        arg_parser = construct_arguments()

        #
        # parse arguments
        #
        if len(argvs) < 3:
            arg_parser.print_help()
            sys.exit( 0 )

        Geno.options = arg_parser.parse_args()
        Geno.dir[ 'genomon' ] = os.path.dirname( os.path.realpath(__file__) )

        #
        # Logging setup
        #
        #  logger which can be passed to multiprocessing ruffus tasks
        if Geno.options.verbose:
            verbose_level = Geno.options.verbose
        else:
            verbose_level = 0

        log, log_mutex = cmdline.setup_logging( __name__,
                                                Geno.options.log_file,
                                                verbose_level)
        #
        # Print header in log
        #
        printheader( argvs[ 0 ], Geno.options )

        #
        # Parse system and job config file
        #
        Geno.conf = ge_cfg( config_file = Geno.options.config_file, log = log )
        Geno.job = ge_job( job_file = Geno.options.job_file,
                           param_file = Geno.options.param_file,
                           log = log )

        #
        # Prepare directory tree for pipeline to run.
        # Copy the input configuration files to results directory
        # Link input_data to project's data directory
        #
        make_directories()
        copy_config_files()
        copy_script_files()

        if not set_env_variables():
            log.error( "Necesesary value in [ENV] is not set in system configuration file." )
            raise
            

        #
        # Initialize RunTask object
        #
        native_param = None
        if Geno.options.mpi and not Geno.options.drmaa:
            run_mode = 'MPI'
        elif Geno.options.drmaa and not Geno.options.mpi:
            run_mode = 'DRMAA'
            native_param = Geno.job.get_job( 'drmaa_native' )
        else:
            run_mode = 'qsub'

        if not Geno.options.abpath:
            drmaa_log_dir = Geno.dir[ 'project_root' ] + '/' + Geno.dir[ 'log' ]
        else:
            drmaa_log_dir = Geno.dir[ 'log' ]

        Geno.RT = RunTask( run_mode = run_mode,
                           drmaa_native = native_param,
                           log_dir = drmaa_log_dir,
                           work_dir = Geno.dir[ 'project_root' ],
                           log = log,
                           ncpus = Geno.options.jobs,
                           qsub_cmd = Geno.job.get_job( 'qsub_cmd' ) )

        #
        # Print information
        #
#        log.info( '# main: process={num}'.format( num = Geno.options.jobs ) )

        #######################################################################
        #
        # Run the defined pipeline
        # Figure out what analysis to run from job configuration file
        #
        job_tasks = Geno.job.get_job( 'tasks' )
        run_flag = False

        if 'WGS' in job_tasks:
            from helpers import wgs_pipeline as pipeline
            run_flag = True

        elif 'RNA' in job_tasks:
            from helpers import rna_pipeline as pipeline
            run_flag = True

        elif 'STAR' in job_tasks:
            from helpers import star_pipeline as pipeline
            run_flag = True

        if not run_flag:
            log.error( "Proper task is not set in job configuration file." )
            raise

#
#       multiprocess
#
#       Optional. The number of processes which should be dedicated to running in parallel independent tasks
#            and jobs within each task. If multiprocess is set to 1, the pipeline will execute in the main process.
#
#       multithread
#
#       Optional. The number of threads which should be dedicated to running in parallel independent tasks
#            and jobs within each task. Should be used only with drmaa.
#            Otherwise the CPython global interpreter lock (GIL) will slow down your pipeline
#
#       verbose            
#
#       Optional parameter indicating the verbosity of the messages sent to logger: (Defaults to level 1 if unspecified)
#
#           level 0 : nothing
#           level 1 : Out-of-date Task names
#           level 2 : All Tasks (including any task function docstrings)
#           level 3 : Out-of-date Jobs in Out-of-date Tasks, no explanation
#           level 4 : Out-of-date Jobs in Out-of-date Tasks, with explanations and warnings
#           level 5 : All Jobs in Out-of-date Tasks, (include only list of up-to-date tasks)
#           level 6 : All jobs in All Tasks whether out of date or not
#           level 10: logs messages useful only for debugging ruffus pipeline code
        pipeline_run(   target_tasks = [ pipeline.last_function ],
                        multiprocess = Geno.options.jobs,
                        exceptions_terminate_immediately = True,
                        log_exceptions = True,
                        checksum_level = 2,
                        logger = log )
                        
#        pipeline_cleanup()

#        pipeline_printout_graph( "flow_{job_type}".format( job_type = job_tasks.keys()[ 0 ] ),
#                                 "jpg",
#                                 [ pipeline.last_function ]
#                )
#        cmdline.run( Geno.options )
        #
        #######################################################################

    except IOError as (errno, strerror):
        log.error( "{0}: I/O error({1}): {2}".format( whoami(), errno, strerror) )

    except ValueError:
        log.error( "{0}: ValueError".format( whoami() ) )

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{0}: Unexpected error".format( whoami() ) )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )

    #
    # Clean up
    #
    cleanup_intermediates()


    sys.exit( 0 )


################################################################################
if __name__ == "__main__":
    main()

