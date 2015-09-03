"""
Utilty functions

"""
import os
import inspect
from datetime import datetime
from glob import glob

from resource import genomon_rc as res

#####################################################################
#
# Private modules
#
########################################
def whoami():
    return inspect.stack()[1][3]

########################################
def whosdaddy():
    return inspect.stack()[2][3]

########################################
def split_file( file_name, output_file_name, split_len ):
    """
    Split files by line number

    """

    input = open( file_name, 'r' )
    count = 0
    at = 0
    dest = None

    for line in input:
        if count % split_len == 0:
            if dest: dest.close()
            dest = open( output_file_name.format( num=at ) )
            at += 1
        dest.write( line )
        count += 1

    dest.close()
    input.close()


########################################
def make_script_file_name( function_name, Geno ):
    """
    Split files by line number

    """

    now = datetime.now()
    shell_script_name = res.file_timestamp_format.format(
                                    name=function_name,
                                    year=now.year,
                                    month=now.month,
                                    day=now.day,
                                    hour=now.hour,
                                    min=now.minute,
                                    msecond=now.microsecond )

    if not Geno.options.abpath:
        script_dir = Geno.dir[ 'project_root' ] + '/' + Geno.dir[ 'script' ]
    else:
        script_dir = Geno.dir[ 'script' ]

    shell_script_full_path = "{script}/{file}.sh".format(
                                    script = script_dir,
                                    file = shell_script_name )

    return shell_script_full_path

########################################
def make_sample_file_name( filename,
                           file_fmt,
                           dir = None,
                           subdir = None,
                           ext = None ):
    """
    {dir}
    {base}
    {ext}
    """
    tmp_name = os.path.splitext( os.path.basename( filename ) )
    basename = tmp_name[ 0 ]
    if ext == None:
        ext = tmp_name[ 1 ]
    if not dir:
        dir = os.path.dirname( filename )

    if subdir:
        ret_name = file_fmt.format(
                        dir = dir,
                        subdir = subdir,
                        base = basename,
                        ext = ext )
    else:
        ret_name = file_fmt.format(
                        dir = dir,
                        base = basename,
                        ext = ext )
    return ret_name

########################################
def replace_reserved_string( dir_tmp, cwd, Geno ):
    """
    Reserved names to replace strings defined in job configuration file
        project_directory   -> defined as project in job configuration file
        sample_date         -> defined sample_date in job configuration file
        sample_group        -> defined sample_group in job configuration file
        analysis_date       -> date of the pipeline to run
    """
    #
    # Replace reserved strings
    #
    dir_replace = None
    if dir_tmp == 'project_directory':
        pass
    elif dir_tmp == 'sample_date':
        dir_replace = str( Geno.job.get_job( 'sample_date' ) )
    elif dir_tmp == 'sample_group':
        dir_replace = str(Geno.job.get_job( 'sample_group' ))
    elif dir_tmp == 'analysis_date':
        dir_replace = str( Geno.job.get_job( 'analysis_date' ) )
        if dir_replace == 'today' or Geno.options.sample_file:
            dir_replace = res.date_format.format( 
                        year = Geno.now.year, 
                        month = Geno.now.month, 
                        day = Geno.now.day ) 
    else:
        dir_replace = dir_tmp
    
    return dir_replace

def make_dir( dir, Geno ):
    if not os.path.exists( dir ):
        os.makedirs( dir )
        os.chmod( dir, Geno.dir_mode )
    return dir

def get_dir ( dir_tree, cwd, dir_name, Geno ):
    """
    return the path to the specified directory by dir_tree

    """
    if isinstance( dir_tree, dict ):
        for dir_tmp in dir_tree.keys():
            cwd_tmp = cwd
            dir_replace = replace_reserved_string( dir_tmp, cwd, Geno )
            if dir_replace:
                 cwd_tmp += '/' + dir_replace

            if dir_tmp == dir_name:
                return cwd_tmp

            dir_returned =  get_dir( dir_tree[ dir_tmp ], cwd_tmp, dir_name, Geno  )

            if None != dir_returned:
                return dir_returned
            
    elif isinstance( dir_tree, list ):
        n = 0
        for dir_tmp in dir_tree:
            if isinstance( dir_tmp, str):
                cwd_tmp = cwd
                dir_replace = replace_reserved_string( dir_tmp, cwd, Geno )
                if dir_replace:
                    cwd_tmp += '/' + dir_replace
            elif isinstance( dir_tmp, dict):
                cwd_tmp = cwd
                dir_replace = replace_reserved_string( dir_tmp.keys()[ 0 ] , cwd, Geno )
                if dir_replace:
                    cwd_tmp += '/' + dir_replace

            if ( ( isinstance( dir_tmp, str) and dir_tmp == dir_name ) or
                 ( isinstance( dir_tmp, dict) and dir_tmp.keys()[0] == dir_name ) ):
                return cwd_tmp
            else:
                if ( isinstance( dir_tree[ n ], dict ) or
                     isinstance( dir_tree[ n ], list ) ):
                    dir_returned =  get_dir( dir_tree[ n ], cwd, dir_name, Geno )
                    if None != dir_returned:
                        return dir_returned
            n = n + 1
    else:
        if isinstance( dir_tree, str) and dir_tree == dir_name:
            cwd_tmp = cwd
            dir_replace = replace_reserved_string( dir_tree, cwd, Geno )
            if dir_replace != '':
                cwd_tmp += '/' + dir_replace
            return cwd_tmp

    return None

def make_input_target( dir_tree, cwd, Geno, sample_name = None, sample_list = None ):

    dir = Geno.job.get_job( 'input_file_dir' ),
    if sample_list:
        subdir_list = sample_list

    elif dir and sample_name:
        subdir_list = glob( "{dir}/{subdir}".format(
                                dir = Geno.job.get_job( 'input_file_dir' ),
                                subdir = sample_name ) )

    for target_dir in res.end_dir_list:
        tmp_dir = get_dir( dir_tree, cwd, target_dir, Geno )
        if tmp_dir:
            Geno.dir[ target_dir ] = tmp_dir

            if target_dir in res.dir_nontask_list:
                make_dir( Geno.dir[ target_dir ], Geno )

            elif ( target_dir in res.dir_task_list and
                 set( res.dir_task_list[ target_dir ] ).intersection( Geno.job.get_job( 'tasks' )[ Geno.job_tasks ] ) ):
                make_dir( Geno.dir[ target_dir ], Geno )
                if sample_list:
                    for subdir_tmp in subdir_list:
                        make_dir( "{dir}/{subdir}".format( 
                                    dir = Geno.dir[ target_dir ],
                                    subdir = subdir_tmp ),
                            Geno )
                elif sample_name :
                    for subdir_tmp in [ x for x in subdir_list if os.path.isdir( x ) ]:
                        make_dir( "{dir}/{subdir}".format(
                                    dir = Geno.dir[ target_dir ],
                                    subdir = os.path.basename( subdir_tmp ) ),
                                  Geno )

def make_script_file( function_name, string, Geno, **kwargs ):
    shell_script_full_path = make_script_file_name( function_name, Geno )
    shell_script_file = open( shell_script_full_path, 'w' )
    shell_script_file.write( string.format( **kwargs))
    shell_script_file.close()
    return shell_script_full_path

def save_status_of_this_process( process_name, output_file, return_code, Geno, use_subdir ):

    Geno.status.save_status( process_name, output_file, return_code, use_subdir = use_subdir )

def get_status_of_this_process( process_name, output_file, Geno, use_subdir ):

    exit_status = Geno.status.check_exit_status(
                    process_name,
                    output_file,
                    use_subdir = use_subdir ) 

    return exit_status