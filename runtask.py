#!/usr/bin/python
"""

Genome Analysis Pipeline
RunTask class


"""

import subprocess
import genomon_rc as rc

#
# Run process
#
class RunTask:

    def __init__( self, enable_mpi = False, ncpus = 0, log = None ):
        """
        Constructor

        """
        self.enable_mpi = enable_mpi
        self.log = log

        if enable_mpi:
            self.ncpus = ncpus
            self.comm = []

            from mpi4py import MPI
            self.MPI = MPI

            self.log.info( "RunTask: mpi={enable_mpi}".format( enable_mpi = enable_mpi ) )
            self.log.info( "RunTask: ncpus={cpus}".format( ncpus = ncpus ) )


    def __del__( self ):
        """
        Destructor

        """
        if self.enable_mpi:
            self.disconnect()


    def runtask( self, job_type, memory, run_cmd ):
        """
        Front end funtion to run task

        """

        if self.enable_mpi:
            self.__runtask_by_mpi( job_type, memory, run_cmd )
        else:
            self.__runtask_by_qsub( job_type, memory, run_cmd )

    def disconnect( self ):
        id = 0
        for comm_tmp in self.comm:
            self.log.info( "Disconnect process {id}".format( id = id ) )
            comm_tmp.Disconnect()
            id += 1

    def __runtask_by_mpi( self, job_type, memory, run_cmd ):

        return_code = 0

        try:
            self.comm.append( self.MPI.COMM_SELF.Spawn( sys.executable, args=['mpi_worker.py', run_cmd ], maxprocs=1 ) )

        except OSError as e:
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "OS error." )
            return_code = 1

        except IOErr as (errno, strerror):
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "IOError {0}{1}",format( errno, strerror ) )
            return_code = 1

        except:
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "Unexpected error." )
            return_code = 1
                
            return_code = 1

        return return_code

    def __runtask_by_qsub( self, job_type, memory, run_cmd ):
        """
        Submit a job by qsub.

        """
        self.log.info( '# runtask_by_qsub' )
        self.log.info( "command = {cmd}".format(cmd = run_cmd) )
        self.log.info( "memory  = {mem}".format(mem = memory) )
        self.log.info( "job     = {job}\n".format(job = job_type) )

        if job_type == 'mjob':
            job_type = ''

        cmd_tmp = rc.qsub_cmd.format(
                            s_vmem  = memory,
                            mem_req = memory[:-1],
                            job_type = job_type,
                            cmd     = run_cmd )

        return_code = 0
        std_out = None
        std_err = None
        try:
            process = subprocess.Popen( cmd_tmp,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE )

            std_out, std_err = process.communicate()

        except OSError as e:
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "OS error." )
            return_code = 1

        except IOErr as (errno, strerror):
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "IOError {0}{1}",format( errno, strerror ) )
            return_code = 1

        except:
            self.log.error( "RunTask.runtaskby_qsub failed." )
            self.log.error( "Unexpected error." )
            return_code = 1
                
        self.log.info( "STDOUT: {stdout}".format( stdout = std_out ) )
        self.log.info( "STDERR: {stderr}".format( stderr = std_err ) )

        return return_code
