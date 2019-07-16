#! /usr/bin/env python

import os
import sys
import datetime
import subprocess

file_timestamp_format = "{name}_{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}{second:0>2d}_{msecond:0>6d}"

class Stage_task(object):

    def __init__(self, qsub_option, use_drmaa_flag):
        self.qsub_option = qsub_option
        self.drmaa = use_drmaa_flag
        self.retry_count = 2


    def task_exec(self, arguments, log_dir, script_dir, max_task=0):
        # Make shell script

        now = datetime.datetime.now()
        shell_script_name = file_timestamp_format.format(
                                 name=self.task_name,
                                 year=now.year,
                                 month=now.month,
                                 day=now.day,
                                 hour=now.hour,
                                 min=now.minute,
                                 second=now.second,
                                 msecond=now.microsecond )
        
        shell_script_full_path = "{script}/{file}.sh".format(script = script_dir, file = shell_script_name)
        shell_script_file = open(shell_script_full_path, 'w')
        shell_script_file.write(self.script_template.format(**arguments))
        shell_script_file.close()

        if self.drmaa:
            import drmaa
        
            s = drmaa.Session()
            s.initialize()
        
            jt = s.createJobTemplate()
            jt.jobName = shell_script_name
            jt.outputPath = ':' + log_dir
            jt.errorPath = ':' + log_dir
            jt.nativeSpecification = self.qsub_option
            jt.remoteCommand = shell_script_full_path
            if sys.version_info.major == 2:
                mode = 0750
            else:
                mode = 0o750
            os.chmod(shell_script_full_path, mode)

            returncode = 0
            returnflag = True
            if max_task == 0:
                for var in range(0, (self.retry_count+1)):
                    jobid = s.runJob(jt)
                    returncode = 0
                    returnflag = True
                    now = datetime.datetime.now()
                    date = now.strftime("%Y-%m-%d %H:%M:%S")
                    print >> sys.stderr, "Your job has been submitted with id: " + jobid + " at Date/Time: " + date
                    retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                    now = datetime.datetime.now()
                    date = now.strftime("%Y-%m-%d %H:%M:%S")
                    print >> sys.stderr, "Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date
                    returncode = retval.exitStatus
                    returnflag = retval.hasExited
                    if returncode == 0 and returnflag: break
                s.deleteJobTemplate(jt)
                s.exit()

            else:
                joblist = s.runBulkJobs(jt,1,max_task,1)
                all_jobids = []
                for var in range(0, (self.retry_count+1)):
                    if len(all_jobids) > 0:
                        joblist = all_jobids
                        all_jobids = []
                    returncode = 0
                    returnflag = True
                    now = datetime.datetime.now()
                    date = now.strftime("%Y-%m-%d %H:%M:%S")
                    print >> sys.stderr, 'Your job has been submitted with id ' + str(joblist) + " at Date/Time: " + date
                    s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
                    for curjob in joblist:
                        print >> sys.stderr, 'Collecting job ' + curjob
                        retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                        now = datetime.datetime.now()
                        date = now.strftime("%Y-%m-%d %H:%M:%S")
                        print >> sys.stderr, "Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date
                        
                        if retval.exitStatus != 0 or not retval.hasExited:
                            returncode = retval.exitStatus
                            returnflag = retval.hasExited
                            if var == self.retry_count: break
                            jobId_list = ((retval.jobId).encode('utf-8')).split(".")
                            taskId = int(jobId_list[1])
                            all_jobids.extend(s.runBulkJobs(jt,taskId,taskId,1))
                       
                    if returncode == 0 and returnflag: break
                s.deleteJobTemplate(jt)
                s.exit()

            if returncode != 0 or not returnflag: 
                raise RuntimeError("Job: " + str(retval.jobId)  + ' failed at Date/Time: ' + date)

        else:
            qsub_commands = ['qsub', '-sync', 'yes']
            if max_task != 0:
                qsub_commands.extend(['-t', '1-'+str(max_task)+':1'])

            qsub_options = self.qsub_option.split(' ')
            returncode = subprocess.call(qsub_commands + qsub_options + [shell_script_full_path])

            if returncode != 0: 
                raise RuntimeError("The batch job failed.")



