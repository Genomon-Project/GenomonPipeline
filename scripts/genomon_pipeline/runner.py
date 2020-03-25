#! /usr/bin/env python

import os
import datetime
import subprocess
import stat

class Runner(object):
    def __init__(self, singularity_script, qsub_option, log_dir, max_task, retry_count):
        self.qsub_option = qsub_option
        self.retry_count = retry_count
        self.singularity_script = os.path.abspath(singularity_script)
        self.jobname = os.path.basename(self.singularity_script).replace(".sh", "").replace("singularity_", "")
        self.log_dir = os.path.abspath(log_dir)
        self.max_task = max_task
        
    def task_exec(self):
        pass

class Drmaa_runner(Runner):
    
    def task_exec(self):
            import drmaa
        
            s = drmaa.Session()
            s.initialize()
             
            jt = s.createJobTemplate()
            jt.jobName = self.jobname
            jt.outputPath = ':' + self.log_dir
            jt.errorPath = ':' + self.log_dir
            jt.nativeSpecification = self.qsub_option
            jt.remoteCommand = self.singularity_script
            os.chmod(self.singularity_script, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IRGRP | stat.S_IXGRP)

            returncode = 0
            returnflag = True
            if self.max_task == 0:
                for var in range(0, (self.retry_count+1)):
                    jobid = s.runJob(jt)
                    returncode = 0
                    returnflag = True
                    now = datetime.datetime.now()
                    date = now.strftime("%Y-%m-%d %H:%M:%S")
                    print ("Your job has been submitted with id: " + jobid + " at Date/Time: " + date)
                    retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                    now = datetime.datetime.now()
                    date = now.strftime("%Y-%m-%d %H:%M:%S")
                    print ("Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date)
                    returncode = retval.exitStatus
                    returnflag = retval.hasExited
                    if returncode == 0 and returnflag: break
                s.deleteJobTemplate(jt)
                s.exit()

            else:
                joblist = s.runBulkJobs(jt,1,self.max_task,1)
                all_jobids = []
                for var in range(0, (self.retry_count+1)):
                    if len(all_jobids) > 0:
                        joblist = all_jobids
                        all_jobids = []
                    returncode = 0
                    returnflag = True
                    now = datetime.datetime.now()
                    date = now.strftime("%Y-%m-%d %H:%M:%S")
                    print ('Your job has been submitted with id ' + str(joblist) + " at Date/Time: " + date)
                    s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
                    for curjob in joblist:
                        print ('Collecting job ' + curjob)
                        retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
                        now = datetime.datetime.now()
                        date = now.strftime("%Y-%m-%d %H:%M:%S")
                        print ("Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date)
                        
                        if retval.exitStatus != 0 or not retval.hasExited:
                            returncode = retval.exitStatus
                            returnflag = retval.hasExited
                            if var == self.retry_count: break
                            jobId_list = retval.jobId.split(".")
                            taskId = int(jobId_list[1])
                            all_jobids.extend(s.runBulkJobs(jt,taskId,taskId,1))
                       
                    if returncode == 0 and returnflag: break
                s.deleteJobTemplate(jt)
                s.exit()

            if returncode != 0 or not returnflag: 
                raise RuntimeError("Job: " + str(retval.jobId)  + ' failed at Date/Time: ' + date)

class Qsub_runner(Runner):
    def task_exec(self):
        qsub_commands = ['qsub', '-sync', 'yes', '-N', self.jobname]
        if self.max_task != 0:
            qsub_commands.extend(['-t', '1-'+str(self.max_task)+':1'])

        qsub_options = []
        if type(self.qsub_option) == type(""):
            qsub_options += self.qsub_option.split(' ')
        returncode = subprocess.call(qsub_commands + qsub_options + [self.singularity_script])

        if returncode != 0: 
            raise RuntimeError("The batch job failed.")

def main(args):
    import yaml
    conf = yaml.safe_load(open(args.conf))
    
    if conf["drmaa"]:
        runner = Drmaa_runner(args.script, conf["qsub_option"], conf["log_dir"], conf["max_task"], conf["retry_count"])
    else:
        runner = Qsub_runner(args.script, conf["qsub_option"], conf["log_dir"], conf["max_task"], conf["retry_count"])
        
    runner.task_exec()
