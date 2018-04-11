"""Configuration of the project enviroment.

The environments defined in this module can be auto-detected.
This helps to define environment specific behaviour in heterogenous
environments.
"""
import flow
import flow.environment
from flow.environment import format_timedelta
from flow.slurm import SlurmJob
from flow.manage import JobStatus,ClusterJob
from datetime import timedelta
from time import sleep
import getpass
import subprocess
import scheduler as epoxpy_flow_schedulers


def _fetch(user=None):
    def parse_status(s):
        s = s.strip()
        if s == 'PD':
            return JobStatus.queued
        elif s == 'R':
            return JobStatus.active
        elif s in ['CG', 'CD', 'CA', 'TO']:
            return JobStatus.inactive
        elif s in ['F', 'NF']:
            return JobStatus.error
        return JobStatus.registered
    if user is None:
        user = getpass.getuser()
    cmd = ['squeue', '-u', user, '-h', '-o "%2t %100j"']
    try:
        result = subprocess.check_output(cmd).decode()
    except subprocess.CalledProcessError as error:
        print('error', error)
        raise
    except FileNotFoundError:
        raise RuntimeError("Slurm not available.")
    lines = result.split('\n')
    for line in lines:
        if line:
            status, name = line.strip()[1:-1].split()
            yield ClusterJob(name, parse_status(status))


class SQueueSlurmScheduler(flow.slurm.SlurmScheduler):
    def jobs(self):
        self._prevent_dos()
        for job in _fetch(user=self.user):
            #print('job yielded by _fetch',job, job.status())
            yield job


class fryEnvironment(flow.environment.SlurmEnvironment):
    hostname_pattern = 'fry'
    cores_per_node = 16
    scheduler_type = SQueueSlurmScheduler

    @classmethod
    def is_present(cls):
        return super(fryEnvironment, cls).is_present()

    @classmethod
    def mpi_cmd(cls, cmd, np):
        return 'srun -np {np} {cmd}'.format(n=np, cmd=cmd)

    @classmethod
    def script(cls, _id, **kwargs):
        nn=1
        walltime = timedelta(days=14)
        js = super(fryEnvironment, cls).script(_id=_id, **kwargs)
        js.writeline('#!/bin/bash')
        js.writeline('#SBATCH --job-name={}'.format(_id))
        js.writeline('#SBATCH -N {}'.format(nn))
        js.writeline('#SBATCH -t {}'.format(format_timedelta(walltime)))
        js.writeline('#SBATCH -n 8')
        js.writeline('#SBATCH -p batch')
        js.writeline('#SBATCH --output={}.o'.format(_id))
        js.writeline('#SBATCH --mail-type=All')
        js.writeline('#SBATCH --mail-user=mattyjones@boisestate.edu')
        js.writeline('#SBATCH --gres=gpu:1')

        js.writeline('on-conda')
        js.writeline('source activate lynx')
        return js

    @classmethod
    def submit(cls, script, flags=None, *args, **kwargs):
        sleep(0.05)
        return super(fryEnvironment, cls).submit(script, flags, *args,
                                                 **kwargs)


class r2Environment(flow.environment.SlurmEnvironment):
    hostname_pattern = 'r2'
    cores_per_node = 16
    scheduler_type = SQueueSlurmScheduler

    @classmethod
    def is_present(cls):
        return super(r2Environment, cls).is_present()

    @classmethod
    def mpi_cmd(cls, cmd, np):
        return 'srun -np {np} {cmd}'.format(n=np, cmd=cmd)

    @classmethod
    def script(cls, _id, **kwargs):
        nn=1
        walltime = timedelta(days=3)
        js = super(r2Environment, cls).script(_id=_id, **kwargs)
        js.writeline('#!/bin/bash')
        js.writeline('#SBATCH --job-name={}'.format(_id))
        #js.writeline('#SBATCH -N {}'.format(nn))
        js.writeline('#SBATCH -n 8')
        js.writeline('#SBATCH -t {}'.format(format_timedelta(walltime)))
        js.writeline('#SBATCH -p gpuq')
        js.writeline('#SBATCH --output={}.o'.format(_id))
        js.writeline('#SBATCH --mail-type=All')
        js.writeline('#SBATCH --mail-user=mattyjones@boisestate.edu')
        js.writeline('#SBATCH --gres=gpu:1')

        js.writeline('on-conda')
        js.writeline('source activate lynx')
        return js

    @classmethod
    def submit(cls, script, flags=None, *args, **kwargs):
        sleep(0.5)
        return super(r2Environment, cls).submit(script, flags, *args, **kwargs)


class kestrelEnvironment(flow.environment.SlurmEnvironment):
    hostname_pattern = 'kestrel'
    cores_per_node = 16
    scheduler_type = SQueueSlurmScheduler

    @classmethod
    def is_present(cls):
        return super(kestrelEnvironment, cls).is_present()

    @classmethod
    def mpi_cmd(cls, cmd, np):
        return 'srun -np {np} {cmd}'.format(n=np, cmd=cmd)

    @classmethod
    def script(cls, _id, **kwargs):
        nn=1
        walltime = timedelta(days=2)
        js = super(kestrelEnvironment, cls).script(_id=_id, **kwargs)
        js.writeline('#!/bin/bash')
        js.writeline('#SBATCH --job-name={}'.format(_id))
        js.writeline('#SBATCH -N {}'.format(nn))
        js.writeline('#SBATCH -n 8')
        js.writeline('#SBATCH -t {}'.format(format_timedelta(walltime)))
        js.writeline('#SBATCH -p batch')
        js.writeline('#SBATCH --output={}.o'.format(_id))
        js.writeline('#SBATCH --mail-type=All')
        js.writeline('#SBATCH --mail-user=mattyjones@boisestate.edu')
        js.writeline('#SBATCH --gres=gpu:1')

        js.writeline('on-conda')
        js.writeline('source activate lynx')

        return js

    @classmethod
    def submit(cls, script, flags=None, *args, **kwargs):
        sleep(0.5)
        return super(kestrelEnvironment, cls).submit(script, flags, *args,
                                                     **kwargs)
