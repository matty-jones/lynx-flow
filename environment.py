"""Configuration of the project enviroment.

The environments defined in this module can be auto-detected.
This helps to define environment specific behaviour in heterogenous
environments.
"""
import flow
import flow.environment
from flow.environment import format_timedelta
from flow.slurm import SlurmJob
from flow.manage import JobStatus, ClusterJob
from datetime import timedelta
from time import sleep
import getpass
import subprocess
import scheduler as rhaco_flow_schedulers


def _fetch(user=None):
    def parse_status(s):
        s = s.strip()
        if s == "PD":
            return JobStatus.queued
        elif s == "R":
            return JobStatus.active
        elif s in ["CG", "CD", "CA", "TO"]:
            return JobStatus.inactive
        elif s in ["F", "NF"]:
            return JobStatus.error
        return JobStatus.registered

    if user is None:
        user = getpass.getuser()
    cmd = ["squeue", "-u", user, "-h", '-o "%2t %100j"']
    try:
        result = subprocess.check_output(cmd).decode()
    except subprocess.CalledProcessError as error:
        print("error", error)
        raise
    except FileNotFoundError:
        raise RuntimeError("Slurm not available.")
    lines = result.split("\n")
    for line in lines:
        if line:
            status, name = line.strip()[1:-1].split()
            yield ClusterJob(name, parse_status(status))


class SQueueSlurmScheduler(flow.slurm.SlurmScheduler):
    def jobs(self):
        self._prevent_dos()
        for job in _fetch(user=self.user):
            # print('job yielded by _fetch',job, job.status())
            yield job


class fryEnvironment(flow.environment.SlurmEnvironment):
    hostname_pattern = "fry"
    cores_per_node = 16
    scheduler_type = SQueueSlurmScheduler

    @classmethod
    def is_present(cls):
        return super(fryEnvironment, cls).is_present()

    @classmethod
    def mpi_cmd(cls, cmd, np):
        return "srun -np {np} {cmd}".format(n=np, cmd=cmd)

    @classmethod
    def script(cls, _id, **kwargs):
        nn = 1
        walltime = timedelta(hours=1)
        js = super(fryEnvironment, cls).script(_id=_id, **kwargs)
        js.writeline("#!/bin/bash")
        js.writeline("#SBATCH --job-name={}".format(_id))
        js.writeline("#SBATCH -N {}".format(nn))
        js.writeline("#SBATCH -t {}".format(format_timedelta(walltime)))
        js.writeline("#SBATCH -n 8")
        js.writeline("#SBATCH -p batch")
        js.writeline("#SBATCH --output={}.o".format(_id))
        js.writeline("#SBATCH --mail-type=All")
        js.writeline("#SBATCH --mail-user=mattyjones@boisestate.edu")
        js.writeline("#SBATCH --gres=gpu:1")

        js.writeline("on-conda")
        js.writeline("source activate rhaco")
        return js

    @classmethod
    def submit(cls, script, flags=None, *args, **kwargs):
        sleep(0.05)
        return super(fryEnvironment, cls).submit(script, flags, *args, **kwargs)


class r2Environment(flow.environment.SlurmEnvironment):
    hostname_pattern = "r2"
    cores_per_node = 16
    scheduler_type = SQueueSlurmScheduler

    @classmethod
    def is_present(cls):
        return super(r2Environment, cls).is_present()

    @classmethod
    def mpi_cmd(cls, cmd, np):
        return "srun -np {np} {cmd}".format(n=np, cmd=cmd)

    @classmethod
    def script(cls, _id, **kwargs):
        nn = 1
        walltime = timedelta(hours=1)
        js = super(r2Environment, cls).script(_id=_id, **kwargs)
        js.writeline("#!/bin/bash")
        js.writeline("#SBATCH --job-name={}".format(_id))
        # js.writeline('#SBATCH -N {}'.format(nn))
        js.writeline("#SBATCH -n 8")
        js.writeline("#SBATCH -t {}".format(format_timedelta(walltime)))
        js.writeline("#SBATCH -p gpuq")
        js.writeline("#SBATCH --output={}.o".format(_id))
        js.writeline("#SBATCH --mail-type=All")
        js.writeline("#SBATCH --mail-user=mattyjones@boisestate.edu")
        js.writeline("#SBATCH --gres=gpu:1")

        js.writeline("on-conda")
        js.writeline("source activate rhaco")
        return js

    @classmethod
    def submit(cls, script, flags=None, *args, **kwargs):
        sleep(0.5)
        return super(r2Environment, cls).submit(script, flags, *args, **kwargs)


class kestrelEnvironment(flow.environment.SlurmEnvironment):
    hostname_pattern = "kestrel"
    cores_per_node = 16
    scheduler_type = SQueueSlurmScheduler

    @classmethod
    def is_present(cls):
        return super(kestrelEnvironment, cls).is_present()

    @classmethod
    def mpi_cmd(cls, cmd, np):
        return "srun -np {np} {cmd}".format(n=np, cmd=cmd)

    @classmethod
    def script(cls, _id, **kwargs):
        nn = 1
        walltime = timedelta(hours=1)
        js = super(kestrelEnvironment, cls).script(_id=_id, **kwargs)
        js.writeline("#!/bin/bash")
        js.writeline("#SBATCH --job-name={}".format(_id))
        js.writeline("#SBATCH -N {}".format(nn))
        js.writeline("#SBATCH -n 8")
        js.writeline("#SBATCH -t {}".format(format_timedelta(walltime)))
        js.writeline("#SBATCH -p batch")
        js.writeline("#SBATCH --output={}.o".format(_id))
        js.writeline("#SBATCH --mail-type=All")
        js.writeline("#SBATCH --mail-user=mattyjones@boisestate.edu")
        js.writeline("#SBATCH --gres=gpu:1")

        js.writeline("on-conda")
        js.writeline("source activate rhaco")

        return js

    @classmethod
    def submit(cls, script, flags=None, *args, **kwargs):
        sleep(0.5)
        return super(kestrelEnvironment, cls).submit(script, flags, *args, **kwargs)


class blueWatersEnvironment(flow.environment.DefaultTorqueEnvironment):
    hostname_pattern = 'h2ologin*'  # TODO: run python -c "import socket; print(socket.gethostname())"
    cores_per_node = 1
    #scheduler_type = rhaco_flow_schedulers.PBSProScheduler

    @classmethod
    def is_present(cls):
        return super(blueWatersEnvironment, cls).is_present()

    @classmethod
    def script(cls, _id, **kwargs):
        nn=1
        kwargs['nn']=None
        js = super(blueWatersEnvironment, cls).script(_id=_id, **kwargs)
        js.writeline('#!/bin/bash')
        js.writeline("#PBS -j {}".format(_id))
        js.writeline('#PBS -l walltime=00:10:00:nodes=1:ppn=16:xk')
        js.writeline('#PBS -q debug')
        js.writeline('#PBS -o {}.o'.format(_id))

        js.writeline('module unload PrgEnv-cray')
        js.writeline('module load PrgEnv-gnu')
        js.writeline('module load bwpy/1.1.0')
        js.writeline('module load cudatoolkit/7.5.18-1.0502.10743.2.1')
        js.writeline('module load ccm')
        js.writeline('bwpy-environ')
        js.writeline('source /u/eot/anderso4/projects/rhaco-virtualenv/bin/activate') 
        js.writeline('export PYTHONPATH="/u/eot/anderso4/software/build-hoomd-on-blue-waters/hoomd_blue/build"')
        return js

    @classmethod
    def submit(cls, script, flags=None, *args, **kwargs):
        sleep(0.5)
        return super(blueWatersEnvironment, cls).submit(script, flags, *args, **kwargs)
