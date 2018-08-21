from flow.manage import Scheduler
from flow.torque import TorqueJob
import tempfile
import subprocess
from flow.manage import ClusterJob, JobStatus
import signac


class DummyTorqueJob(ClusterJob):
    def __init__(self, job):
        self.job = job

    def _id(self):
        return self.job.get_id()

    def __str__(self):
        return str(self._id())

    def name(self):
        return self.job.get_id()

    def status(self):
        return JobStatus.registered


class PBSProScheduler(Scheduler):
    submit_cmd = ['qsub']

    def __init__(self, user=None, **kwargs):
        super(PBSProScheduler, self).__init__(**kwargs)
        self.user = user

    def jobs(self):
        self._prevent_dos()
        project = signac.get_project('.')
        for job in project:
            yield DummyTorqueJob(job)

    def submit(self, script, after=None, pretend=False, hold=False, flags=None, *args, **kwargs):
        if flags is None:
            flags = []
        elif isinstance(flags, str):
            flags = flags.split()

        submit_cmd = self.submit_cmd + flags

        if after is not None:
            submit_cmd.extend(
                ['-W', 'depend="afterok:{}"'.format(after.split('.')[0])])

        if hold:
            submit_cmd += ['-h']

        if pretend:
            print("# Submit command: {}".format(' '.join(submit_cmd)))
            print(script.read())
            print()
        else:
            with tempfile.NamedTemporaryFile() as tmp_submit_script:
                tmp_submit_script.write(script.read().encode('utf-8'))
                tmp_submit_script.flush()
                output = subprocess.check_output(
                    submit_cmd + [tmp_submit_script.name])
            jobsid = output.decode('utf-8').strip()
            return jobsid

    @classmethod
    def is_present(cls):
        try:
            subprocess.check_output(['qsub', '--version'], stderr=subprocess.STDOUT)
        except (IOError, OSError):
            return False
        else:
            return True
