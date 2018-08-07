import signac
import flow
import environment
import shutil
import os


class RhacoProject(flow.FlowProject):
    def __init__(self, *args, **kwargs):
        super(RhacoProject, self).__init__(*args, **kwargs)
        env = flow.get_environment()
        print("Environment selected as =", env)
        self.add_operation(
            name="generate",
            cmd=lambda job: "source generate.sh {}".format(job),
            pre=[RhacoProject.parent_job],
            post=[RhacoProject.generated],
        )
        self.add_operation(
            name="simulate",
            cmd=lambda job: "ccmlogin /u/eot/anderso4/projects/rhaco-flow/simulate.sh {}".format(job),
            # cmd=lambda job: "source simulate.sh {}".format(job),
            # cmd=lambda job: "python -u operations.py simulate {}".format(job),
            pre=[RhacoProject.generated],
            post=[RhacoProject.simulated],
        )

    @flow.staticlabel()
    def generated(job):
        if job.sp.job_type == "child":
            # Get current project
            project = signac.get_project("./")
            # Find all jobs with the same statepoint as the parent
            parent_jobs = project.find_jobs(job.sp.parent_statepoint)
            if len(parent_jobs) == 1:
                parent_job = parent_jobs.next()
            else:
                raise SystemError(
                    "Found {} parent jobs, instead of one. Check"
                    " the workspace for inconsistencies.".format(len(parent_jobs))
                )
            parent_completed = RhacoProject.generated(parent_job)
            if parent_completed:
                # Copy the generated morphology
                shutil.copyfile(
                    os.path.join(parent_job._wd, "output.hoomdxml"),
                    os.path.join(job._wd, "output.hoomdxml"),
                )
                # Also copy the generate stdout
                # shutil.copyfile(
                #    os.path.join(parent_job._wd, "generate_stdout.log"),
                #    os.path.join(job._wd, "generate_stdout.log"),
                # )
            else:
                return False
        return job.isfile("output.hoomdxml")

    @flow.staticlabel()
    def simulated(job):
        if job.sp.job_type == "parent":
            return True
        else:
            return job.isfile("output_final.gsd")

    @flow.staticlabel()
    def parent_job(job):
        if job.sp.job_type == "parent":
            return True
        elif job.sp.job_type == "child":
            return False


if __name__ == "__main__":
    RhacoProject().main()
