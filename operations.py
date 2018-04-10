import subprocess


def generate(job):
    with job:
        subprocess.Popen(["lynx-create-morph", "-gd", str(job.sp.gas_density), "-z", str(job.sp.z), "-d", str(job.sp.dims)])


def generate(job):
    with job:
        subprocess.Popen(["lynx-run-hoomd", "-r", str(1E2)])


if __name__ == "__main__":
    import flow
    flow.run()
