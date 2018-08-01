import subprocess
import inspect
import project
import signac
import time as T


def generate(job):
    generate_flags = {
        "stoichiometry": "-s",
        "dimensions": "-d",
        "template": "-t",
        "crystal_separation": "-c",
        "z_reactor_size": "-z",
        "reactant_composition": "-rc",
        "reactant_num_mol": "-rn",
        "reactant_density": "-rd",
        "forcefield": "-f",
        "integrate_crystal": "-i",
        "signac": "-sig",
    }
    with job, open("generate_stdout.log", "w+") as generate_stdout:
        # Always run with the signac flag -sig
        job_command = ["rhaco-create-morph", "-sig"]
        for flag in job.sp.keys():
            if flag in generate_flags:
                job_command += [str(generate_flags[flag]), str(job.sp[flag])]
        print("Executing job command:", job_command)
        generate = subprocess.Popen(job_command, stdout=generate_stdout)
        generate.wait()


def simulate(job):
    simulate_flags = {"temperature": "-t", "run_time": "-r", "timestep": "-s"}
    with job, open("hoomd_stdout.log", "w+") as hoomd_stdout:
        job_command = ["rhaco-run-hoomd"]
        for flag in job.sp.keys():
            if flag in simulate_flags:
                job_command += [str(simulate_flags[flag]), str(job.sp[flag])]
        job_command += ["output.hoomdxml"]
        print("Executing job command:", job_command)
        simulate = subprocess.Popen(job_command, stdout=hoomd_stdout)
        simulate.wait()


if __name__ == "__main__":
    import flow

    flow.run()
