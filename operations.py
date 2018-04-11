import subprocess
import inspect


def pipeline(job):
    print("Running full pipeline for", job)
    generate(job)
    simulate(job)


def generate(job):
    generate_flags = {'stoichiometry': '-s', 'dimensions': '-d',
                      'template': '-t', 'crystal_separation': '-c',
                      'z_reactor_size': '-z', 'gas_composition': '-gc',
                      'gas_num_mol': '-gn', 'gas_density': '-gd',
                      'forcefield': '-f', 'integrate_crystal': '-i',
                      'signac': '-sig'}
    with job:
        # Always run with the signac flag -sig
        job_command = ["lynx-create-morph", "-sig"]
        for flag in job.sp.keys():
            if flag in generate_flags:
                job_command += [str(generate_flags[flag]), str(job.sp[flag])]
        generate = subprocess.Popen(job_command)
        generate.wait()


def simulate(job):
    simulate_flags = {'temperature': '-t', 'run_time': '-r', 'timestep': '-s'}
    with job:
        job_command = ["lynx-run-hoomd"]
        for flag in job.sp.keys():
            if flag in simulate_flags:
                job_command += [str(simulate_flags[flag]), str(job.sp[flag])]
        job_command += ["output.hoomdxml"]
        simulate = subprocess.Popen(job_command)
        simulate.wait()


if __name__ == "__main__":
    import flow
    flow.run()
