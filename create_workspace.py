import copy
import sys
import signac
import itertools
from collections import OrderedDict


def get_gen_parameters():
    parameters = OrderedDict()
    # Generate Parameters
    parameters["stoichiometry"] = ["{'Mo':1,'V':0.3,'Nb':0.15,'Te':0.15}"]
    parameters["dimensions"] = ["10x10x1"]
    parameters["template"] = ["M1UnitCell.pdb"]
    parameters["crystal_separation"] = [25.0]
    parameters["z_reactor_size"] = [10.0]
    parameters["reactant_composition"] = ["{'C2H6':1}"]
    parameters["reactant_density"] = [0.1]
    parameters["forcefield"] = ["FF_opls_uff"]
    parameters["job_type"] = ["parent"]
    return list(parameters.keys()), list(itertools.product(*parameters.values()))


def get_sim_parameters():
    parameters = OrderedDict()
    # Simulate Parameters
    parameters["temperature"] = [433]
    parameters["run_time"] = [1E3]
    parameters["timestep"] = [1E-3]
    parameters["tau"] = [1E-1]
    parameters["job_type"] = ["child"]
    return list(parameters.keys()), list(itertools.product(*parameters.values()))


if __name__ == "__main__":
    project = signac.init_project("FirstParSweep")
    gen_param_names, gen_param_combinations = get_gen_parameters()
    sim_param_names, sim_param_combinations = get_sim_parameters()
    # Create the generate jobs
    for gen_params in gen_param_combinations:
        parent_statepoint = dict(zip(gen_param_names, gen_params))
        parent_job = project.open_job(parent_statepoint)
        parent_job.init()
        for sim_params in sim_param_combinations:
            child_statepoint = copy.deepcopy(parent_statepoint)
            child_statepoint.update(dict(zip(sim_param_names, sim_params)))
            child_statepoint["parent_statepoint"] = parent_job.statepoint()
            project.open_job(child_statepoint).init()
    project.write_statepoints()
