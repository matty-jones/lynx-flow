import copy
import sys
import signac
import itertools
import numpy as np
from collections import OrderedDict


def get_gen_parameters():
    parameters = OrderedDict()
    # Generate Parameters
    parameters["dimensions"] = ["16x16x1"]
    parameters["template"] = ["Ag_surface.pdb"]
    # parameters["crystal_separation"] = [40.0]
    parameters["reactant_composition"] = ["{'Ag_5nm':'pos'}"]
    # parameters["reactant_num_mol"] = [2]
    parameters["reactant_position"] = ["[[-2.8, 0.0, 4.3], [2.8, 0.0, 4.3]]"]
    parameters["forcefield"] = ["['Ag.eam.fs']"]
    parameters["crystal_x"] = ["0.7810"]
    parameters["crystal_y"] = ["0.7810"]
    parameters["crystal_z"] = ["0.6440"]
    parameters["job_type"] = ["parent"]
    return list(parameters.keys()), list(itertools.product(*parameters.values()))


def get_sim_parameters():
    parameters = OrderedDict()
    # Simulate Parameters
    parameters["temperature"] = [1000]
    parameters["run_time"] = [1E8]
    parameters["omit_lj"] = ["Ag-Ag"]
    parameters["tau"] = [5]
    parameters["energy_scale_unit"] = ["23.06"]
    parameters["job_type"] = ["child"]
    return list(parameters.keys()), list(itertools.product(*parameters.values()))


if __name__ == "__main__":
    project = signac.init_project("PureAg_Sintering_1E8_DiffTau")
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
