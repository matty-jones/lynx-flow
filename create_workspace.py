import copy
import sys
import signac
import itertools
from collections import OrderedDict


def get_gen_parameters():
    parameters = OrderedDict()
    # Generate Parameters
    parameters["dimensions"] = ["30x15x1"]
    parameters["template"] = ["corundum.pdb"]
    parameters["crystal_separation"] = [40.0]
    parameters["reactant_composition"] = ["{'Ag_5nm':1}"]
    parameters["reactant_num_mol"] = [2]
    parameters["reactant_position"] = ["[[-2.8, 0.0, 5.8], [-2.8, 0.0, 5.8]]"]
    parameters["forcefield"] = ["['corundum.eam.fs', 'FF_Ag_Sapphire.xml']"]
    parameters["crystal_x"] = ["0.4759"]
    parameters["crystal_y"] = ["0.4759"]
    parameters["crystal_z"] = ["1.299"]
    parameters["job_type"] = ["parent"]
    return list(parameters.keys()), list(itertools.product(*parameters.values()))


def get_sim_parameters():
    parameters = OrderedDict()
    # Simulate Parameters
    parameters["temperature"] = [1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600]
    parameters["run_time"] = [1E7]
    parameters["omit_lj"] = ["Ag-Ag"]
    parameters["energy_scale_unit"] = ["23.06"]
    parameters["job_type"] = ["child"]
    return list(parameters.keys()), list(itertools.product(*parameters.values()))


if __name__ == "__main__":
    project = signac.init_project("AgNP_Corundum_TempSweep")
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
