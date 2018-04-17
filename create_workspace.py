import sys
import signac
import itertools
import copy
from collections import OrderedDict


def get_gen_parameters():
    parameters = OrderedDict()
    # Generate Parameters
    parameters['stoichiometry'] = ["{'Mo':1,'V':0.3,'Nb':0.15,'Te':0.15}"]
    parameters['dimensions'] = ['10x10x1']
    parameters['template'] = ['M1UnitCell.pdb']
    parameters['crystal_separation'] = [25.0]
    parameters['z_reactor_size'] = [10, 15]
    parameters['gas_composition'] = ["{'C2H6':1}"]
    parameters['gas_density'] = [0.001356]
    parameters['forcefield'] = ['FF_opls_uff']
    parameters['stage'] = ['parent']
    return list(parameters.keys()),list(itertools.product(
        *parameters.values()))


def get_sim_parameters():
    parameters = OrderedDict()
    # Simulate Parameters
    parameters['temperature'] = [633, 733]
    parameters['run_time'] = [1E6]
    parameters['timestep'] = [1E-3]
    parameters['stage'] = ['child']
    return list(parameters.keys()),list(itertools.product(
        *parameters.values()))


if __name__ == "__main__":
    project = signac.init_project('FirstParSweep')
    gen_param_names, gen_param_combinations = get_gen_parameters()
    sim_param_names, sim_param_combinations = get_sim_parameters()
    # Create the generate jobs
    for gen_params in gen_param_combinations:
        parent_statepoint = dict(zip(gen_param_names, gen_params))
        project.open_job(parent_statepoint).init()
        for sim_params in sim_param_combinations:
            child_statepoint = copy.deepcopy(parent_statepoint)
            child_statepoint.update(dict(zip(sim_param_names, sim_params)))
            child_statepoint['parent_statepoint'] = parent_statepoint
            project.open_job(child_statepoint).init()
    project.write_statepoints()
