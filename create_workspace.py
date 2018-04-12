import sys
import signac
import itertools
from collections import OrderedDict


def get_parameters():
    parameters = OrderedDict()
    # Generate Parameters
    parameters['stoichiometry'] = ["{'Mo':1,'V':0.3,'Nb':0.15,'Te':0.15}"]
    parameters['dimensions'] = ['10x10x1', '10x10x2', '10x10x3']
    parameters['template'] = ['M1UnitCell.pdb']
    parameters['crystal_separation'] = [25.0]
    parameters['z_reactor_size'] = [10, 15, 20, 25]
    parameters['gas_composition'] = ["{'C2H6':1}"]
    parameters['gas_density'] = [0.001356]
    parameters['forcefield'] = ['FF_opls_uff']
    # Simulate Parameters
    parameters['temperature'] = [233, 333, 433, 533, 633, 733, 833]
    parameters['run_time'] = [1E6]
    parameters['timestep'] = [1E-3]
    return list(parameters.keys()),list(itertools.product(
        *parameters.values()))


if __name__ == "__main__":
    project = signac.init_project('FirstParSweep')
    param_names, param_combinations = get_parameters()
    for params in param_combinations:
        job = project.open_job(dict(zip(param_names, params))).init()
    project.write_statepoints()
