import sys
import signac
import itertools
from collections import OrderedDict


def get_parameters():
    parameters = OrderedDict()
    parameters['stoichiometry'] = ["{'Mo':1,'V':0.3,'Nb':0.15,'Te':0.15}"]
    parameters['dimensions'] = ['2x2x2', '3x3x3', '4x4x4', '5x5x5']
    parameters['template'] = ['templateM1.pdb']
    parameters['crystal_separation'] = [25.0]
    parameters['z_reactor_size'] = [10]#[10, 15, 20, 25]
    parameters['gas_composition'] = ["{'C2H6':1}"]
    parameters['gas_density'] = [0.1]#[0.1, 0.2, 0.3, 0.4, 0.5]
    parameters['forcefield'] = ['FF_opls_uff']
    return list(parameters.keys()), list(itertools.product(*parameters.values()))



if __name__ == "__main__":
    project = signac.init_project('FirstParSweep')

    param_names, param_combinations = get_parameters()

    for params in param_combinations:
        job = project.open_job(dict(zip(param_names, params))).init()

    project.write_statepoints()
