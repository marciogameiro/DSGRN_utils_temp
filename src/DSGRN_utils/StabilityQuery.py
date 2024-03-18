### StabilityQuery.py
### MIT LICENSE 2022 Marcio Gameiro

import DSGRN
import DSGRN_utils

def StabilityQuery(network, param_indices=None, level=3):
    parameter_graph = DSGRN.ParameterGraph(network)
    if param_indices == None:
        param_indices = range(parameter_graph.size())
    param_stability = {}
    for par_index in param_indices:
        parameter = parameter_graph.parameter(par_index)
        morse_graph = DSGRN_utils.ConleyMorseGraph(parameter, level=level)
        attractors = [v for v in morse_graph.vertices() if not morse_graph.adjacencies(v)]
        n_stable = len(attractors)
        if n_stable not in param_stability:
            param_stability[n_stable] = []
        param_stability[n_stable].append(par_index)
    return param_stability
