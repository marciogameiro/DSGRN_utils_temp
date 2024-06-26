### ComputeMorseGraph.py
### MIT LICENSE 2023 Marcio Gameiro

import DSGRN
import pychomp
import DSGRN_utils

def ConleyMorseGraph(parameter, prune_grad=True, level=0):
    # Compute the multivalued map (state transition graph)
    stg = DSGRN_utils.RookRulesCubicalComplex(parameter, level=level)
    # Compute the flow graded complex
    (scc_dag, graded_complex) = pychomp.FlowGradedComplex(stg.complex(), stg.diagram())
    # Compute the connection matrix of the graded complex
    connection_matrix = pychomp.ConnectionMatrix(graded_complex)
    # Compute the Morse graph
    morse_graph = DSGRN_utils.MorseGraph(stg, scc_dag, graded_complex, connection_matrix, prune_grad=prune_grad)
    return morse_graph, stg, graded_complex
