### ComputeMorseGraph.py
### MIT LICENSE 2023 Marcio Gameiro

import DSGRN
import pychomp
import DSGRN_utils

def ConleyMorseGraph(parameter, level=0):
    # Compute the multivalued map (state transition graph)
    stg = DSGRN_utils.BlowupGraph(parameter, level=level)
    # Compute the flow graded complex
    (scc_dag, graded_complex) = pychomp.FlowGradedComplex(stg.complex(), stg.diagram())
    # Compute the connection matrix of the graded complex
    connection_matrix = pychomp.ConnectionMatrix(graded_complex)
    # Compute the Morse graph
    morse_graph = DSGRN_utils.MorseGraph(stg, scc_dag, graded_complex, connection_matrix)
    return morse_graph