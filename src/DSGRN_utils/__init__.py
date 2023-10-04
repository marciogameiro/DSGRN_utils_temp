# __init__.py
# Marcio Gameiro
# 2023-05-25
# MIT LICENSE

import DSGRN
import pychomp2 as pychomp

from DSGRN_utils.BlowupGraph import *
from DSGRN_utils.DrawPosetGraph import *
from DSGRN_utils.Poset_E import *
from DSGRN_utils.SetCM import *
from DSGRN_utils.SaveDatabaseJSON import *

def ConleyMorseGraph(parameter, parameter_graph, level=0):
    std = BlowupGraph(parameter, parameter_graph, level=level)
    (dag, fibration) = pychomp.FlowGradedComplex(std.complex(), std.diagram())
    connection_matrix = pychomp.ConnectionMatrix(fibration)
    conleyindices = connection_matrix.count()
    fringenode = fibration.value(std.complex().size() - 1)

    def non_trivial_scc(v):
        scc_v = [c for c in std.digraph.vertices() if fibration.value(c) == v]
        return len(scc_v) > 1 or any(c in std.digraph.adjacencies(c) for c in scc_v)

    CMG = InducedPoset_E(dag, lambda v : non_trivial_scc(v) and v != fringenode)
    return connection_matrix, CMG

def PlotMorseGraph(morse_graph):
    connection_matrix, CMG = morse_graph
    G = DrawPosetGraph(connection_matrix, CMG).graphviz()
    return graphviz.Source(G)
