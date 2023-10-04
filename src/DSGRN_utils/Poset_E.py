# Plot2D.py  # 2020-09-12
# Base on Shaun Poset.py
# MIT LICENSE 2020 Ewerton R. Vieira

from pychomp2.DirectedAcyclicGraph import *

import graphviz

class Poset_E:  # Represent a Poset

    def ListOfNivelSets(self):  # return the list of increasing Nivel sets of a given poset

        P0 = {v for v in self.vertices() if len(self.children(v)) == 0}

        def NivelSet(P_i):  # return the parent set of a given set
            P_iplus = set()
            for a in P_i:
                #print(str(a)+' -> '+ str(CMG.parents(a)))
                P_iplus = P_iplus.union(self.parents(a))
            return P_iplus - P_i

        P = [P0]
        j = 0
        while P[j] != set():
            P.append(NivelSet(P[j]))
            j += 1

        return P

    def __init__(self, graph):
        """
        Create a Poset P from a DAG G such that x <= y in P iff there is a path from x to y in G
        """
        self.vertices_ = set(graph.vertices())
        self.descendants_ = graph.transitive_closure()
        self.ancestors_ = self.descendants_.transpose()
        self.children_ = graph.transitive_reduction()
        self.parents_ = self.children_.transpose()

    def __iter__(self):
        """
        Allows for the semantics
          [v for v in poset]
        """
        return iter(self.vertices())

    def vertices(self):
        """
        Return the set of elements in the poset
        """
        return self.vertices_

    def parents(self, v):
        """
        Return the immediate predecessors of v in the poset
        """
        return self.parents_.adjacencies(v)

    def children(self, v):
        """
        Return the immediate successors of v in the poset
        """
        return self.children_.adjacencies(v)

    def ancestors(self, v):
        """
        Return the set { u : u > v }
        """
        return self.ancestors_.adjacencies(v)

    def upset(self, v):
        """
        Return the set { u : u >= v }
        """
        return self.ancestors(v).union({v})

    def descendants(self, v):
        """
        Return the set { u : u < v }
        """
        return self.descendants_.adjacencies(v)

    def downset(self, v):
        """
        Return the set { u : u <= v }
        """
        return self.descendants(v).union({v})

    def interval(self, p, q):
        """
        Return the minimal interval that has p and q
        """
        if p < q:
            return self.downset(q).intersection(self.upset(p))
        else:
            return self.downset(p).intersection(self.upset(q))

    def isConvexSet(self, I):  # return true if I is a convex subset in poset P
        for a in I:
            for b in I:
                if a < b:
                    if not self.interval(a, b).issubset(set(I)):
                        return False
        return True

    def less(self, u, v):
        """
        Return True if u < v, False otherwise
        """
        return u in self.ancestors(v)

    def maximal(self, subset):
        """
        Return the set of elements in "subset" which are maximal
        """
        return frozenset({u for u in subset if not any(self.less(u, v) for v in subset)})

    def _repr_svg_(self):
        """
        Return svg representation for visual display
        """
        return graphviz.Source(self.children_.graphviz())._repr_svg_()

#from pychomp.DirectedAcyclicGraph import *

def InducedPoset_E(G, predicate):
    result = DirectedAcyclicGraph()
    S = set([v for v in G.vertices() if predicate(v)])
    for v in S:
        result.add_vertex(v)
    for v in S:
        for u in G.descendants(v):
            if u in S and u != v:
                result.add_edge(v, u)
    return Poset_E(result)
