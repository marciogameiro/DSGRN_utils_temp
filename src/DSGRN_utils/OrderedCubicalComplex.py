# OrderedCubicalComplex.py  # 2020-09-12
# MIT LICENSE 2020 Ewerton R. Vieira

from pychomp2 import *
import operator
import functools
# from itertools import product

class OrderedCubicalComplex:  # Ordered Cubical Complex

    def buildingOcc(self):
        for a in range(self.cc.size() - self.cc.size(self.D)):
            self.occ_.update({a: self.cc.star({a}) - {a}})

    def isLess(self, s, t):  # given s and t return true if s<t
        if t in self.occ_.get(s):
            return True
        return False

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == self.cc.size() - self.cc.size(self.D) - 1:
            raise StopIteration
        self.index = self.index + 1
        return self.occ_.get(self.index)

    def __call__(self, i):
        return self.occ_.get(i)

    def size(self):
        return self.cc.size() - self.cc.size(self.D)

    # Given a cubical complex cc return a dictonary that represents an order on cc.

    def Edges(self):  # build the set of (s, t) where s<t
        for s in range(self.size()):
            for t in self.occ_.get(s):
                self.edges.add((s, t))
        return self.edges

    def __init__(self, cc):
        """ initialize complex, graph, supporting data structures """
        self.cc = cc
        self.D = self.cc.dimension()
        self.index = 0

        self.occ_ = {}
        self.buildingOcc()

        self.edges = set()  # base set to create the set of (s, t) where s<t
