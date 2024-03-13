# RookRulesCubicalComplex.py  # 2021-12-01
# MIT LICENSE 2021 Bernardo Rivas, Ewerton R. Vieira

import DSGRN
import re, time
from pychomp import *
from itertools import product
from itertools import combinations
from DSGRN._dsgrn import *
from DSGRN_utils.OrderedCubicalComplex import *


class BlowupGraph:
    def complex(self):
        return self.fc

    def get_vertices(self):
        return self.digraph.vertices_

    def get_adjacencies(self):
        return self.digraph.adjacency_lists_

    def diagram(self):
        return self.digraph.adjacencies

    def cc_to_fc(self, a):  # given a cell in cc return the correspond topcell in fc (blowup/face cubical complex)
        return self.D_cc_to_fc.get(a)

    def fc_to_cc(self, a):  # given a cell in cc return the correspond topcell in fc (blowup/face cubical complex)
        return self.D_fc_to_cc.get(a)

    def get_rmap(self,s): # given a cell in cc return the corresponding regulation map
        # coord = [ 1 if n == 0 else self.D-1 if n == self.D else n for n in self.cc.coordinates(a) ]
        coord = self.cc.coordinates(s)
        limit =  [ len(self.network.outputs(d))+1 for d in range(self.D) ]

        pos = [ coord[i] if (coord[i]>0 and coord[i]<limit[i]) else int(coord[i]*(limit[i]-2)/limit[i]+1) for i in range(self.D)]

        return { n : self.rmap[n][pos[n]-1] for n in self.InessentialDirections(s) }

    def Dict_cc_fc(self):
        """ create two dictonary {a:b} and {b:a} where  a is a cell in cc and b is the correspond topcell in fc (blowup/face cubical complex)"""
        D_ab = {}  # Dictonary cc to fc
        D_ba = {}  # Dictonary fc to cc
        B = self.fc.boxes()
        pv = [1]

        def prodd(j):
            while j > 0:
                return B[j] * prodd(j - 1)
            return B[0]

        # Compute PV := [1, boxes[0], boxes[0]*boxes[1], ..., boxes[0]*boxes[1]*...*boxes[D-1]]
        for i in range(len(B) - 1):  # Maybe range(len(B)) check for bugs
            pv.append(prodd(i))

        initial = self.fc.size() - self.fc.size(self.D)  # count from the top dimensional cells

        for a in range(self.cc.size()):
            baryc_a = self.cc.barycenter(a)
            b = initial
            for ii in range(len(pv)):  # convert coordenate to value of a cell (begin from initial (top cells))
                b += baryc_a[ii] * pv[ii]

            D_ab.update({a: b})
            D_ba.update({b: a})

        return D_ab, D_ba

    def build_rmap(self):
        """
            Returns a list with the order of thresholds for each variable.
            Example:
                partial order: T[X->X] < T[X->Y], T[Y->Y] < T[Y->X]
                output: [[1,2],[2,1]]
            Can be optimized.
        """
        partialorders = self.parameter.partialorders('T')
        thresholds = re.findall(r'\[.*?\]', partialorders)
        for c in '[->]':
            thresholds = [s.replace(c,"") for s in thresholds]
        output = []
        for i in range(self.network.size()):
            order = []
            S = self.network.name(i)
            for word in thresholds:
                for j in range(self.network.size()):
                    T = self.network.name(j)
                    if word == S+T:
                        order.append(j)
            output.append(order)
        return output

    def closedrightfringe(self, i):
        return any(self.cc.rightfringe(k) for k in self.cc.star({i}))

    def blowupinfinity(self, i):
        """
        determine if cell i (in cubical complex) is in the fringe,
        i. e., associate with the repeller at infinity
        """
        return self.closedrightfringe(i)

    def blowupfringe(self, i):
        """
        determine if cell i (in blowup complex) is in the fringe
        """
        return self.blowupinfinity(self.fc_to_cc(i))

    def InessentialDirections(self, s):  # i.e. inessential directions
        """
        return the set of normal variables for a cubical cell s
        """
        shape = self.cc.cell_shape(s)
        return [d for d in range(self.D) if shape & (1 << d) == 0]

    def EssentialDirections(self, s):  # i.e. essential directions
        """
        return the set of tangent variables for a cubical cell s
        """
        shape = self.cc.cell_shape(s)
        return [d for d in range(self.D) if shape & (1 << d) != 0]

    def shapevector(self,s):
        """
        return the shape vector for a cubical cell s
        """
        return [1 if n in self.EssentialDirections(s) else 0 for n in range(self.D)]

    def RelativePositionVector(self, s, t):
        """
        Return the relative position vector of two cubical cells;
        this is equal to the difference of combinatorial position of two cells,
        regarding vertices to be at (even,even,...) integer coordinates and
        higher dimensional cells to have odd coordinates for their dimensions with
        extent

        Disclaimer: Relative position p(xi,mu) points from xi to mu.\nIn the paper it always points from the larger to the smaller
        """
        return [y-x for (x, y) in zip(self.cc.barycenter(s), self.cc.barycenter(t))]


    def Absorbing(self, cctopcell, d, direction):
        """
        Return True if the codim 1 face of cctopcell
          collapsed in dimension d in direction "direction"
          is absorbing. Here direction = -1 means left
          wall and direction = 1 means right wall.
        """
        # Absorbing means the flow goes towards the face.
        # The cc topcell indexing is not compatible with labelling indexing
        # due to the extra wrapped layer, so we compute the index idx
        coords = self.cc.coordinates(cctopcell)
        # "Wrap layer" boxes have all walls absorbing.
        if any(coords[i] == self.limits[i] for i in range(self.D)):
            return True
        idx = sum(c * self.pv[i] for (i, c) in enumerate(coords))
        labelbit = 1 << (d + (self.D if direction == 1 else 0))  # bug fix
        # labelbit = 1 << (d + (self.D if direction == -1 else 0))
        return self.labelling[idx] & labelbit != 0

    def flowdirection(self, cctopcell, s, t):
        """
        Given a cubical domain, and a dual-order wall (s,t) (i.e. s is a face
        of t in cc) return 1 if a proof of transversality from s->t is known
        valid in dom and return -1 if a proof of transversality from t->s is
        known valid in dom. It return 0 otherwise.
        """
        sshape = self.cc.cell_shape(s)
        scoords = self.cc.coordinates(s)
        tshape = self.cc.cell_shape(t)
        tcoords = self.cc.coordinates(t)
        xorshape = sshape ^ tshape  # 1's in positions where s is 0 and t is 1.
        absorbing = False
        repelling = False
        for d in range(self.D):
            if xorshape & (1 << d):
                direction = 1 if scoords[d] > tcoords[d] else -1  # -1 left, 1 right ## bug fix
                # direction = -1 if scoords[d] > tcoords[d] else 1  # -1 left, 1 right
                if self.Absorbing(cctopcell, d, direction):
                    absorbing = True
                else:
                    repelling = True
        if absorbing and not repelling:
            return -1
        elif repelling and not absorbing:
            return 1
        else:
            return 0

    def DoubleArrow(self,S,T):
        """
        Return True if (S,T) and (T,S) are edges in DirectedAcyclicGraph
        """
        return S in self.digraph.adjacencies(T) and T in self.digraph.adjacencies(S)

    def Arrow(self,S,T):
        """
        Return true if (S,T) edge is in DirectedAcyclicGraph
        """
        return T in self.digraph.adjacencies(S)

    def RookField(self, kappa):
        """Return the rook field of a cubical top-cell kappa as a list
        of tuples, where each tuple gives the sign of the flow at the
        lef and right walls, respectively.
        """
        sign = lambda kappa, d, direc: direc if self.Absorbing(kappa, d, direc) else -direc
        return [(sign(kappa, d, -1), sign(kappa, d, 1)) for d in range(self.D)]

    def Phi(self, sigma, kappa):
        """
        Return the rook field extension from a top-cell kappa to a lower dimensionl cell sigma.
        """
        rf = self.RookField(kappa)        # Rook field of top cell
        J_i = self.InessentialDirections(sigma) # Inessential directions
        p = self.RelativePositionVector(sigma, kappa)
        # Set flow to zero if essential direction and opposite signs at walls
        phi = lambda left, right, d: (left if p[d] == 1 else right) if d in J_i else 0
        return [left if left == right else phi(left, right, d) for d, (left, right) in enumerate(rf)]

    def Psi(self, sigma, tau, kappa):
        """
        Return the rook field extension from sigma to tau with respect to kappa
        """
        PsiSigmaKappa = self.Phi(sigma,kappa)
        PsiTauKappa = self.Phi(tau, kappa)
        return [ PsiSigmaKappa[n] if n in self.EssentialDirections(tau) else PsiTauKappa[n] for n in range(self.D)]

    def Psi_n(self,sigma,n):
        """
        Return the n-th entry of the rook field at sigma.
        Well defined if n is in G(sigma) or sigma top-cell.
        """
        kappa = self.cc.topstar(sigma)
        kappa = kappa[0]
        PsiSigmaKappa = self.Phi(sigma, kappa)
        return PsiSigmaKappa[n]

    def isOpaque(self, s):
        """
        Return True if sigma is an opaque cubical cell.
        """
        if self.closedrightfringe(s):
            return False
        if self.InessentialDirections(s) == []:
            return True
        rmap = self.RegulationMap(s)
        if rmap is None or rmap == {}:
            return False
        else:
            return set(rmap.keys())==set(rmap.values())

    def OpaqueCells(self):
        """
        Return a list of opaque cells of the cubical complex from analysis on the digraph vertices
        """
        return [self.fc_to_cc(V) for V in self.digraph.vertices() if (self.isOpaque(self.fc_to_cc(V)) and (not self.closedrightfringe(self.fc_to_cc(V))))]

    def OpaqueDirections(self,s):
        """
        Return a list of opaque directions of sigma
        """
        rmap = self.RegulationMap(s)
        gradient = self.GradientDirections(s)
        inessential = self.InessentialDirections(s)
        if inessential != []:
            return []
        else:
            return [ n for n in inessential if (n in set(rmap.values()) and n not in gradient) ]

    def GradientDirections(self,sigma):
        """
        Return a list of gradient directions of sigma.
        """
        phi_tstar = [self.Phi(sigma, kappa) for kappa in self.cc.topstar(sigma) if not self.cc.rightfringe(sigma)]
        return [ n for (n,phi) in enumerate(zip(*phi_tstar)) if (set(phi) == {-1} or set(phi) == {1}) ]

    def NeutralDirections(self,sigma):
        """
        Return a list of neutral directions of sigma
        """
        phi_tstar = [self.Phi(sigma, kappa) for kappa in self.cc.topstar(sigma) if not self.cc.rightfringe(sigma)]
        return  [ n for (n,phi) in enumerate(zip(*phi_tstar)) if (0 in set(phi)) ]

    def isEquilibriumCell(self, sigma):
        """
        Return True if sigma is an equilibrium cubical cell.
            def equilibrium cell: - opaque cell
                                  - no gradient direction
        """
        if self.closedrightfringe(sigma):
            return False
        if not self.isOpaque(sigma):
            return False
        if set(self.GradientDirections(sigma)) != set():
            return False
        return True

    def JiCapJe(self,s,t):
        """
        Return J_i(s) \cap J_e(t)
        """
        return [ n for n in self.EssentialDirections(t) if n in self.InessentialDirections(s) ]

    def isThereACycle(self,s,t):
        """
        Return True if there is n \in J_i(s) \cap J_e(t) such that o_s^k(n) = n.
        """
        for cycle in self.CycleDecomposition(self.RegulationMap(s)):
            if len(cycle) > 1 and any(n in cycle for n in self.JiCapJe(s,t)):
                return True
        return False

    def isEquilibriumNonTopCell(self, s):  # New, Marcio
        """
        Return True if s is an equilibrium non top cell
        """
        if self.cc.cell_dim(s) == self.D:
            return False
        return self.isEquilibriumCell(s)

    def EquilibriumCells(self):
        """
        Returning list of equilibrium cells
        """
        return [s for s in self.cc if (not self.closedrightfringe(s)) and self.isEquilibriumCell(s)]

    def EquilibriumNonTopCells(self):  # New, Marcio
        """
        Returning list of equilibrium non top cells
        """
        return [s for s in self.cc if (not self.closedrightfringe(s)) and self.isEquilibriumNonTopCell(s)]

    def Feedback(self,cycle):
        """
        return {  1   if cycle has an even number of repressors
               { -1   otherwise
        """
        edges = list(zip(cycle[:-1],cycle[1:]))+[(cycle[-1],cycle[0])]
        return functools.reduce(operator.mul, [ 1 if self.network.interaction(i,j) else -1 for (i,j) in edges])

    def CrossingNumber(self, q, cycle):
        edges = list(zip(cycle[:-1],cycle[1:]))+[(cycle[-1],cycle[0])]
        return functools.reduce(operator.add, [ 1 if ((1 if self.network.interaction(i,j) else -1) != q[i]*q[j]) else 0 for (i,j) in edges])

    def RegulationMap(self,s):
        """
        Return a dictionary with key-value pairs
        { (i,j) : s has inessential variable i which regulates variable j}
        It gets the local inducement map and checks if it actively regulates the variable
        """
        phi_tstar = [self.Phi(s, kappa) for kappa in self.cc.topstar(s) if not self.cc.rightfringe(s)]
        pos_tstar = [self.RelativePositionVector(s,kappa) for kappa in self.cc.topstar(s) if not self.cc.rightfringe(s)]

        phi = [ [x for x in v] for v in zip(*phi_tstar) ]
        pos = [ [x for x in v] for v in zip(*pos_tstar) ]

        rmap = self.get_rmap(s).copy()
        # The condition below checks if i is impacting rmap[i]. If not, we want to remove it.
        remove_list = [i for i in rmap if (sum([x*y for x,y in zip(phi[rmap[i]],pos[i])]) == 0) ]
        for key in remove_list:
            del rmap[key]
        return rmap

    def CycleDecomposition(self,perm):
        """
        Given a key-value map "perm" defining a permutation,
        return a list of cycles, where each cycle is represented
        as a list of elements
        e.g.  CycleDecomposition({ 1:1, 3:4, 4:3 }) may return
              [[1],[3,4]]  (or any equivalent answer)

        Added: It drops cycles candidates that aren't cycles
        e.g.  CycleDecomposition({ 1:2, 2:1, 3:1 }) may return
              [[1,2]]   (or any equivalent answer)
        """
        def ExtractCycle(G):
            cycle = list(G.popitem())
            while cycle[-1] in G:
                cycle.append(G.pop(cycle[-1]))
            if cycle[-1] == cycle[0]:
                return cycle[:-1]
            else:
                return []

        # This gets rid of values that are not actually cyclic.
        # It will run n times, where n is the longest almost cycle:
        # For example: a_0 -> a_1 -> ... -> a_n and ( a_n 0 1 2 ... s)
        def NewFormat(dict):
            output = {}
            for k,v in dict.items():
                if v in dict.keys():
                    output[v] = dict[v]
            return output

        rmap = perm.copy()
        new_rmap = NewFormat(rmap)
        while rmap != new_rmap:
            rmap = new_rmap
            new_rmap = NewFormat(rmap)
        cycles = []
        while rmap:
            cycles.append(ExtractCycle(rmap))
        # print("cycles =",cycles)
        return cycles

    def OpacityCycles(self,s):
        """
        Return list of regulatory cycles for an opaque cubical cell s if it exists
        Otherwise it returns []
        """
        rmap = self.RegulationMap(s)
        return [ cycle for cycle in self.CycleDecomposition(rmap) if cycle != [] ]

    def Rule0(self):
        """
        Make a digraph with the following vertices and edges
            V = { s : s in Cubical Complex }
            E = { (s,t) : s<t, t<s, s=t}
        then remove (s,t) when t in fringe.
        """
        for S in self.fc(self.D):
            if not self.blowupfringe(S):
                self.digraph.add_edge(S, S)
        for (s, t) in self.edges:
            S,T = self.cc_to_fc(s),self.cc_to_fc(t)
            if self.blowupfringe(S) or not self.blowupfringe(T):
                self.digraph.add_edge(S, T)
            if self.blowupfringe(T) or not self.blowupfringe(S):
                self.digraph.add_edge(T, S)

    def Rule11(self):
        """
        Removes s->t or t->s if it opposes the flow direction.
        """
        for (s, t) in self.edges:
            S = self.cc_to_fc(s)
            T = self.cc_to_fc(t)
            if self.blowupfringe(S) or self.blowupfringe(T):
                continue
            flowdata = [self.flowdirection(cctopcell, s, t) for cctopcell in self.cc.topstar(t)]
            if all(k == 1 for k in flowdata):
                self.digraph.remove_edge(T, S)
            if all(k == -1 for k in flowdata):
                self.digraph.remove_edge(S, T)

    def Rule12(self):
        """
        Removes s<->t if
            - if s is not equilibrium and
            - dim(t) > dim(s) + 1
        This is obsolete as of 01/30/2022. Rule12 is embbeded in Rule22 as a subcase of checking for dim(t)>dim(s)+1
        """
        for (s, t) in self.edges:
            if self.blowupinfinity(s) or self.blowupinfinity(t):
                continue
            S,T = self.cc_to_fc(s), self.cc_to_fc(t)
            if not self.isThereACycle(s,t):
                if abs(self.cc.cell_dim(t) - self.cc.cell_dim(s)) > 1:
                    self.digraph.remove_edge(S, T)
                    self.digraph.remove_edge(T, S)

    def Rule21(self):
        counter = 0
        for S in self.fc(self.D):
            if self.blowupfringe(S):
                continue
            s = self.fc_to_cc(S)
            if self.GradientDirections(s):
                self.digraph.remove_edge(S, S)
                counter = counter+1

    def Rule22(self):
        counter = 0
        for i in range(2,self.D+1):
            for J_e in combinations(range(self.D),self.D-i):
                coords = [ list(c) for c in product(*( range(0,n) for n in self.domains)) ]
                shape = sum([ 2**n for n in J_e ])
                for coord in coords:
                    s = self.cc.cell_index(coord,shape)
                    if self.blowupinfinity(s):
                        continue
                    S = self.cc_to_fc(s)
                    if self.Arrow(S,S): 
                        continue
                    cycles = self.CycleDecomposition(self.RegulationMap(s))
                    remove_ST,remove_TS = [],[]
                    for T in self.digraph.adjacencies(S): # S->T
                        t = self.fc_to_cc(T)
                        if s == t or s > t :
                            continue
                        if self.Arrow(T,S): # S->T exists by default
                            if not self.isThereACycle(s,t):
                                # We're shifting (s,t) through the directions
                                # From here on: s<t, s<->t and G(t)\cap J_i(t) \neq \emptyset
                                s_,t_ = s, t
                                ps_,pt_ = s,t
                                GradientInessential = [ v for v in self.GradientDirections(t) if v in self.InessentialDirections(t) ]
                                for n in GradientInessential: 
                                    if self.Psi_n(s,n) == 1:
                                        s_,t_ = self.cc.right(s_,n),self.cc.right(t_,n)
                                        ps_,pt_ = self.cc.left(ps_,n),self.cc.left(pt_,n)
                                    else: #self.Psi_n(s,n) == -1
                                        s_,t_ = self.cc.left(s_,n),self.cc.left(t_,n)
                                        ps_,pt_ = self.cc.right(ps_,n),self.cc.right(pt_,n)
                                S_,T_ = self.cc_to_fc(s_),self.cc_to_fc(t_)
                                PS_,PT_ = self.cc_to_fc(ps_),self.cc_to_fc(pt_)
                                if self.Arrow(PT_,PS_) and self.Arrow(T_,S_):
                                    remove_ST.append(T)
                                elif self.Arrow(PS_,PT_) and self.Arrow(S_,T_):
                                    remove_TS.append(T)
                            # The following is from the original code and was adapted to not run Case III of Rule 2 in N=3
                            # # otherwise it is different and we want it opposing s_ -> t_
                            #     elif self.Arrow(T_,S_):
                            #         remove_TS.append(T)
                            #     elif self.Arrow(S_,T_):
                            #         remove_ST.append(T)
                                # counter = counter+1
                                elif self.Arrow(T_,S_) and i == 2:
                                    remove_TS.append(T)
                                elif self.Arrow(S_,T_) and i == 2:
                                    remove_ST.append(T)
                                counter = counter+1

                    for T in remove_ST:
                        self.digraph.remove_edge(S,T)
                    for T in remove_TS:
                        self.digraph.remove_edge(T,S)

    def Rule31(self):
        counter = 0
        for (s,t) in self.edges:
            # Check if cell is a boundary cell. If so, continue.
            if s == t:
                continue
            if self.blowupinfinity(s) or self.blowupinfinity(t):
                continue

            # Check for double-arrow
            S,T = self.cc_to_fc(s), self.cc_to_fc(t)
            if not self.DoubleArrow(S,T):
                continue

            # Decompose the local inducement into its cycle-decomposition
            cycles = self.CycleDecomposition(self.RegulationMap(s))
            cyclic = [ n for cycle in cycles for n in cycle ]

            # Check if cell is an equilibrium cell
            if not cycles:
                continue
            if set(cyclic) != set(self.InessentialDirections(s)):
                continue
         
            if self.isThereACycle(s,t):
                # Find the relevant cycle that contain Ext(\xi,\xi'),
                cycles = [ cycle for cycle in cycles if set(self.JiCapJe(s,t)).issubset(set(cycle))]

                # If there are no such cycles, it means that Ext(\xi,\xi') is contained in two different cycles. Remove both S->T and S<-T in that case.
                if cycles == []: 
                    self.digraph.remove_edge(S,T)
                    self.digraph.remove_edge(T,S)

                # Verify which top cells are unstable using the lap-number
                unstable = {}
                for kappa in self.cc.topstar(s):
                    unstable_kappa = { str(cycle) : (self.CrossingNumber(self.RelativePositionVector(s,kappa),cycle) < (len(cycle)/2)) for cycle in cycles }
                    unstable.update({ kappa : unstable_kappa })

                # Remove S->T if T is unstable, i.e., T is in the "intersection" of unstable cells
                removed = False
                for kappa in self.cc.topstar(t):
                    for k,cycle in enumerate(cycles):
                        if any(n in cycle for n in self.JiCapJe(s,t)):
                            if not unstable[kappa][str(cycle)]:
                                self.digraph.remove_edge(S,T)
                                counter = counter+1
                                removed = True
                # Otherwise, it has a stable top-dim cell nearby. So it is stable for our purposes. Remove S<-T. 
                if not removed:
                    self.digraph.remove_edge(T,S)
                    counter = counter+1

    def checkfordoublearrow(self):
        num_double_edges = 0
        # num_double_edges_in_equilibrium = 0
        for v in self.digraph.vertices():
            for u in self.digraph.adjacencies(v):
                if self.blowupfringe(u):
                    continue
                if v in self.digraph.adjacencies(u) and v < u:
                # if v in self.digraph.adjacencies(u) and v != u:
                    # print("from v",self.fc_to_cc(v))
                    # print("to u",self.fc_to_cc(u))
                    num_double_edges += 1
                    # if self.isEquilibriumCell(self.fc_to_cc(v)) or self.isEquilibriumCell(self.fc_to_cc(u)):
                    #     num_double_edges_in_equilibrium += 1
        # print('Number of double edges:', num_double_edges)
        return num_double_edges

    def __init__(self, param, pg, level=0):
        net_spec = param.network().specification()
        par_index = pg.index(param)
        network = DSGRN.Network(net_spec, edge_blowup='none')
        parameter_graph = DSGRN.ParameterGraph(network)
        p = parameter_graph.parameter(par_index)
        """ initialize complex, graph, supporting data structures """
        # Pretty sure you don't need all the below since it was from old code
        self.parameter = p  # parameter
        self.network = self.parameter.network()
        self.D = self.network.size()
        self.domains = [len(self.network.outputs(d)) + 1 for d in range(self.D)]
        # self.domains = self.network.domains()  # a list of numbers of thresholds
        self.cc = CubicalComplex([x + 1 for x in self.domains])  # original cubical complex
        self.digraph = DirectedAcyclicGraph()
        self.labelling = self.parameter.labelling()
        self.limits = [len(self.network.outputs(d))+ 1 for d in range(self.D)]

        self.pv = [1]
        for i in self.limits:
            self.pv.append(self.pv[-1] * i)

        self.occ = OrderedCubicalComplex(self.cc)  # ordered cubical complex
        self.edges = self.occ.Edges()  # set of (s, t) where s<t

        self.fc = CubicalComplex([2 * x + 2 for x in self.domains])  # blowup/face cubical complex
        for cell in self.fc(self.D):  # combinatorial dynamical system
            self.digraph.add_vertex(cell)

        self.D_cc_to_fc, self.D_fc_to_cc = self.Dict_cc_fc()

        self.w = {}

        # STEP 1: Build local inducement map
        self.rmap = self.build_rmap()
        # rmap is a list of vectors given by the relation
        # rmap[i][j] = T[i -> rmap[i][j]]

        # STEP 2: Build complex
        start = time.time()
        self.Rule0()
        end = time.time()

        # STEP 3: Apply rules
        # Define F_1
        if level >= 1:
            start = time.time()
            self.Rule11() # flow data
            end = time.time()
            start = time.time()
            self.Rule12()
            end = time.time()
            start = time.time()
            self.Rule21() # remove self edge
            end = time.time()
        # Define F_2
        if level >= 2:
            start = time.time()
            self.Rule22()
            end = time.time()
        # Define F_3
        if level >= 3:
            start = time.time()
            self.Rule31()
            end = time.time()
