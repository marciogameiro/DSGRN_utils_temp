# SetCM.py
# MIT LICENSE 2020 Ewerton R. Vieira & Marcio Gameiro

import numpy as np

class SetCM:
    """
    Overview:
        Compute the set of Connection Matrices
    Inputs:
        Connection Matrix (pyCHomP)
        Poset induced by Connection Matrix
        Fringenode = fibration.value(std.complex().size()-1) tem que retirar isso
    """

    def CM_to_Dict(self):
        """
        For CM, return the associated Dictonary
        """
        Delta_0 = {}
        for d in range(self.CM.complex().dimension() + 1):
            for c in self.CM.complex()(d):  # for each cell c in CM
                # cell in the original chain complex that corresponds to cell c in CM
                cell = self.CM.value(c)
                # all cell in the boundary of c
                if cell != self.fringenode:
                    bd_cell = {b for b in self.CM.complex().boundary({c})}
                    Delta_0.update({c: bd_cell})
        return Delta_0

    def CMcell_to_Orig(self):
        j = 0
        J = {}  # dict of number of cells in same morse node
        M = {}  # dict CM cells to original cells with index in the case of same morse node
        for c in range(self.CM.complex().size()):
            cell = self.CM.value(c)
            if cell != self.fringenode:
                if cell in J.keys():
                    j += J.get(cell) + 1
                    J.update({cell: j})
                else:
                    J.update({cell: 0})
                M.update({c: (cell, J.get(cell))})
        return M

    def CM_to_origComp(self, D):
        """
        For D (in dictonary form) return new D_ with cell that correspond to the original complex
        """
        D_ = {}  # new CM in dict form
        keys = enumerate(D.keys())
        M = self.CMcell_to_Orig()

        for _, d in keys:
            cell = self.CM.value(d)
            if cell != self.fringenode:
                bd_cell = {M.get(b) for b in D.get(d)}
                D_.update({M.get(d): bd_cell})

        return D_

    def Elementary_TM(self):
        """
        For CM return the set of Elementary Transition Matrices(E)
        each E is viewed as an tuple (a,b,c) where a is the dimension of the 0 degree map from c to b
        """
        E_list = set()
        C_index = self.CM.count()  # List of Conley indices

        for i, _ in enumerate(self.CM.complex()):
            for j, _ in enumerate(self.CM.complex()):
                if i != j:
                    I = self.CM.value(i)
                    J = self.CM.value(j)
                    I_index = C_index.get(I)  # Conley index of I

                    if I != self.fringenode:
                        if I == J:  # true means that there is more then one cell in CM.complex() for I
                            for k, a in enumerate(I_index):
                                if a > 1:
                                    # Update the automorphism j to i
                                    E_list.update({(k, i, j)})
                        if J != self.fringenode:

                            if self.CMG.less(J, I):  # true if I<J

                                J_index = C_index.get(J)
                                # k is the dimension of the cell a
                                for k, a in enumerate(I_index):

                                    # Non trivial cell in dimension k
                                    if a * J_index[k] != 0:
                                        # Elementary transition matrix with off diagonal (i,j)
                                        E_list.update({(k, i, j)})
        return E_list

    # Given D return EDE^-1 where E is elementary transition matrix
    def conjugation(self, a, b, c, D):
        D.update({c: D.get(b) ^ D.get(c)})
        for _, i in enumerate(self.CM.complex()(a + 1)):
            if self.CM.value(i) == self.fringenode:
                continue
            Di = D.get(i)

            if {b, c}.issubset(Di):
                D.update({i: Di - {b}})
            elif Di - {b, c} == Di:
                continue
            else:
                D.update({i: Di | {b}})

        return D

    def Compute_all_CM(self, D, setCM=None):
        if not setCM:
            setCM = []
        setCM.append(str(D))
        for a, b, c in self.Elementary_TM():
            D_p = self.conjugation(a, b, c, D)
            if str(D) not in setCM:
                self.Compute_all_CM(D_p, setCM)
        return setCM

    def Get_CM(self):
        return self.All_CM

    def matrix_form(self, D):
        k0 = len([i for i in self.CM.complex()(0)])
        k1 = len([i for i in self.CM.complex()(self.CM.complex().dimension())])
        lj = len(D.keys())  # colum size without 0 dimension cells
        li = len(D.keys()) - k1  # row size without topcells
        A = np.zeros((li, lj))  # smaller matrix version of D
        for i, a in enumerate(D.keys()):
            for j, b in enumerate(D.keys()):
                if all([i <= li, j >= k0, a in D.get(b)]):
                    A[i, j] = 1
        return A

    def exclude_cm(self):
        """
        Comment: Band-aid to exclude matrices that have \del(a) jumping 2 or more dimenions
        Check each M in self.All_CM and delete it whenever there is an element c_n such that d(c_n) = c_i, i \neq n-1
        """
        # The dimension of the highest cell is (including the fringenode)
        d = self.CM.complex().dimension()
        # List with amount of sets in each level set
        count = self.CM.complex().count()
        rcount = [ i for i in reversed(count)]
        # List of sets in each level_set[i]
        level_set = [ [i for i in self.CM.complex()(n)] for n in range(d) ]

        remove = []
        for k,D in enumerate(self.All_CM):
            E = eval(D)
            M = self.matrix_form(E)
            for i in range(1,d-1):
                beg1, end1 = 0, sum(count[:i])
                beg2 = sum(count)-sum(rcount[:(i+1)])
                end2 = sum(count)-sum(rcount[:i])
                if np.any(M[beg1:end1,beg2:end2]):
                    remove.append(k)
        # print("list of matrices to remove:",remove)
        for i in reversed(remove):
            del self.All_CM[i]
        # Now you might have checked the boundary but there might be 2 identical matrices (because of reasons), so we check for repetited matrices
        remove = []
        n = len(self.All_CM)
        for i in range(n):
            for j in range(i+1,n):
                if np.allclose(self.matrix_form(eval(self.All_CM[i])),self.matrix_form(eval(self.All_CM[j]))):
                    remove.append(j)
        # Need to make sure remove is sorted (Marcio)
        for i in reversed(sorted(remove)): # remove in reversed order to preserve the index
            del self.All_CM[i]

        # # Failed attempt at exclusion:
        # for D in self.All_CM:
        #     remove = False
        #     M = eval(D)
        #     for i in range(1,d):
        #         for cell in level_set[i]:
        #             for boundary in M[cell]:
        #                 if boundary not in level_set[i-1]:
        #                     remove = True
        #     if remove:
        #         self.All_CM.remove(D)

    def __init__(self, CM, CMG, fringenode):
        self.CM = CM  # Connection Matrix
        self.CMG = CMG  # Poset
        # fringe cell that the blowup face complex has (DSGRN+)
        self.fringenode = fringenode
        self.All_CM = self.Compute_all_CM(self.CM_to_Dict())

        # Bernardo's bandaid
        self.exclude_cm()
