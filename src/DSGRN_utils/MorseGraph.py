### MorseGraph.py
### MIT LICENSE 2023 Marcio Gameiro

import pychomp

def MorseGraph(stg, scc_dag, graded_complex, connection_matrix):
    """Construct the Morse graph"""

    def nontrivial_scc(v):
        """Check if the SCC with grading v is nontrivial"""
        if v == fringe_node:
            return False
        if v in conley_indices:
            return True
        # Get list of cells in the SCC with grading v
        scc_v = [c for c in stg.digraph.vertices() if graded_complex.value(c) == v]
        # Check if there are common gradient directions
        common_grad_dirs = set()
        for c in scc_v:
            # Get gradient directions for cell c
            grad_dirs = stg.GradientDirections(stg.fc_to_cc(c))
            # Initialize common_grad_dirs
            if not common_grad_dirs:
                common_grad_dirs.update(grad_dirs)
            # Intersect with previous gradient directions
            common_grad_dirs.intersection_update(grad_dirs)
            # Nontrivial if no common gradient directions
            if not common_grad_dirs:
                return True
        return False

    def vertex_label(v):
        """Vertex label for Morse graph"""
        return str(verts_index[v]) + " : " + str(tuple(conley_indices[v]))

    # Get the grading for fringe cells
    fringe_node = graded_complex.value(stg.complex().size() - 1)
    # Get the cell complex dimension
    D = graded_complex.complex().dimension()
    # Get the nontrivial Conley indices
    conley_indices = connection_matrix.count()
    # Get the Morse graph vertices
    morse_graph_verts = set([v for v in scc_dag.vertices() if nontrivial_scc(v)])
    # Make an indexing of the Morse graph vertices
    verts_index = {v: k for k, v in enumerate(sorted(morse_graph_verts))}
    # Update conley_indices with the trivial Conley indices
    for v in morse_graph_verts:
        if v not in conley_indices:
            conley_indices[v] = [0] * (D + 1)
    # Define the Morse graph
    morse_graph = pychomp.DirectedAcyclicGraph()
    # Add vertices and edges
    for v in morse_graph_verts:
        morse_graph.add_vertex(v, label=vertex_label(v))
    for v in morse_graph_verts:
        for u in scc_dag.descendants(v):
            if (u == v) or (u not in morse_graph_verts):
                continue
            morse_graph.add_edge(v, u)
    # Morse graph is the transitive reduction of morse_graph
    return morse_graph.transitive_reduction()
