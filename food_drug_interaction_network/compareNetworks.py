import networkx as nx

def collapseNodes(G, mapping):

    """
    given a graph G, and a mapping of nodes G to other nodes
    make a new network G2 where the neighbors of G2 are the union
    of the neighbors of all nodes in G that map to that node in G2
    """
    G2 = nx.Graph()

    for u,v in G.edges():

        if u in mapping and v in mapping:

            u_map = mapping[u]
            v_map = mapping[v]

            if u_map != v_map:

                G2.add_edge(u_map, v_map)

    return G2

def removeInterSetEdges(G,seta,setb):

    """
    
    From DrugBank network, remove edges from compounds in the same set
    as in the bipartite network, unless they are in both sets
    
    """
    to_remove = []
    for u,v in G.edges():

        u_in_a = u in seta
        u_in_b = u in setb
        v_in_a = v in seta
        v_in_b = v in setb
        
        if (u_in_a and u_in_b) or (v_in_a and v_in_b):

            pass

        elif (u_in_a and v_in_a) or (u_in_b and v_in_b):

            to_remove.append((u,v))

    G.remove_edges_from(to_remove)

def removeNodes(G1,G2):

    """
    Remove nodes not in both networks
    """

    to_remove_1 = []
    to_remove_2 = []

    for node in G1.nodes():

        if node not in G2:

            to_remove_1.append(node)

    for node in G2.nodes():

        if node not in G1:

            to_remove_2.append(node)

    G1.remove_nodes_from(to_remove_1)
    G2.remove_nodes_from(to_remove_2)

def compareEdges(G1,G2):

    """
    Compute intersection, both differences of 
    edge sets
    """

    in_G1_only = 0
    in_G2_only = 0
    in_both = 0

    for u,v in G1.edges():

        if G2.has_edge(u,v):

            in_both += 1

        else:

            in_G1_only += 1

    for u,v in G2.edges():

        if not G1.has_edge(u,v):

            in_G2_only += 1

    return in_both, in_G1_only, in_G2_only

