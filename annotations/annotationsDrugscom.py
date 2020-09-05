import networkx as nx
import requests
import re

def _get_all_parents(G, n, res = None):

    """
    Get parents of a node from a DAG
    Important, if it isnt a DAG it may get stuck on an infinite loop in a cycle!
    """

    if res is None:

        res = []

    for parent in G.predecessors(n):

        res.append(parent)
        _get_all_parents(G, parent, res)

    return res


def get_all_parents(G,n):
    """
    Checks for DAG at beginning
    """
    if not nx.is_directed_acyclic_graph(G):

        raise TypeError("Graph MUST be direcetd and acyclic")


    return _get_all_parents(G,n)

def classifyCp(cp):

    regexp = r"Drug class: ([^<]*)</a>"
    searchurl = "https://www.drugs.com/mtm/"+ cp +".html"
    response = requests.get(searchurl)
    response.raise_for_status()

    return re.findall(regexp,response.text.replace("\n", ""))
