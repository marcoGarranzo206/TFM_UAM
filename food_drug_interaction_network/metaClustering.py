import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.spatial import distance as ssd
import sys
sys.setrecursionlimit(90000)

class DisjointSet:

    """
    Barebones implementation of disjoint set
    The one in pypi ran into an unresolved error for high number of union operations
    Each array index represents an element, and the value, its parent
    2 nodes are in the same set if they share the same root (node with itself as a parent)
    
    union: join 2 sets
    find: find the set (root) of node i. traverse array until root. Along the way
    change all nodes parents to the root for faster finds in future

    """

    def __init__(self,n = None):

        
        self.size = n
        self.parent = [i for i in range(n)]
        self.rank = [0 for _ in range(n)]

    def find(self,i):

        if i > self.size - 1:

            raise ValueError(f"{i} exceeds set length")

        if i != self.parent[i]:

            self.parent[i] = self.find(self.parent[i])

        return self.parent[i]

    def union(self, i,j):

        i_id = self.find(i)
        j_id = self.find(j)
        if i_id == j_id:

            return
        if self.rank[i_id] > self.rank[j_id]:

            self.parent[j_id] = i_id

        else:

            self.parent[i_id] = j_id
            if self.rank[i_id] == self.rank[j_id]:

                self.rank[j_id] = self.rank[j_id] + 1

class metaCluster:


    def __init__(self, cluster_runs, labels):

        """
        labels: names for the points
        cluster_runs. Matrix of (n_cluster_runs, n_points)
        where the i,j element is the community assigned at point j
        in cluster run i
        """
        self.labels = labels
        self.coMatrix = self._getCoMatrix(cluster_runs)

    def majorityVote(self,t):

        clusters = DisjointSet(len(self.labels))

        for i in range(self.coMatrix.shape[0]):

            for j in range(self.coMatrix.shape[1]):

                if self.coMatrix[i,j] > t:

                    clusters.union(i,j)
        return np.array([clusters.find(i) for i in range(len(self.labels))])

    def _getCoMatrix(self,cluster_runs):

        coMatrix = np.zeros((cluster_runs.shape[1], cluster_runs.shape[1]))
        for row in cluster_runs:

            unique_vals = np.unique(row)
            for val in unique_vals:

                idxs = np.where(row == val)[0]
                meshed = np.meshgrid(idxs, idxs)
                coMatrix[meshed[0].reshape(1,-1), meshed[1].reshape(1,-1)] += 1
                

        coMatrix /= cluster_runs.shape[0]
        return coMatrix

    def HAC(self, t, distance):

        D = 1 - self.coMatrix
        linkage = sch.linkage(ssd.squareform(D), method = distance, metric = "euclidian")
        return sch.fcluster(linkage, t, "distance")
