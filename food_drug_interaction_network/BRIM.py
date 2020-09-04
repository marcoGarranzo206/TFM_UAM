import networkx as nx
import numpy as np

def Q(R_T, B, S,m):
    
    return np.trace((R_T @ B @ S))/m

def BRIM_loop(B,c, S, m):

    """
    Given modularity matrix B (or a rank approximation of it) and a number of modules c,
    run an iteration of BRIM algorithm described in 
    
    There are two sets (bipartite graph)
    Assign modules to the first set (say A) randomly.
    Based on those assingments, assign modules to the other set (say B)
    Based on this new B assignment, reassign the first set to new modules
    
    Do until convergence: NOT IMPLEMENTED
    
    """
    
    Q_old = -1
    while True:
        
        #assign to R_T based on S
        T = B @ S
        R_T = np.zeros((c, B.shape[0]))
        R_modules = np.argmax(T, axis=1).reshape(1,-1)
        R_T[R_modules , range(B.shape[0])] = 1
        #assign to S based on R_T
        
        T = R_T @ B 
        S = np.zeros((B.shape[1], c))
        S_modules = np.argmax(T, axis=0)
        S[range(B.shape[1]), S_modules ] = 1
        Q_new = Q(R_T,B,S,m)

        if Q_new > Q_old:

            Q_old = Q_new

        else:

            break

    return R_T, B, S, Q_new

def BNullBipartiteConfigPos(g, resolution = 1):
        
    """
    Computes B_tilde of a bipartite graph following a null configuration model
    Prob of edge (i,j) is 0 if they belong to the same set, degree_i*degree_j/m otherwise
    will end up with a square matrix of two diagonal 0 blocks, and 2 non diagonal blocks
    who are each others transpose
    
    As such, save space computing only one of the non 0 blocks
    """
    seta, setb = nx.bipartite.sets(g)

    if len(seta) > len(setb):

        s1 = seta
        s2 = setb

    else:

        s1 = setb
        s2 = seta
        
    
    A = nx.bipartite.biadjacency_matrix(g,row_order= s2, column_order=s1)
    ki = np.sum(A, axis=1)
    dj = np.sum(A, axis = 0)
    m = np.sum(A)
    B = A - resolution*((ki@dj)/m)
        
    return B, m
        

def BNullBipartiteConfigNeg(g, resolution = 1):
    
    """
    Computes B_tilde of a bipartite graph following a null configuration model that can include
    both pos and neg edges as described by gomez et al in 
    Analysis of community structure in networks of correlated data (2009)
    
    Prob of any edge (i,j) is 0 if they belong to the same set
    Prob of pos edge (i,j) is pos_degree_i*pos_degree_j/w_pos
    Prob of neg edge (i,j) is neg_degree_i*neg_degree_j/w_pos
    
    B = A - prob_pos + prob_neg
    
    will end up with a square matrix of two diagonal 0 blocks, and 2 non diagonal blocks
    who are each others transpose
    
    As such, save space computing only one of the non 0 blocks
    """
    
    seta, setb = nx.bipartite.sets(g)
    
    
    if len(seta) > len(setb):

        s1 = seta
        s2 = setb

    else:

        s1 = setb
        s2 = seta
 
    A = nx.bipartite.biadjacency_matrix(g,row_order= s2, column_order=s1).toarray()
    A_pos = np.zeros(A.shape)
    A_pos[A > 0] = 1
    
    A_neg = np.zeros(A.shape)
    A_neg[A < 0] = 1 
   
    ki_pos = np.sum(A_pos, axis=0).reshape(-1,1)
    dj_pos = np.sum(A_pos, axis = 1).reshape(1,-1)
    
    ki_neg = np.sum(A_neg, axis=0).reshape(-1,1)
    dj_neg = np.sum(A_neg, axis = 1).reshape(1,-1)
    
    
    w_pos = np.sum(A_pos)
    w_neg = np.sum(A_neg)
    m = w_pos + w_neg
    
    B = A - resolution*((ki_pos@dj_pos).T/w_pos) + (1/resolution)*((ki_neg@dj_neg).T/w_neg)
    return B, m, w_pos, w_neg
    

def BRIM(g,c,resolution = 1,null = "config", seed = None, assingments = None):
    """
    Find modules using brim alg for a fixed number of c (up to c)
    """
    if null in ("config", "configuration", "configPos", "configuartionPos"):
    
        B, m = BNullBipartiteConfigPos(g, resolution)
        
    elif null in ("configNeg", "configurationNeg"):
        
        B, m, wpos, wneg = BNullBipartiteConfigNeg(g, resolution)
    
    return _BRIM(B,c,m,seed = seed, assingments = assingments)
    
def _BRIM(B, c, m,seed = None, assingments = None):


    S = np.zeros((B.shape[1], c))
    if assingments is None:

        if seed is not None:
            
            np.random.seed(seed)

        S[range(B.shape[1]), np.random.choice(c, B.shape[1])] = 1

    else:

        S[range(B.shape[1]), assingments] = 1

    return BRIM_loop(B,c,S,m)

def BRIM_bisec(g,resolution = 1, c_max = None,null = "config", B = None, m = None):
    
    """
    Find modules using brim alg without a fixed number of c, using binary search method described
    in BRIM paper (not exactly for now)
    
    start with low c, calc Q
    double c until Q decreases
    go to middle point between new c and olde c
    repeat until convergence
    """
    
    if B is None: #for repeated runs avoids having to compute B

        if null in ("config", "configuration", "configPos", "configuartionPos"):
        
            B, m = BNullBipartiteConfigPos(g, resolution)
            
        elif null in ("configNeg", "configurationNeg"):
            
            B, m, wpos, wneg = BNullBipartiteConfigNeg(g, resolution)
        
    if c_max is None:
        
        c_max = min(B.shape) 
    
   
    c = 2
    half = B.shape[1]//2
    assingments_s = np.ones(B.shape[1]).astype(int)
    Q_old = -100
    prev_c = 1
    c_max = B.shape[1]
    while True:


        assingments_s[np.random.choice(range(B.shape[1]), half, replace = False )] = np.random.randint(low = 0, high = c, size = half)
        R_T, B, S, Q_new = _BRIM(B,c,m,assingments = assingments_s)
        assingments_s = np.where(S != 0 )[1]
        
        if Q_new > Q_old:

            R_T_max, B_max, S_max, Q_max = R_T, B, S, Q_new
            if c == B.shape[1]:

                return R_T_max, B_max, S_max, Q_max
        
            prev_c = c
            c = min(2*c, B.shape[1] - 1)
            Q_old = Q_new

        else:

            break

    high = c
    low = prev_c
    c = (high + low)//2
    c_prev = c
    Q_upper = Q_new
    Q_lower = Q_old

    while True:

        if high < low or c_prev == c:

            return R_T_max, B_max,S_max, Q_max
        
        assingments_s[np.random.choice(range(B.shape[1]), half, replace = False )] = np.random.randint(low = 0, high = c, size = half)
        R_T, B, S, Q_new = _BRIM(B,c,m,assingments = assingments_s)
        assingments_s = np.where(S != 0 )[1]
        
        if Q_new > Q_max:

            R_T_max, B_max, S_max, Q_max = R_T, B, S, Q_new

        if  Q_new > Q_upper and Q_lower > Q_upper:

            #possible Q local maximum by searching between
            #low and high c
            #is a heuristic, even though Q_new > Q_upper, there could be
            # a maxima in between
            Q_upper = Q_new
            high = c
            c = (high + low)//2

        else:

            # euther possible maxima in the other end using the above heuristic
            #or stuck in a minimum or maximim between both c, in which case arbitrarily
            #go the other end too

            Q_lower = Q_upper
            low = c
            c = (high + low)//2

        c_prev = c
        
